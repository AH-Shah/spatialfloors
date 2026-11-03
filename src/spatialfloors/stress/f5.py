"""F5 stress panel (PLAN §13) as a reusable context: built once per (unit, perturbation), profiles any fine-label call.
  agreement   held-out-marker AUROC (primary), per-class AUROC and the worst class, IF-gated macro-F1 (GSE1)
  references  marker threshold: per class, call the n_pred cells with the highest score on the class's top-5 reference
              INPUT genes (how easy the metric is at that call rate)
  degeneracy  predicted share per class and "other" (collapse / dumping)
  confounds   depth-stratified AUROC (mean over library-size quintiles), share of called cells in the top depth quintile"""
import json
from pathlib import Path
import anndata as ad, numpy as np, scipy.sparse as sp
from sklearn.metrics import roc_auc_score
from spatialfloors.data.genes import filter_reason
from spatialfloors.floors import f5
from spatialfloors.metrics.f5 import if_gate2, macro_f1_gated, heldout_auroc, CLASSES, MIN_PRED
from spatialfloors.stress import perturb

MK0 = json.loads(Path("manifests/panels/index.json").read_text())["f5_heldout_markers"]; MK = MK0


def load(u):
    """Unit inputs: all counts (csr), gene names, input genes (panel minus held-out markers, in Pelka), NegPrb means, IF gate."""
    global MK; MK = f5.for_unit(u, MK0)  # new cancers: own reference and markers (PLAN §17); CRC unchanged
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); p = json.loads(Path(f"manifests/panels/{u}.json").read_text())
    held = {x for v in MK.values() for x in v}; pg = set(f5.pelka_genes()); genes = [x for x in p["genes"] if x not in held and x in pg]
    negc = [x for x in a.var_names if x.startswith(("NegPrb", "Negative"))]
    gate = if_gate2(a.obs.IF_PanCK.values, a.obs.IF_CD68_CK8_18.values)[0] if u.startswith("gse1_") else None
    return dict(unit=u, X=sp.csr_matrix(a.X), var=np.asarray(a.var_names).astype(str), genes=genes, fov=a.obs.fov.astype(str).values,
                neg=np.asarray(a[:, negc].X.mean(1)).ravel() if negc else None, gate=gate)


def inputs(d, name="none", seed=0):
    """Perturbed unit: (Xq input counts, input gene list, Xf scoring counts, scoring gene names, neg, gate, cells)."""
    X, c, _ = perturb.apply(name, d["X"], d["fov"], seed) if name.startswith(("thin", "sub", "fov")) else (d["X"], np.arange(d["X"].shape[0]), None)
    genes = list(d["genes"]) if not name.startswith("genes") else [d["genes"][i] for i in perturb.subsample(len(d["genes"]), int(name[5:]) / 100, seed)]
    vi = {g: i for i, g in enumerate(d["var"])}; Xq = sp.csr_matrix(X[:, [vi[g] for g in genes]])
    keep = filter_reason(d["var"], "human") == ""; Xf, vf = sp.csr_matrix(X[:, keep]), list(d["var"][keep])
    negc = [i for i, g in enumerate(d["var"]) if g.startswith(("NegPrb", "Negative"))]; neg = np.asarray(X[:, negc].mean(1)).ravel() if negc else None
    return Xq, genes, Xf, vf, neg, None if d["gate"] is None else d["gate"][c], c


def auroc(pos, sc):
    return float(roc_auc_score(pos, sc)) if MIN_PRED <= pos.sum() <= len(pos) - MIN_PRED else .5


class Panel:
    def __init__(self, Xq, genes, Xf, vf, gate, prof, names):
        """prof, names: reference mean profiles (genes x fine labels) on `genes`, for the marker-threshold reference."""
        self.Xf, self.vf, self.gate = Xf, vf, gate; lib = np.maximum(np.asarray(Xf.sum(1)).ravel(), 1)
        self.dq = np.digitize(np.log(lib), np.quantile(np.log(lib), [.2, .4, .6, .8]))
        self.hscore = {c: np.log1p(np.asarray(Xf[:, [vf.index(g) for g in MK[c] if g in vf]].sum(1)).ravel() / lib * 1e4) for c in CLASSES}
        P = np.asarray(prof); names = list(names); ilib = np.maximum(np.asarray(Xq.sum(1)).ravel(), 1); self.iscore = {}
        for c in CLASSES:  # top-5 input genes by log-fold of the class's mean profile over the other classes' profiles
            inc = [i for i, n in enumerate(names) if f5.COARSE.get(n) == c]; out = [i for i in range(len(names)) if i not in inc]
            top = np.argsort(np.log((P[:, inc].mean(1) + 1e-3) / (P[:, out].mean(1) + 1e-3)))[::-1][:5]
            self.iscore[c] = np.log1p(np.asarray(Xq[:, top].sum(1)).ravel() / ilib * 1e4)

    def __call__(self, fine):
        c = f5.coarse(fine); h, per = heldout_auroc(c, self.Xf, self.vf, MK); r = dict(heldout_auroc=h, worst_class=float(min(v["auroc"] for v in per.values())))
        r.update({f"share_{k}": float((c == k).mean()) for k in (*CLASSES, "other")})
        ex, ds, td = [], [], []
        for k in CLASSES:
            pos = c == k; n = int(pos.sum()); thr = np.zeros(len(c), bool); thr[np.argsort(self.iscore[k])[::-1][:n]] = True
            r[f"auroc_{k}"] = per[k]["auroc"]; ex.append(per[k]["auroc"] - auroc(thr, self.hscore[k]))
            ds.append(float(np.mean([auroc(pos[self.dq == q], self.hscore[k][self.dq == q]) for q in range(5)]))); td.append(float((self.dq[pos] == 4).mean()) if n else np.nan)
        r.update(excess_over_threshold=float(np.mean(ex)), depth_strat_mean=float(np.mean(ds)), top_depth_share=float(np.nanmean(td)))
        if self.gate is not None: r["gated_macro_f1"] = macro_f1_gated(self.gate, c, ("epithelial", "myeloid"))[0]
        return r
