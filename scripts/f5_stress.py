"""F5 development stress panel (PLAN §13): does a transfer win the held-out-marker AUROC by typing cells, or by calling
only the easiest / deepest cells? Development crops only; predictions from runs/p3*/labels/f5_<unit>_*.npz.
  agreement   held-out-marker AUROC (primary), per-class AUROC and the worst class, IF-gated macro-F1 (GSE1)
  references  marker threshold: per class, call the n_pred cells with the highest score on the class's top-5 reference
              (Pelka profile) INPUT genes (no model; how easy the metric is at that call rate); random calls at the
              same shares (~.5)
  degeneracy  predicted share per class and "other" (collapse / dumping); call-rate-matched AUROC (each class scored on
              its top-n cells by the model's calls vs the marker-threshold reference at the same n)
  confounds   depth-stratified AUROC (mean over library-size quintiles: a model calling deep cells wins on CP10k markers),
              share of called cells in the top depth quintile
Usage (env floors): python scripts/f5_stress.py <unit> -> runs/p3/f5_stress_<unit>.json"""
import glob, json, sys
from pathlib import Path
import anndata as ad, numpy as np, scipy.sparse as sp
from sklearn.metrics import roc_auc_score
from spatialfloors.data.genes import filter_reason
from spatialfloors.floors import f5
from spatialfloors.metrics.f5 import if_gate2, macro_f1_gated, heldout_auroc, CLASSES, MIN_PRED
from spatialfloors.provenance import stamp

u = sys.argv[1]; a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); p = json.loads(Path(f"manifests/panels/{u}.json").read_text())
MK = f5.for_unit(u, json.loads(Path("manifests/panels/index.json").read_text())["f5_heldout_markers"]); held = {x for v in MK.values() for x in v}
genes = [x for x in p["genes"] if x not in held and x in set(f5.pelka_genes())]
keep = filter_reason(a.var_names, "human") == ""; Xf, vf = sp.csr_matrix(a.X[:, keep]), list(map(str, a.var_names[keep])); lib = np.maximum(np.asarray(Xf.sum(1)).ravel(), 1)
dq = np.digitize(np.log(lib), np.quantile(np.log(lib), [.2, .4, .6, .8])); log = lambda *s: print(*s, flush=True)
gate = if_gate2(a.obs.IF_PanCK.values, a.obs.IF_CD68_CK8_18.values)[0] if u.startswith("gse1_") else None
hscore = {c: np.log1p(np.asarray(Xf[:, [vf.index(g) for g in MK[c] if g in vf]].sum(1)).ravel() / lib * 1e4) for c in CLASSES}
prof, names = f5.reference(genes, rows="all"); names = list(names); P = np.asarray(prof)  # genes x fine labels
Xi = sp.csr_matrix(a[:, genes].X); ilib = np.maximum(np.asarray(Xi.sum(1)).ravel(), 1)
iscore = {}
for c in CLASSES:  # top-5 input genes by log-fold of the class's mean profile over the other classes' profiles
    inc = [i for i, n in enumerate(names) if f5.COARSE.get(n) == c]; out = [i for i, n in enumerate(names) if i not in inc]
    lfc = np.log((P[:, inc].mean(1) + 1e-3) / (P[:, out].mean(1) + 1e-3)); top = np.argsort(lfc)[::-1][:5]
    iscore[c] = np.log1p(np.asarray(Xi[:, top].sum(1)).ravel() / ilib * 1e4)


def auroc(pos, sc):
    return float(roc_auc_score(pos, sc)) if MIN_PRED <= pos.sum() <= len(pos) - MIN_PRED else .5


def panel(fine):
    c = f5.coarse(fine); h, per = heldout_auroc(c, Xf, vf, MK); r = dict(heldout_auroc=h, worst_class=float(min(v["auroc"] for v in per.values())))
    r["share"] = {k: float((c == k).mean()) for k in (*CLASSES, "other")}
    for k in CLASSES:
        pos = c == k; n = int(pos.sum()); thr = np.zeros(len(c), bool); thr[np.argsort(iscore[k])[::-1][:n]] = True
        r[k] = dict(auroc=per[k]["auroc"], n_pred=n, marker_threshold_auroc=auroc(thr, hscore[k]),
                    depth_strat_auroc=float(np.mean([auroc(pos[dq == q], hscore[k][dq == q]) for q in range(5)])),
                    top_depth_share=float((dq[pos] == 4).mean()) if n else float("nan"))
        r[k]["excess_over_threshold"] = r[k]["auroc"] - r[k]["marker_threshold_auroc"]
    r["depth_strat_mean"] = float(np.mean([r[k]["depth_strat_auroc"] for k in CLASSES])); r["threshold_mean"] = float(np.mean([r[k]["marker_threshold_auroc"] for k in CLASSES]))
    if gate is not None: r["gated_macro_f1"] = macro_f1_gated(gate, c, ("epithelial", "myeloid"))[0]
    return r


rep = dict(unit=u, predictions={})
for f in sorted(set(glob.glob(f"runs/p3*/labels/f5_{u}_*.npz") + glob.glob(f"runs/p3/relabel/labels/f5_{u}_*.npz"))):
    z = np.load(f, allow_pickle=True)
    for k in z.files:
        name = f"{Path(f).stem.split('_')[-1]}:{k}"; rep["predictions"][name] = panel(z[k])
        r = rep["predictions"][name]; log(name, {x: round(r[x], 3) for x in ("heldout_auroc", "worst_class", "depth_strat_mean", "threshold_mean") + (("gated_macro_f1",) if gate is not None else ())}, {x: round(v, 2) for x, v in r["share"].items()})
Path(f"runs/p3/f5_stress_{u}.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
