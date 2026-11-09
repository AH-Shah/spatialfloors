"""Exploratory (after results-v1; external review 2026-11-09, item F5.1/F5.2; no decision weight).
F5 held-out-marker AUROC at a matched call rate. Per F5 crop and class c, the locked floor (InSituType ML + spatial
prior) calls c on its top-n cells by class posterior (sum of the fine-type posteriors of floors.f5.spatial_prior, same
arithmetic), n = the number of cells the comparator (locked scANVI / Tangram) called c. Per-class AUROC by
metrics.f5.heldout_auroc. D = comparator - floor at the comparator's call rate; unit means and one-sided 95% bounds as
tier 1. No model is refitted; reads runs/p4 labels, ist_<crop>/logliks.csv and the cached crops.
Usage (env floors): python scripts/x_f5_matched.py runs/p4 <confirmation crops> -> runs/p4/X_f5_matched.{json,md}"""
import sys, json
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
from scipy.spatial import cKDTree
from spatialfloors import stats
from spatialfloors.data.genes import filter_reason
from spatialfloors.floors import f5 as F5
from spatialfloors.metrics.f5 import CLASSES, heldout_auroc
from spatialfloors.provenance import commit

R = Path(sys.argv[1]); U = [u for u in sys.argv[2:] if not u.startswith(("allen", "spatch"))]
L = json.loads(Path("manifests/locked_v3.json").read_text())["locked"]; cfg = lambda m: L["F5"][m]["human"]["config"]
assert cfg("floors") == "insitutype_spatial"
donors = {f"gse303070_{d['slide'].replace('Tumor', 'Tumour')}_{d['pref']}": d["donor_cells"] for d in json.loads(Path("manifests/gse303070_donors.json").read_text())["check_vs_fov_csv"]}
unit_of = lambda c: donors.get(c, c)
MK = json.loads(Path("manifests/panels/index.json").read_text())["f5_heldout_markers"]; COARSE0 = dict(F5.COARSE)


def posterior(ll, xy, k=15, beta=1., sweeps=5):
    """floors.f5.spatial_prior, returning the posterior instead of its argmax (identical arithmetic)."""
    ll = np.asarray(ll, float); s2 = np.sort(ll, 1); z0 = ll / max(float(np.median(s2[:, -1] - s2[:, -2])), 1e-6)
    nb = cKDTree(xy).query(xy, k + 1)[1][:, 1:]; P = np.exp(z0 - z0.max(1, keepdims=True)); P /= P.sum(1, keepdims=True)
    for _ in range(sweeps):
        z = z0 + beta * P[nb].mean(1); z -= z.max(1, keepdims=True); P = np.exp(z); P /= P.sum(1, keepdims=True)
    return P


def auroc_top(score, n, c, X, vn, mk):
    pred = np.full(len(score), "x", dtype=object); pred[np.argsort(-score, kind="stable")[:n]] = c
    return heldout_auroc(pred, X, vn, mk, classes=(c,))[1][c]["auroc"]


rows = []
for u in U:
    F5.COARSE.clear(); F5.COARSE.update(COARSE0); mk = F5.for_unit(u, MK)
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); keep = filter_reason(a.var_names, "human") == ""
    X, vn, xy = sp.csr_matrix(a.X[:, keep]), list(map(str, a.var_names[keep])), a.obs[["x_um", "y_um"]].values.astype(float); del a
    LL = pd.read_csv(R / f"ist_{u}" / "logliks.csv"); P = posterior(LL.values, xy); names = np.asarray(LL.columns)
    floor_lab = np.load(R / f"labels/f5_{u}_floors.npz", allow_pickle=True)["insitutype_spatial"].astype(str)
    assert (names[P.argmax(1)] == floor_lab).all(), f"{u}: posterior argmax differs from the saved floor labels"
    cls = F5.coarse(names); Pc = {c: P[:, cls == c].sum(1) for c in CLASSES}; fl = F5.coarse(floor_lab)
    f_reg = heldout_auroc(fl, X, vn, mk)[1]
    for m in ("scanvi", "tangram"):
        lab = F5.coarse(np.load(R / f"labels/f5_{u}_{m}.npz", allow_pickle=True)[cfg(m)].astype(str)); m_reg = heldout_auroc(lab, X, vn, mk)[1]
        reg = json.loads((R / f"f5_{u}_{m}.json").read_text())["grid"][cfg(m)]["per_class"]
        for c in CLASSES:
            assert abs(m_reg[c]["auroc"] - reg[c]["auroc"]) < 1e-9, f"{u} {m} {c}: recomputed AUROC differs from the run's"
            n_m, n_f = int((lab == c).sum()), int((fl == c).sum())
            rows.append(dict(crop=u, comparator=m, cls=c, n_comp=n_m, n_floor=n_f, auroc_comp=m_reg[c]["auroc"], auroc_floor=f_reg[c]["auroc"],
                             auroc_floor_ownrate=auroc_top(Pc[c], n_f, c, X, vn, mk), auroc_floor_matched=auroc_top(Pc[c], n_m, c, X, vn, mk)))
    print(u, "done", flush=True)

df = pd.DataFrame(rows); margin = L["F5"]["margin"]; out = dict(commit=commit(), margin=margin, per_crop_class=rows, tests={})
md = ["# X_f5_matched: F5 AUROC at the comparator's call rate (exploratory, after results-v1; no decision weight)", "",
      f"commit {commit()}; margin {margin}. Floor top-n by class posterior, n = comparator's calls of that class.", "",
      "| Comparison | Classes | Units | D | One-sided 95% bounds | Mean first | Mean second |", "|---|---|---|---|---|---|---|"]
for m in ("scanvi", "tangram"):
    for cl in ("mean", *CLASSES):
        for second in ("auroc_floor", "auroc_floor_ownrate", "auroc_floor_matched"):
            s = df[df.comparator == m] if cl == "mean" else df[(df.comparator == m) & (df.cls == cl)]
            g = s.groupby("crop")[["auroc_comp", second]].mean()
            units, D, strata = stats.aggregate_units(dict(zip(g.index, g.auroc_comp - g[second])), {c: unit_of(c) for c in g.index})
            t = stats.family_test(D, strata, margin=margin); key = f"{m} - {second} [{cl}]"
            out["tests"][key] = dict(n_units=len(D), est=t["est"], lower=t["lower"], upper=t["upper"], noninferior=t["noninferior"],
                                     strata_means={s_: float(D[np.asarray(strata) == s_].mean()) for s_ in np.unique(strata)})
            md.append(f"| {m} − {second} | {cl} | {len(D)} | {t['est']:+.4f} | {t['lower']:+.4f}, {t['upper']:+.4f} | {g.auroc_comp.mean():.3f} | {g[second].mean():.3f} |")
sh = df.groupby(["comparator", "cls"])[["n_comp", "n_floor"]].sum()
md += ["", "Calls summed over crops (comparator / floor): " + "; ".join(f"{m} {c} {int(r.n_comp)} / {int(r.n_floor)}" for (m, c), r in sh.iterrows())]
md += ["", "Strata means of D (mean over classes): " + "; ".join(f"{k}: " + ", ".join(f"{s_} {v:+.4f}" for s_, v in out['tests'][k]['strata_means'].items())
                                                       for k in out["tests"] if k.endswith("[mean]"))]
(R / "X_f5_matched.json").write_text(json.dumps(out, indent=1)); (R / "X_f5_matched.md").write_text("\n".join(md) + "\n"); print("\n".join(md))
