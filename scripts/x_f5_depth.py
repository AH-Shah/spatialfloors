"""Exploratory (PROTOCOL §6, F5 stress panel on the confirmation units): the F5 score with library size held fixed.
Per crop and locked configuration (floor, scANVI, Tangram): the registered AUROC (recomputed and checked
against the registered score), the depth-stratified AUROC (the same AUROC within each library-size quintile, mean over
quintiles: depth is held fixed) and the share of each class's calls in the deepest quintile. Then D = learned - floor
per unit and the registered test (stats.family_test, platform strata, F5 margin) on the depth-stratified score.
Usage (env floors, cwd = the runs root holding runs/p2/cache and runs/p4):
  python scripts/x_f5_depth.py <out dir>    -> <out dir>/X_f5_depth.{json,md}"""

import json, sys
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd
from spatialfloors import stats
from spatialfloors.floors import f5
from spatialfloors.metrics.f5 import CLASSES, heldout_auroc, unit_inputs
from spatialfloors.provenance import commit

PUB = Path(__file__).resolve().parents[1]
OUT = Path(sys.argv[1])
L = json.loads((PUB / "manifests/locked_v3.json").read_text())["locked"]["F5"]
CFG = {m: L[m]["human"]["config"] for m in ("floors", "scanvi", "tangram")}
MK_ALL = json.loads((PUB / "manifests/panels/index.json").read_text())["f5_heldout_markers"]
U = pd.read_csv(PUB / "registration/units/units.tsv", sep="\t")
CROPS = [r.crop for r in U.itertuples() if r.role == "confirmation" and "F5" in r.families.split(",")]
donors = {
    f"gse303070_{d['slide'].replace('Tumor', 'Tumour')}_{d['pref']}": d["donor_cells"]
    for d in json.loads((PUB / "manifests/gse303070_donors.json").read_text())["check_vs_fov_csv"]
}
unit_of = lambda c: donors.get(c, c)

per_crop = {}
for u in CROPS:
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad")
    X, vn, _, _ = unit_inputs(u, a)
    MK = f5.for_unit(u, MK_ALL)
    lib = np.maximum(np.asarray(X.sum(1)).ravel(), 1)
    dq = np.digitize(np.log(lib), np.quantile(np.log(lib), [0.2, 0.4, 0.6, 0.8]))  # library-size quintile 0..4
    r = {}
    for m, cfg in CFG.items():
        c = np.asarray(f5.coarse(np.load(f"runs/p4/labels/f5_{u}_{m}.npz", allow_pickle=True)[cfg])).astype(str)
        h, per = heldout_auroc(c, X, vn, MK)
        reg = json.loads(Path(f"runs/p4/f5_{u}_{m}.json").read_text())["grid"][cfg]["heldout_auroc"]
        assert abs(h - reg) < 1e-9, (u, m, h, reg)  # same score as the registered run
        strat = [heldout_auroc(c[dq == q], X[dq == q], vn, MK)[0] for q in range(5)]
        r[m] = dict(
            heldout_auroc=h,
            depth_strat=float(np.mean(strat)),
            by_quintile=strat,
            top_quintile_share={k: float((dq[c == k] == 4).mean()) if (c == k).any() else None for k in CLASSES},
        )
    per_crop[u] = r
    print(u, {m: (round(v["heldout_auroc"], 3), round(v["depth_strat"], 3)) for m, v in r.items()}, flush=True)

out = dict(commit=commit(), configs=CFG, per_crop=per_crop, pairs={})
md = [
    "# F5 with depth held fixed (exploratory)",
    "",
    "Depth-stratified = held-out-marker AUROC within each library-size quintile, mean over quintiles. "
    "Registered scores recomputed exactly on every crop.",
    "",
    "| Pair | Score | Units | Floor | Learned | D (est) | one-sided 95% bounds | Non-inferior | Floor superior |",
    "|---|---|---|---|---|---|---|---|---|",
]
for m in ("scanvi", "tangram"):
    out["pairs"][m] = {}
    for key in ("heldout_auroc", "depth_strat"):
        D_c = {u: per_crop[u][m][key] - per_crop[u]["floors"][key] for u in CROPS}
        units, D, strata = stats.aggregate_units(D_c, {c: unit_of(c) for c in CROPS})
        t = stats.family_test(D, strata, margin=L["margin"])
        fl = stats.aggregate_units({u: per_crop[u]["floors"][key] for u in CROPS}, {c: unit_of(c) for c in CROPS})[1]
        le = stats.aggregate_units({u: per_crop[u][m][key] for u in CROPS}, {c: unit_of(c) for c in CROPS})[1]
        out["pairs"][m][key] = {k: t[k] for k in ("est", "lower", "upper", "noninferior", "floor_superior")} | dict(
            n_units=len(D), floor_mean=float(fl.mean()), learned_mean=float(le.mean()), per_unit=dict(zip(units, D))
        )
        s = out["pairs"][m][key]
        md.append(
            f"| F5 {m} | {'registered' if key == 'heldout_auroc' else 'depth-stratified'} | {s['n_units']} | "
            f"{s['floor_mean']:.3f} | {s['learned_mean']:.3f} | {s['est']:+.4f} | {s['lower']:+.4f}, {s['upper']:+.4f} | "
            f"{'yes' if s['noninferior'] else 'no'} | {'yes' if s['floor_superior'] else 'no'} |"
        )
md += [
    "",
    "Share of each class's calls in the deepest library-size quintile (mean over crops; 0.20 = depth-neutral):",
    "",
]
md += ["| Model | " + " | ".join(CLASSES) + " |", "|---|" + "---|" * len(CLASSES)]
for m in CFG:
    sh = [np.nanmean([per_crop[u][m]["top_quintile_share"][k] or np.nan for u in CROPS]) for k in CLASSES]
    md.append(f"| {m} | " + " | ".join(f"{x:.2f}" for x in sh) + " |")
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "X_f5_depth.json").write_text(json.dumps(out, indent=1, default=float))
(OUT / "X_f5_depth.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
print("X_F5_DEPTH_DONE", flush=True)
