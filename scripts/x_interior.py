"""Exploratory, post hoc (specified 2026-11-04 while P4 ran, before any confirmation output was read; no decision weight):
do the crops' cut edges penalise the spatial models? Each non-whole crop is an L-inf square cut from a larger slide
(manifests/crops/index.json, whole_slide False), so a cell near a cut edge loses its neighbours outside the crop, for every
model alike, but models with a longer spatial reach lose more. Re-score on the interior cells only, i.e. a scoring halo:
cells >= R um from every edge of the crop's bounding box, R in RS (whole-sample crops keep every cell). No refit.
  F2  balanced ARI of the saved partitions (CellCharter / Novae vs the smoothed-expression floor; all spatial)
  F5  held-out-marker AUROC of the saved fine labels (scANVI / Tangram are non-spatial; the floor's spatial prior is not)
Not re-scored: F1 (non-spatial models and metric), F3 (per-cell predictions are not saved; refits ~60 CPU h), F4 (whole
slides). D = learned - floor per crop at R = 0 (all scored cells, equals P5) and each R -> unit means (stats.aggregate_units)
-> mean over units. Usage (env floors): python scripts/x_interior.py <runs dir> <unit ...> -> <runs dir>/X_interior.{json,md}"""
import json, os, sys
from pathlib import Path
import anndata as ad, numpy as np
from spatialfloors import stats
from spatialfloors.floors import f5
from spatialfloors.metrics.f2 import scored_regions, ari_bal
from spatialfloors.metrics.f5 import unit_inputs, score_classes
from spatialfloors.provenance import stamp

R = Path(sys.argv[1]); U = sys.argv[2:]; RS = (0, 100, 250, 500)
L = json.loads(Path(os.environ.get("SF_LOCKED", "manifests/locked_v3.json")).read_text())["locked"]
sp_ = lambda u: "mouse" if u.startswith("allen") else "human"; cfg = lambda f, m, u: L[f][m][sp_(u)]["config"]
idx = json.loads(Path("manifests/crops/index.json").read_text())
donors = {f"gse303070_{d['slide'].replace('Tumor', 'Tumour')}_{d['pref']}": d["donor_cells"] for d in json.loads(Path("manifests/gse303070_donors.json").read_text())["check_vs_fov_csv"]}
MK0 = json.loads(Path("manifests/panels/index.json").read_text())["f5_heldout_markers"]
D = {}  # (fam, model) -> crop -> {R: D}
N = {}  # crop -> {R: share of the crop's cells kept}


def interior(u, xy, r):
    if r == 0 or idx[u].get("whole_slide"): return np.ones(len(xy), bool)
    lo, hi = xy.min(0), xy.max(0); return np.minimum(xy - lo, hi - xy).min(1) >= r


def ari_saved(path, key, m, reg):
    z = np.load(path, allow_pickle=True); Lb = np.asarray(z[key]).astype(str); Lb = Lb[None] if Lb.ndim == 1 else Lb
    y = reg[m]; return float(np.mean([ari_bal(y, l[m]) for l in Lb])) if any((y == c).sum() >= 200 for c in np.unique(y)) else np.nan


for u in U:
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); xy = a.obs[["x_um", "y_um"]].values.astype(float); IN = {r: interior(u, xy, r) for r in RS}
    N[u] = {r: round(float(IN[r].mean()), 3) for r in RS}
    fl = R / f"labels/f2_{u}_floor.npz"
    if fl.exists() and cfg("F2", "floor", u) in np.load(fl, allow_pickle=True).files:
        reg = np.asarray(a.obs.region.astype(str).values, dtype=object).astype(str); m0, _ = scored_regions(reg, idx[u]["present"])
        flv = {r: ari_saved(fl, cfg("F2", "floor", u), m0 & IN[r], reg) for r in RS}
        for mo in ("cellcharter", "novae"):
            lf = R / f"labels/f2_{u}_{mo}.npz"
            try: c = cfg("F2", mo, u)
            except KeyError: continue
            if lf.exists() and c in np.load(lf, allow_pickle=True).files:
                D.setdefault(("F2", mo), {})[u] = {r: ari_saved(lf, c, m0 & IN[r], reg) - flv[r] for r in RS}
    ff = R / f"labels/f5_{u}_floors.npz"
    if ff.exists() and cfg("F5", "floors", u) in np.load(ff, allow_pickle=True).files:
        MK = f5.for_unit(u, MK0); Xf, vf, gate, cls = unit_inputs(u, a)
        def auroc(path, key, m):
            fine = np.asarray(np.load(path, allow_pickle=True)[key]).astype(str)[m]
            return score_classes(f5.coarse(fine), Xf[m], vf, MK, gate[m] if gate is not None else None, cls)["heldout_auroc"]
        flv = {r: auroc(ff, cfg("F5", "floors", u), IN[r]) for r in RS}
        for mo in ("scanvi", "tangram"):
            lf = R / f"labels/f5_{u}_{mo}.npz"
            if lf.exists() and cfg("F5", mo, u) in np.load(lf, allow_pickle=True).files:
                D.setdefault(("F5", mo), {})[u] = {r: auroc(lf, cfg("F5", mo, u), IN[r]) - flv[r] for r in RS}
    print(u, N[u], {f"{k[0]}|{k[1]}": {r: round(v[u][r], 4) for r in RS} for k, v in D.items() if u in v}, flush=True)

out, md = dict(radii_um=RS, kept_share=N, pairs={}), ["| pair | units | " + " | ".join(f"mean D, R = {r} um" for r in RS) + " |", "|---|---|" + "---|" * len(RS)]
for (f, mo), d in D.items():
    per = {}
    for r in RS:
        ok = {c: v[r] for c, v in d.items() if np.isfinite(v[r])}; us, Du, _ = stats.aggregate_units(ok, lambda c: donors.get(c, c)); per[r] = dict(units=us, D=[float(x) for x in Du], mean=float(np.mean(Du)))
    out["pairs"][f"{f}|{mo}"] = dict(per_crop=d, by_radius=per)
    md.append(f"| {f} {mo} | {len(per[0]['units'])} | " + " | ".join(f"{per[r]['mean']:+.4f}" for r in RS) + " |")
(R / "X_interior.json").write_text(json.dumps(stamp(out), indent=1, default=float))
(R / "X_interior.md").write_text("Interior-only re-score (exploratory, scripts/x_interior.py). D = learned - floor; R = scoring halo from the crop's cut edges.\n\n" + "\n".join(md) + "\n")
print("\n".join(md))
