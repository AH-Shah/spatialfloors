"""Exploratory scale curves: windows of one GSE1 confirmation slide at several cell counts, cut from the slide's lean
whole-slide cache (scripts/xw_prepare.py). Window = data.crops.choose_crop(xy, regions, n_target=n) on the whole slide
(the registered crop rule at another n_target); sizes with n_slide <= 1.1 n are skipped (the whole slide is the
whole-slide arm's unit gse1_<sid>w). Each window keeps the whole slide's gene panel and count split (rows of
gse1_<sid>w.h5ad and gse1_<sid>w_split.npz), so only the number of cells changes along a curve.
Writes, per window u = gse1_<sid>s<n>: runs/p2/cache/<u>.h5ad, runs/p2/cache/<u>_split.npz, manifests/crops/<u>.csv.gz,
manifests/crops/index.json[u], manifests/panels/<u>.json (the slide panel), manifests/panels/index.json[u].
Usage (env floors, round-1 runs root as working directory): python scripts/x_scale_prepare.py <sid> [n ...]"""

import json, sys
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
from spatialfloors.data.crops import choose_crop
from spatialfloors.provenance import stamp

SIZES = (12500, 25000, 50000, 100000, 200000, 400000)
sid = sys.argv[1]
sizes = [int(x) for x in sys.argv[2:]] or SIZES
w = f"gse1_{sid}w"
a = ad.read_h5ad(f"runs/p2/cache/{w}.h5ad")
Xa = sp.load_npz(f"runs/p2/cache/{w}_split.npz").tocsr()
pan = json.loads(Path(f"manifests/panels/{w}.json").read_text())
xy = a.obs[["x_um", "y_um"]].values.astype(float)
reg = a.obs.region.astype(object).where(a.obs.region.notna(), None).values
reg = np.array([r if isinstance(r, str) and r not in ("", "nan") else None for r in reg], dtype=object)
for n in sizes:
    if a.n_obs <= 1.1 * n:
        print(sid, n, "skipped: slide has", a.n_obs, "cells", flush=True)
        continue
    u = f"gse1_{sid}s{n}"
    idx, rec = choose_crop(xy, reg, n_target=n)
    idx = np.sort(idx)
    b = a[idx].copy()
    b.uns["unit"] = dict(name=u, slide=w, n_target=n)
    b.write_h5ad(f"runs/p2/cache/{u}.h5ad", compression="gzip")
    sp.save_npz(f"runs/p2/cache/{u}_split.npz", Xa[idx])
    pd.Series(a.obs_names[idx].values, name="cell_id").to_csv(f"manifests/crops/{u}.csv.gz", index=False)
    rec.update(dataset="GSE1", role="exploratory_scale", region_key="typ", note=f"window of {w} (scripts/x_scale_prepare.py)")
    p = dict(pan, unit=u, n_cells=int(len(idx)), window_of=w)
    for f, key, val in (
        (Path("manifests/crops/index.json"), u, rec),
        (Path("manifests/panels/index.json"), u, {k: p[k] for k in ("mode", "n_input", "n_panel", "removed", "n_cells")}),
    ):
        d = json.loads(f.read_text())
        d[key] = stamp(val)
        f.write_text(json.dumps(d, indent=1, default=float))
    Path(f"manifests/panels/{u}.json").write_text(json.dumps(stamp(p), indent=1, default=float))
    print(sid, n, u, rec["n_crop"], rec.get("crop_shares"), flush=True)
print("X_SCALE_PREPARE_DONE", sid)
