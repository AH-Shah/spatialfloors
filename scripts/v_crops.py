"""Crop manifests for the validation slides and the Bruker units by the registered crop rule (data.crops.choose_crop);
appended to manifests/crops/index.json with role="validation" (existing entries untouched).
Usage (env floors): python scripts/v_crops.py <unit> ...   units: gse1_231, gse1_231w, bruker_*"""

import json, sys
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
from spatialfloors.data.crops import choose_crop
from spatialfloors.data.units import RAW, bruker_obs
from spatialfloors.provenance import stamp

OUT = Path("manifests/crops")
idx = json.loads((OUT / "index.json").read_text())
for u in sys.argv[1:]:
    if u == "gse1_231w":  # crop check: every cell of the slide
        a = ad.read_h5ad(RAW / "GSE1/cosmx_wt_231.h5ad", backed="r")
        ids = a.obs.index
        reg = a.obs["typ"].astype(str).str.extract(r"(REF|TVA|CRC)", expand=False)
        pd.Series(np.asarray(ids), name="cell_id").to_csv(OUT / f"{u}.csv.gz", index=False)
        idx[u] = dict(
            n_slide=len(ids),
            n_crop=len(ids),
            whole_slide=True,
            slide_shares=reg.value_counts(normalize=True).round(4).to_dict(),
            present=sorted(reg.dropna().unique()),
            dataset="GSE1",
            role="validation",
            region_key="typ",
            note="crop check: whole slide of gse1_231",
        )
        a.file.close()
        print(u, len(ids))
        continue
    if u == "gse1_231":  # same patient as development gse1_232 (P0 exclusion); used for validation only, flagged
        a = ad.read_h5ad(RAW / "GSE1/cosmx_wt_231.h5ad", backed="r")
        o = a.obs
        ids = o.index
        xy = o[["CenterX_global_px", "CenterY_global_px"]].values.astype(float)
        reg = o["typ"].astype(str).str.extract(r"(REF|TVA|CRC)", expand=False).values
        extra = dict(dataset="GSE1", role="validation", region_key="typ", note="same patient as development gse1_232")
        a.file.close()
    else:
        o = bruker_obs(u)
        ids = o.index
        xy = o[["x_px", "y_px"]].values.astype(float)
        reg = None  # no regions: densest contiguous square (P0 rule)
        extra = dict(
            dataset="Bruker CosMx WTx", role="validation" if u == "bruker_breast" else "confirmation", region_key=None
        )  # roles_v2.json
    sel, rec = choose_crop(xy, reg)
    pd.Series(np.asarray(ids[sel]), name="cell_id").to_csv(OUT / f"{u}.csv.gz", index=False)
    idx[u] = dict(rec, **extra)
    print(u, {k: idx[u].get(k) for k in ("n_slide", "n_crop", "crop_shares", "min_share")}, flush=True)
(OUT / "index.json").write_text(json.dumps(stamp(idx), indent=1, default=float))
