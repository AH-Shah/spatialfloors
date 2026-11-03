"""F2 development region truth (PLAN §19, user 2026-11-03 "find another slide"): the other two Allen Zhuang animal 4
sections (C57BL6J-4.001, 4.002; the development animal, never a confirmation unit) become F2 development units, crops by
the §4 rule (data.crops.choose_crop on coordinates + major_brain_region only), as scripts/p3_allen4_dev.py did for 4.003.
Appends to manifests/crops/index.json with role "development" (F2 v3 selection only; F1 / F3 / F5 choices keep the
original development set). Usage (env floors): python scripts/v3_allen4_more.py"""
import json
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd
from spatialfloors.data.crops import choose_crop
from spatialfloors.data.units import ALLEN
from spatialfloors.provenance import stamp

A = ad.read_h5ad(ALLEN["animal4"], backed="r"); o = A.obs; man = Path("manifests/crops"); ix = json.loads((man / "index.json").read_text())
for sid in ("C57BL6J-4.001", "C57BL6J-4.002"):
    name = f"allen_animal4_{sid}"
    if name in ix: print("exists", name); continue
    m = (o.brain_section_label.astype(str) == sid).values; xy = np.asarray(A.obsm["X_spatial_coords"][np.flatnonzero(m)], float)[:, :2]
    idx, rec = choose_crop(xy, o.major_brain_region.values[m].astype(str))
    pd.Series(np.asarray(A.obs_names[m][idx]), name="cell_id").to_csv(man / f"{name}.csv.gz", index=False)
    ix[name] = stamp(dict(rec, dataset="Allen Zhuang animal 4", role="development", region_key="major_brain_region", note="PLAN §19: F2 v3 selection only"))
    print("saved", name, rec.get("n_crop"), rec.get("present"), rec.get("crop_shares"), flush=True)
(man / "index.json").write_text(json.dumps(ix, indent=1, default=str))
