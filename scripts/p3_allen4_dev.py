"""Mouse development unit (review 2026-10-28): Allen Zhuang animal 4 (3 lateral sagittal sections, unused since P0 replaced
it by Yao 638850 for confirmation: no thalamus) becomes the MERFISH development unit for F1-F3 (and F2: >= 2 present
regions). Section = the one with most present regions, 'n/a' excluded (ties: most cells); crop rule as data.crops (coordinates + region
labels only). Appends to manifests/crops/index.json. Usage (env floors): python scripts/p3_allen4_dev.py"""
import json
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd
from spatialfloors.data.crops import choose_crop
from spatialfloors.data.units import ALLEN
from spatialfloors.metrics.f2 import JUNK
from spatialfloors.provenance import stamp

A = ad.read_h5ad(ALLEN["animal4"], backed="r"); o = A.obs; best = None
for sid in sorted(o.brain_section_label.astype(str).unique()):
    m = (o.brain_section_label.astype(str) == sid).values; xy = np.asarray(A.obsm["X_spatial_coords"][np.flatnonzero(m)], float)[:, :2]
    idx, rec = choose_crop(xy, o.major_brain_region.values[m].astype(str)); key = (len([r for r in rec.get("present", []) if r not in JUNK]), int(m.sum()))  # real regions only
    print(sid, int(m.sum()), rec.get("present"), rec.get("crop_shares"), flush=True)
    if best is None or key > best[0]: best = (key, sid, A.obs_names[m][idx], rec)
_, sid, ids, rec = best; name = f"allen_animal4_{sid}"; man = Path("manifests/crops")
pd.Series(np.asarray(ids), name="cell_id").to_csv(man / f"{name}.csv.gz", index=False)
ix = json.loads((man / "index.json").read_text()); ix[name] = stamp(dict(rec, dataset="Allen Zhuang animal 4", role="development", region_key="major_brain_region"))
(man / "index.json").write_text(json.dumps(ix, indent=1, default=str))  # the index keeps its P0 stamp; the new record carries its own; print("saved", name, rec.get("n_crop"), rec.get("crop_shares"))
