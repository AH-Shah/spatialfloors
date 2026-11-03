"""Cache a crop as h5ad (runs/p2/cache/<unit>.h5ad) so later P2 steps skip the 13 GB raw-slide read; prints geometry.
Usage (env floors): python scripts/p2_cache_crop.py <unit>"""
import sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from spatialfloors.data.units import load_crop

u = sys.argv[1]; f = Path(f"runs/p2/cache/{u}.h5ad")
a = load_crop(u); a.obs["region"] = a.obs.region.astype("category")  # None -> NaN (h5py cannot write None)
a.uns["unit"] = {k: v for k, v in a.uns["unit"].items() if v is not None}; a.write_h5ad(f, compression="gzip")
xy = a.obs[["x_um", "y_um"]].values; d = cKDTree(xy).query(xy, 16)[0]
print(u, a.n_obs, "fov", a.obs.fov.nunique(), "extent_um", np.ptp(xy, 0).round(), "nn1_um", np.median(d[:, 1]).round(2), "nn15_um", np.median(d[:, 15]).round(2))
