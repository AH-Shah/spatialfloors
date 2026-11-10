"""Tier 2 (PROTOCOL §5): F5 floor with spatial-prior strength beta {.5, 2} on the unit's tier-1 InSituType
log-likelihoods ($SF_T1/ist_<unit>/logliks.csv); beta 0 (InSituType) and beta 1 (the locked floor) are tier-1 outputs.
SF_BETAS overrides the betas (check: beta 1 reproduces the locked floor).
Usage (env floors): python scripts/t2_f5_prior.py <unit> -> $SF_OUT/f5_<unit>_prior.json"""

import json, os, sys
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd
from spatialfloors.floors import f5
from spatialfloors.metrics.f5 import unit_inputs, score_classes
from spatialfloors.provenance import stamp

u = sys.argv[1]
OUT = Path(os.environ.get("SF_OUT", "runs/p4/t2/f5_floor_prior"))
OUT.mkdir(parents=True, exist_ok=True)
MK = f5.for_unit(u, json.loads(Path("manifests/panels/index.json").read_text())["f5_heldout_markers"])
a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad")
Xf, vf, gate, cls = unit_inputs(u, a)
xy = a.obs[["x_um", "y_um"]].values.astype(float)
LL = pd.read_csv(Path(os.environ.get("SF_T1", "runs/p4")) / f"ist_{u}" / "logliks.csv")
rep = dict(unit=u, model="prior", grid={})
LAB = {}
for b in map(float, os.environ.get("SF_BETAS", ".5,2").split(",")):
    k = f"prior_b{b:g}"
    LAB[k] = f5.spatial_prior(LL.values, list(LL.columns), xy, beta=b)
    rep["grid"][k] = score_classes(f5.coarse(LAB[k]), Xf, vf, MK, gate, cls)
    print(u, k, round(rep["grid"][k]["heldout_auroc"], 4), flush=True)
(OUT / "labels").mkdir(exist_ok=True)
np.savez_compressed(OUT / "labels" / f"f5_{u}_prior.npz", **LAB)
(OUT / f"f5_{u}_prior.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
