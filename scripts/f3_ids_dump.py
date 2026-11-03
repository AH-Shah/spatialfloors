"""Freeze the shared F3 receiver ids (floors.f3.prepare: Leiden on Expr50, env floors) per unit -> runs/p2/f3_ids/<unit>.npy,
so models run in other envs (NCEM: external/ncem_env, older scanpy / igraph) get the same receivers, lineages and masks.
Usage (env floors): python scripts/f3_ids_dump.py <unit ...>"""
import json, sys
from pathlib import Path
import anndata as ad, numpy as np
from spatialfloors.floors import f3

Path("runs/p2/f3_ids").mkdir(parents=True, exist_ok=True)
for u in sys.argv[1:]:
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); g = a.var_names.get_indexer(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"])
    D = f3.prepare(a.X[:, g], a.obs[["x_um", "y_um"]].values, a.obs.fov.values); np.save(f"runs/p2/f3_ids/{u}.npy", D["ids"]); print(u, D["K"], len(f3.receivers(D)), flush=True)
