"""F2 development stress panel (PLAN §13; panel in spatialfloors.stress.f2): does a partition find the regions, or a coarse
geometry that the metric happens to reward? Development crops only. Floor variants computed here (k-means / full GMM,
3 seeds); learned models from runs/p3*/labels/*.npz.
Usage (env cellcharter, OMP_NUM_THREADS=4): python scripts/f2_stress.py <unit> -> runs/p3/f2_stress_<unit>.json"""
import glob, json, sys, time
from pathlib import Path
import numpy as np
from spatialfloors.floors.simple import expr_pca, neighbour_mean
from spatialfloors.metrics.f2 import kmeans_labels, gmm_labels
from spatialfloors.provenance import stamp
from spatialfloors.stress.f2 import load, Panel

u = sys.argv[1]; d = load(u); P = Panel(d); k = P.k; SEEDS = (0, 1, 2); log = lambda *s: print(*s, flush=True)
rep = dict(unit=u, k=k, k_fine=P.kf, n_boundary=int(P.bnd.sum()), n_interior=int(P.inn.sum()), blob_ari=P.blob_ari, partitions={})
rep["partitions"]["blob_xy_kmeans"] = P(P.blob[None]); log("blob", {key: round(v, 3) for key, v in rep["partitions"]["blob_xy_kmeans"].items()})


def zb(Mx): Mx = Mx - Mx.mean(0); return Mx / max(float(np.sqrt((Mx ** 2).sum(1).mean())), 1e-6)


E = expr_pca(d["X"], 30); feats = {"expr30": E}
for ell in (25., 100., 200.):
    C = neighbour_mean(d["xy"], d["fov"], E, ell)[0]; feats[f"ctx_ell{int(ell)}"] = C; feats[f"floor_ell{int(ell)}_lam0.8"] = np.hstack([np.sqrt(.2) * zb(E), np.sqrt(.8) * zb(C)])
for name, F in feats.items():
    for cn, clus in (("kmeans", kmeans_labels), ("gmm", gmm_labels)):
        t = time.time(); key = f"{name}_{cn}"; L = np.stack([clus(F, k, s) for s in SEEDS])
        rep["partitions"][key] = dict(**P(L), **P.detail(F, clus), s=round(time.time() - t, 1)); log(key, {kk: round(v, 3) for kk, v in rep["partitions"][key].items()})
for f in sorted(set(glob.glob(f"runs/p3*/labels/f2_{u}_*.npz") + glob.glob(f"runs/p3/relabel/labels/f2_{u}_*.npz"))):
    if "floor" in f: continue
    z = np.load(f, allow_pickle=True)
    for c in z.files:
        key = f"{Path(f).stem.split('_')[-1]}_{c}"; rep["partitions"][key] = P(z[c]); log(key, {kk: round(v, 3) for kk, v in rep["partitions"][key].items()})
Path(f"runs/p3/f2_stress_{u}.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
