"""Per-region Hungarian recall over seeds: floor GMM anchors (5 seeds) vs CellCharter (P3 configs + embedding seeds 1, 2)
and Novae. Which region is the minimum, and is it the same region for floor and learned? Usage: python scripts/f2_region_recall.py <unit>"""
import glob, sys
import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn import metrics as M
from spatialfloors.floors.simple import expr_pca, neighbour_mean
from spatialfloors.metrics.f2 import gmm_labels
from spatialfloors.stress.f2 import load, Panel

u = sys.argv[1]; d = load(u); P = Panel(d); regs = np.unique(P.y)
zb = lambda X: (X - X.mean(0)) / max(float(np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean())), 1e-6)


def rec(lab):
    Ct = M.cluster.contingency_matrix(P.y, np.asarray(lab).astype(str)[P.m]); r, c = linear_sum_assignment(-Ct); out = np.zeros(len(regs))
    for i, j in zip(r, c): out[i] = Ct[i, j] / Ct[i].sum()
    return out


rows = {}; E = zb(expr_pca(d["X"], 30))
for ell, lam in (() if "learned" in sys.argv else ((100., .5), (100., .8), (200., .8), (400., .8))):
    C = zb(neighbour_mean(d["xy"], d["fov"], E, ell)[0]); F = np.hstack([np.sqrt(1 - lam) * E, np.sqrt(lam) * C]).astype(np.float32)
    rows[f"floor ell{ell:g} lam{lam:g}"] = [rec(gmm_labels(F, P.k, s)) for s in range(5)]
for f in glob.glob(f"runs/p3/relabel/labels/f2_{u}_*.npz") + glob.glob(f"runs/landscape/labels/f2_{u}_*_none_s*.npz"):
    z = np.load(f, allow_pickle=True)
    for c in z.files:
        L = z[c] if z[c].ndim == 2 else z[c][None]; rows.setdefault(f.split("/")[-1].replace(f"f2_{u}_", "")[:-4] + ":" + c, []).extend(rec(l) for l in L)
print("row".ljust(44), " ".join(r[:9].rjust(9) for r in regs), " min (region, per seed)")
for n, R in rows.items():
    R = [x for x in R if x is not None]
    if not R: continue
    A = np.array(R); print(n[:44].ljust(44), " ".join(f"{v:9.2f}" for v in A.mean(0)), " ", ", ".join(regs[a.argmin()][:8] + f" {a.min():.2f}" for a in A))
