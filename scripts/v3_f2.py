"""F2 floor v3 candidates on one unit (PLAN §19; development region truth: gse1_232 crop + Allen animal 4 sections
4.001 / 4.002 / 4.003). GMM full at true k, seeds metrics.f2.SEEDS (as the v2 runner), each raw and with a mean-field
Potts smoothing of its responsibilities (15-NN spatial graph, log-likelihoods tempered by the median top-1 / top-2
margin, beta 2, 5 sweeps). Scored by stress.f2.Panel (balanced ARI, boundary / interior ARI, min region recall, seed
stability) plus per-region recall (Hungarian match, seed 0).
  v2             domain_features(50 PCs, neighbour_mean ell 400 um, lambda .8) (locked)
  fixed400       v2 design with the operator fixed: grid Gaussian sigma 400 um across FOVs, no PCA (100 dims as v2)
  ms25_100_400   multi-scale grid Gaussian sigma 25 / 100 / 400 um, PCA 60
  ms50_200_800   sigma 50 / 200 / 800 um, PCA 60
  ms25_100_400+c100, ms50_200_800+c200   + Leiden-id composition at the middle scale
Usage (env floors, OMP_NUM_THREADS=3): [SF_OUT=runs/v3] python scripts/v3_f2.py <unit> -> <out>/f2_<unit>.json, labels npz"""
import json, os, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from sklearn import metrics as M
from sklearn.mixture import GaussianMixture
from spatialfloors.floors.simple import domain_features, domain_features_v3
from spatialfloors.metrics.f2 import SEEDS
from spatialfloors.provenance import stamp
from spatialfloors.stress.f2 import load, Panel

u = sys.argv[1]; t0 = time.time(); log = lambda *s: print(f"[{time.time() - t0:6.0f}s]", *s, flush=True); BETA = 2.
OUT = Path(os.environ.get("SF_OUT", "runs/v3")); OUT.mkdir(parents=True, exist_ok=True); ONLY = [c for c in os.environ.get("SF_CONFIG", "").split(",") if c]
d = load(u); P = Panel(d); k = P.k; xy, X, fov = d["xy"], d["X"], d["fov"]; A = cKDTree(xy).query(xy, 16)[1][:, 1:]
CFG = {"v2": lambda: domain_features(X, xy, fov, 50, 400., lam=.8)[0],
       "fixed400": lambda: domain_features_v3(X, xy, sigmas=(400.,), n_out=0)[0],
       "ms25_100_400": lambda: domain_features_v3(X, xy, sigmas=(25., 100., 400.))[0],
       "ms50_200_800": lambda: domain_features_v3(X, xy, sigmas=(50., 200., 800.))[0],
       "ms25_100_400+c100": lambda: domain_features_v3(X, xy, sigmas=(25., 100., 400.), comp_sigma=100.)[0],
       "ms50_200_800+c200": lambda: domain_features_v3(X, xy, sigmas=(50., 200., 800.), comp_sigma=200.)[0],
       # round 2 (after round 1 lost thin, cell-type-defined regions, e.g. Fiber_tracts on 4.001 .45 -> .27): more weight on
       # own expression, composition at a fine and a mid scale, less compression
       "ms50_200_800+c200_lam.6": lambda: domain_features_v3(X, xy, sigmas=(50., 200., 800.), comp_sigma=200., lam=.6)[0],
       "ms50_200_800+c50_200": lambda: domain_features_v3(X, xy, sigmas=(50., 200., 800.), comp_sigma=(50., 200.))[0],
       "ms50_200_800+c50_200_lam.6": lambda: domain_features_v3(X, xy, sigmas=(50., 200., 800.), comp_sigma=(50., 200.), lam=.6)[0],
       "ms50_200_800+c200_pc100": lambda: domain_features_v3(X, xy, sigmas=(50., 200., 800.), comp_sigma=200., n_out=100)[0]}


def potts(gm, F):
    lw = gm._estimate_weighted_log_prob(F); s2 = np.sort(lw, 1); lw = lw / max(float(np.median(s2[:, -1] - s2[:, -2])), 1e-6)
    R = np.exp(lw - lw.max(1, keepdims=True)); R /= R.sum(1, keepdims=True)
    for _ in range(5):
        z = lw + BETA * R[A].mean(1); z -= z.max(1, keepdims=True); R = np.exp(z); R /= R.sum(1, keepdims=True)
    return R.argmax(1)


def per_region(lab):
    C = M.cluster.contingency_matrix(P.y, np.asarray(lab).astype(str)[P.m]); r, c = linear_sum_assignment(-C); regs = np.unique(P.y)
    rec = {regs[i]: float(C[i, j] / C[i].sum()) for i, j in zip(r, c)}; return {str(x): dict(share=float((P.y == x).mean()), recall=rec.get(x, 0.)) for x in regs}


rep = dict(unit=u, k=k, seeds=list(SEEDS), beta=BETA, partitions={}); LAB = {}
for name, fn in CFG.items():
    if ONLY and name not in ONLY: continue
    t = time.time(); F = fn(); L, Ls = [], []
    for s in SEEDS:
        gm = GaussianMixture(k, covariance_type="full", reg_covar=1e-4, random_state=s).fit(F); L.append(gm.predict(F)); Ls.append(potts(gm, F))
    for tag, LL in ((name, L), (name + "+potts", Ls)):
        r = P(np.stack(LL)); r["per_region"] = per_region(LL[0]); rep["partitions"][tag] = r; LAB[tag] = np.stack(LL).astype(np.int16)
        log(tag, F.shape[1], {x: round(r[x], 3) for x in ("ari_bal", "ari", "ari_boundary", "ari_interior", "min_region_recall", "seed_stability")},
            {kk: round(v["recall"], 2) for kk, v in r["per_region"].items()}, f"{time.time() - t:.0f}s")
np.savez_compressed(OUT / f"f2_{u}_labels.npz", **LAB)
(OUT / f"f2_{u}.json").write_text(json.dumps(stamp(rep), indent=1, default=float)); print("V3F2_DONE", u, flush=True)
