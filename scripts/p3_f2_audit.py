"""F2 floor audit (PLAN §12 floor-optimisation rule; user 2026-10-28: CellCharter beats the floor on dev, find the exact
source). Development crops only. Separates three candidate sources of the gap:
  (1) information: supervised ceiling (logistic regression, held-out FOVs or quadrants -> balanced accuracy) per feature set
  (2) clustering: k-means vs full-covariance GMM on the same features, all 10 seeds listed (CellCharter uses a GMM)
  (3) neighbourhood: floor Gaussian-kernel mean vs CellCharter's aggregation (Delaunay, long links > 99th pct removed,
      mean of each exact hop 1..L, concatenated) applied to Expr PCA instead of scVI (architecture-matched floor)
plus data structure: region purity of FOVs, neighbour region agreement, Leiden-id composition gap between regions.
Usage (env cellcharter): python scripts/p3_f2_audit.py <unit> -> runs/p3/f2_audit_<unit>.json"""
import json, sys, time
from pathlib import Path
import anndata as ad, numpy as np, scipy.sparse as sp
from scipy.spatial import Delaunay, cKDTree
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import adjusted_rand_score, balanced_accuracy_score
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from spatialfloors.data.units import index
from spatialfloors.floors import bleedsafe as bs
from spatialfloors.floors.simple import expr_pca, lognorm, neighbour_mean
from spatialfloors.metrics.f2 import scored_regions, SEEDS
from spatialfloors.metrics.f3 import quadrant_blocks
from spatialfloors.provenance import stamp

u = sys.argv[1]; a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); p = json.loads(Path(f"manifests/panels/{u}.json").read_text())
X = a.X[:, a.var_names.get_indexer(p["genes"])]; xy = a.obs[["x_um", "y_um"]].values.astype(float); fov = a.obs.fov.astype(str).values
m, regl = scored_regions(a.obs.region.astype(str).values, index()[u]["present"]); k = len(np.unique(regl)); rep = dict(unit=u, k=k, n_scored=int(m.sum()))
log = lambda *s: print(*s, flush=True)

# data structure
ids = bs.leiden_ids(bs.expr50(lognorm(X), 0), seed=0); nb = cKDTree(xy).query(xy, 16)[1][:, 1:]; R = np.where(m, a.obs.region.astype(str).values, "")
pur = [max(np.unique(R[(fov == f) & m], return_counts=True)[1]) / ((fov == f) & m).sum() for f in np.unique(fov) if ((fov == f) & m).sum()]
agree = float(np.mean([(R[nb[i]] == R[i]).mean() for i in np.flatnonzero(m)[::10]]))
comp = {r: np.bincount(ids[R == r], minlength=ids.max() + 1) / (R == r).sum() for r in np.unique(regl)}
tv = float(np.mean([.5 * np.abs(comp[r1] - comp[r2]).sum() for i, r1 in enumerate(comp) for r2 in list(comp)[i + 1:]]))
rep["structure"] = dict(n_fov=len(np.unique(fov)), fov_purity_median=float(np.median(pur)), fov_pure_share=float(np.mean(np.array(pur) > .95)),
                        nn15_region_agreement=agree, n_ids=int(ids.max() + 1), composition_tv_mean_pair=tv)
log("structure", rep["structure"])


def hop_agg(E, L):
    """CellCharter aggregate_neighbors (mean) on a Delaunay graph with links > 99th percentile removed."""
    tri = Delaunay(xy); e = np.vstack([tri.simplices[:, [i, j]] for i, j in ((0, 1), (1, 2), (0, 2))]); e = np.unique(np.sort(e, 1), axis=0)
    d = np.linalg.norm(xy[e[:, 0]] - xy[e[:, 1]], axis=1); e = e[d <= np.percentile(d, 99)]; n = len(xy)
    A = sp.csr_matrix((np.ones(2 * len(e)), (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])), (n, n)); reach = sp.identity(n, format="csr", dtype=bool); front = reach
    out = [E]
    for _ in range(L):
        nxt = ((front.astype(np.float32) @ A) > 0).astype(bool); nxt = nxt - nxt.multiply(reach); nxt.eliminate_zeros()
        c = np.asarray(nxt.sum(1)).ravel(); out.append(np.asarray(sp.diags(1 / np.maximum(c, 1)) @ (nxt.astype(np.float32) @ E))); reach = reach + nxt; front = nxt
    return np.hstack(out).astype(np.float32)


def zblock(M): M = M - M.mean(0); return M / max(float(np.sqrt((M ** 2).sum(1).mean())), 1e-6)


E30 = expr_pca(X, 30); Efc = bs.group_center(E30, fov).astype(np.float32)
feats = {"expr30": E30, "expr30_fovcentred": Efc}
for ell in (25., 100., 200., 400.):
    C = neighbour_mean(xy, fov, E30, ell)[0]; feats[f"ctx_ell{int(ell)}"] = C
    for lam in (.5, .8): feats[f"floor_ell{int(ell)}_lam{lam}"] = np.hstack([np.sqrt(1 - lam) * zblock(E30), np.sqrt(lam) * zblock(C)])
for L in (1, 3): feats[f"hopagg_L{L}"] = hop_agg(E30, L)

blocks = quadrant_blocks(xy)[0] if len(np.unique(fov)) < 4 else fov
rep["features"] = {}
for name, F in feats.items():
    t = time.time(); Fs = StandardScaler().fit_transform(F[m]); y = regl; b = blocks[m]; pr = np.empty(len(y), dtype=object)
    for g in np.unique(b):  # leave-one-block-out (FOVs; quadrants when the crop is one FOV)
        te = b == g
        if len(np.unique(y[~te])) < 2: pr[te] = y[~te][0]; continue
        pr[te] = LogisticRegression(max_iter=300, C=1.).fit(Fs[~te][::3], y[~te][::3]).predict(Fs[te])
    sup = float(balanced_accuracy_score(y, pr.astype(str)))
    km = [adjusted_rand_score(y, KMeans(k, n_init=10, random_state=s).fit_predict(F)[m]) for s in SEEDS]
    gm = [adjusted_rand_score(y, GaussianMixture(k, covariance_type="full", random_state=s, reg_covar=1e-4).fit(F[::2]).predict(F)[m]) for s in SEEDS]
    rep["features"][name] = dict(dim=int(F.shape[1]), supervised_bacc=sup, kmeans_ari=float(np.mean(km)), kmeans_seeds=np.round(km, 3).tolist(),
                                 gmm_ari=float(np.mean(gm)), gmm_seeds=np.round(gm, 3).tolist(), s=round(time.time() - t, 1))
    log(name, {k_: v for k_, v in rep["features"][name].items() if not k_.endswith("seeds")}, "gmm seeds", np.round(gm, 2).tolist())
Path(f"runs/p3/f2_audit_{u}.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
