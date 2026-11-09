"""Memory-lean twins for whole slides (exploratory whole-slide arm, results/WHOLESLIDE_PLAN.md; active only with
SF_LEAN=1, so no registered computation changes). The registered floors build a dense cells x panel matrix
(741k x 2,000 float64 = 11.9 GB, plus copies) and the v3 F3 rings a dense cells x 1,500 neighbour table (17.8 GB
transient on 741k cells). Here:
- PCA runs sklearn's own algorithm for PCA(svd_solver="randomized", random_state=seed) (sklearn 1.9.0
  _randomized_svd: RandomState normal start, n_oversamples 10, n_iter 7 with LU normalisation, QR, SVD of Q^T A,
  svd_flip on Vt), with every product against the centred matrix computed over row chunks that are rebuilt by the
  registered row operations, in the registered dtype (float64 for expr_embed, float32 for expr_pca). Results equal
  the registered ones up to floating-point summation order. Checked by tests/test_lean.py and scripts/xw_check.py.
- Rings are queried in row chunks and kept as a sparse cells x cells 0/1 matrix (bleedsafe.kernel_mean accepts it)."""
import os
import numpy as np
import scipy.sparse as sp

CHUNK = 20000


def on():
    return os.environ.get("SF_LEAN") == "1"


def _rpca(rows_fn, n, p, n_pc, dtype, seed=0, n_oversamples=10):
    """rows_fn(a, b) -> dense rows a:b of the matrix the registered code hands to PCA. Returns PCA(n_pc,
    svd_solver="randomized", random_state=seed).fit_transform of it, as float32, without building it."""
    from scipy import linalg
    from sklearn.utils.extmath import svd_flip
    mu = np.zeros(p)
    for a in range(0, n, CHUNK): mu += rows_fn(a, min(n, a + CHUNK)).sum(0, dtype=np.float64)
    mu = (mu / n).astype(dtype); A = lambda a, b: rows_fn(a, b) - mu
    AQ = lambda Q: np.concatenate([A(a, min(n, a + CHUNK)) @ Q for a in range(0, n, CHUNK)])
    ATQ = lambda P: sum(A(a, min(n, a + CHUNK)).T @ P[a:a + CHUNK] for a in range(0, n, CHUNK))
    k = n_pc + n_oversamples; n_iter = 7 if n_pc < .1 * min(n, p) else 4
    Q = np.random.RandomState(seed).normal(size=(p, k)).astype(dtype, copy=False)
    lu = lambda M: linalg.lu(M, permute_l=True, check_finite=False)[0]
    for _ in range(n_iter): Q = lu(AQ(Q)); Q = lu(ATQ(Q))
    Q = linalg.qr(AQ(Q), mode="economic", check_finite=False)[0]
    Uhat, s, Vt = linalg.svd(ATQ(Q).T, full_matrices=False, lapack_driver="gesdd"); U = Q @ Uhat
    U, Vt = svd_flip(U[:, :n_pc], Vt[:n_pc], u_based_decision=False)
    return (U * s[:n_pc]).astype(np.float32)


def expr_embed(X, n_pc=100, norm="pearson", scale="z", seed=0):
    """floors.simple.expr_embed for norm 'pearson' (analytic Pearson residuals, theta 100, clip sqrt(n), float64);
    scale none (centred) or z (z-scored, clip 10)."""
    assert norm == "pearson", norm
    X = sp.csr_matrix(X, dtype=np.float32); n, p = X.shape; s = np.maximum(np.asarray(X.sum(1)).ravel(), 1)
    pg = np.asarray(X.sum(0)).ravel() / s.sum(); c = np.sqrt(n)

    def resid(a, b):
        Y = X[a:b].toarray().astype(np.float64); mu = s[a:b, None] * pg[None]; Y -= mu; Y /= np.sqrt(mu + mu ** 2 / 100 + 1e-12); np.clip(Y, -c, c, out=Y); return Y
    m = sum(resid(a, min(n, a + CHUNK)).sum(0) for a in range(0, n, CHUNK)) / n
    if scale == "z":
        sd = np.sqrt(sum(((resid(a, min(n, a + CHUNK)) - m) ** 2).sum(0) for a in range(0, n, CHUNK)) / n); sd = np.where(sd > 0, sd, 1)
        fn = lambda a, b: np.clip((resid(a, b) - m) / sd, -10, 10)
    else:
        fn = lambda a, b: resid(a, b) - m  # registered: Y -= mu before PCA
    return _rpca(fn, n, p, min(n_pc, p - 1), np.float64, seed)


def expr_pca(X, n_pc=50, seed=0):
    """floors.simple.expr_pca: log1p-CP10k over the given genes (float32, as lognorm returns it), genes z-scored
    (population sd), clip 10, randomized PCA."""
    from .simple import lognorm
    Y = lognorm(X).tocsr(); n, p = Y.shape; dt = Y.dtype; m = np.asarray(Y.astype(np.float64).mean(0)).ravel()
    sd = np.sqrt(sum(((Y[a:a + CHUNK].toarray().astype(np.float64) - m) ** 2).sum(0) for a in range(0, n, CHUNK)) / n)
    m, sd = m.astype(dt), np.where(sd > 0, sd, 1).astype(dt)
    return _rpca(lambda a, b: np.clip((Y[a:b].toarray() - m) / sd, -10, 10).astype(dt, copy=False), n, p, min(n_pc, p), dt, seed)


def knn_gauss(xy, ell_um, knn, conds):
    """(nb int32, w float32) of the registered k-NN Gaussian kernels (floors.f3.kernel, floors.simple.neighbour_mean):
    the knn nearest other cells (query knn + 1, first column dropped), w = exp(-d^2 / 2 ell^2) where every condition
    holds, else 0; conds = [(labels, same)]: neighbour label == own label (same) or != (not same). Queried in row
    chunks: same neighbours and weights, int32 indices."""
    from scipy.spatial import cKDTree
    xy = np.asarray(xy, float); n = len(xy); tree = cKDTree(xy); nb = np.empty((n, knn), np.int32); w = np.empty((n, knn), np.float32)
    for a in range(0, n, CHUNK):
        b = min(n, a + CHUNK); d, ii = tree.query(xy[a:b], knn + 1); d, ii = d[:, 1:], ii[:, 1:]; keep = np.ones(ii.shape, bool)
        for lab, same in conds: keep &= (lab[ii] == lab[a:b, None]) if same else (lab[ii] != lab[a:b, None])
        nb[a:b] = ii; w[a:b] = (np.exp(-d ** 2 / (2 * ell_um ** 2)) * keep).astype(np.float32)
    return nb, w


def ring(D, lo, hi, kmax):
    """floors.f3.ring as a sparse cells x cells 0/1 matrix (same selection: lo <= d < hi among the kmax nearest, same
    lineage, FOV and quadrant); returned as (W, None) so bleedsafe.kernel_mean(*ring, F, rows) works unchanged."""
    from scipy.spatial import cKDTree
    xy = D["xy"]; n = len(xy); tree = cKDTree(xy); ind, ptr = [], [0]
    for a in range(0, n, CHUNK):
        b = min(n, a + CHUNK); dd, ii = tree.query(xy[a:b], kmax, distance_upper_bound=hi, workers=4); ok = (ii < n) & (dd >= lo); ii = np.where(ok, ii, 0)
        me = np.arange(a, b)[:, None]
        ok &= (D["lin"][ii] == D["lin"][me]) & (D["fov"][ii] == D["fov"][me]) & (D["blocks"][ii] == D["blocks"][me])
        ind.append(ii[ok].astype(np.int32)); ptr.append(ok.sum(1))
    ind = np.concatenate(ind); indptr = np.concatenate([[0], np.cumsum(np.concatenate(ptr[1:]))]).astype(np.int64)
    return sp.csr_matrix((np.ones(len(ind), np.float32), ind, indptr), shape=(n, n)), None


def kernel_mean_sparse(W, F, rows):
    out = np.zeros((len(rows), F.shape[1]), np.float32)
    for a in range(0, len(rows), CHUNK):
        Wr = W[rows[a:a + CHUNK]]; out[a:a + CHUNK] = (Wr @ F) / np.maximum(np.asarray(Wr.sum(1)), 1e-9)
    return out
