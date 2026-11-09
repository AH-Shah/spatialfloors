"""Bleed-safe niche programs (CALCULUS §60): label-free cross-type state coupling at ~40 µm, beyond the neighbourhood type mix.

Pipeline (all label-free; x = log1p-CP10k cells x genes, csr):
  ids      Leiden on unit rows of Expr50 (or any given identity)
  lineage  ids joined by centroid correlation >= .7 (sister clusters of one cell type share a lineage)
  mask     local foreign bleed mask: id r keeps gene g if it expresses g at least as much as its actual foreign-lineage neighbours
  Z        own state per id: PCA-10 of the masked genes, renormalised on those genes (foreign bleed cannot enter via CP10k), FOV-centred,
           depth-regressed, genes z-scored within the id
  w, T     Gaussian kernel (l 5 cell spacings, 200 NN, same FOV, other lineages); T = kernel id mix
  Sc       per sender id: kernel mean of its Z, denominator + tau_w (rare far senders shrink instead of dominating)
  ctx_bs   fair neighbour-expression baseline (BANKSY/Context class): per receiver id, SVD-50 of the genes it does not dominate,
           renormalised on those genes, kernel-averaged
  programs per receiver id: shrinkage CCA of res_T(Z) vs res_T(Sc); held-out rho over 3 FOV folds; baselines scored the same way
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

ELL, KNN, DZ, NAX, REG, TAU_W, LIN_CORR = 5.0, 200, 10, 3, .1, 1.0, .7


# ------------------------------------------------------------------ helpers
def renorm_log(x):
    """log1p-CP10k restricted to a gene subset, renormalised on that subset (rows sum to 1e4 before log1p)."""
    y = x.copy().tocsr(); y.data = np.expm1(y.data); s = np.asarray(y.sum(1)).ravel()
    y = sp.csr_matrix(sp.diags(1e4 / np.maximum(s, 1e-9)) @ y); y.data = np.log1p(y.data); return y


def group_center(M, g):
    return M - pd.DataFrame(M).groupby(np.asarray(g)).transform("mean").values


def cca(X, Y, k=NAX, reg=REG):
    mx, my = X.mean(0), Y.mean(0); X, Y = X - mx, Y - my; n = len(X)
    def isq(C):
        C = (1 - reg) * C + reg * np.trace(C) / len(C) * np.eye(len(C)); e, V = np.linalg.eigh(C); return (V / np.sqrt(np.maximum(e, 1e-12))) @ V.T
    Wx, Wy = isq(X.T @ X / n), isq(Y.T @ Y / n); U, S, Vt = np.linalg.svd(Wx @ (X.T @ Y / n) @ Wy, full_matrices=False)
    return mx, my, Wx @ U[:, :k], Wy @ Vt.T[:, :k]


def folds(groups, k=3, seed=0):
    u = np.random.default_rng(seed).permutation(np.unique(groups)); f = {g: i % k for i, g in enumerate(u)}; return np.array([f[g] for g in groups])


def heldout_rho(X, Y, fo, k=1):
    """mean over FOV folds of held-out canonical correlations (first k)."""
    out = []
    for q in np.unique(fo):
        mx, my, A, B = cca(X[fo != q], Y[fo != q], k)
        out.append([np.corrcoef((X[fo == q] - mx) @ A[:, j], (Y[fo == q] - my) @ B[:, j])[0, 1] for j in range(k)])
    return np.mean(out, 0)


# ------------------------------------------------------------------ identity, lineage, mask
def expr50(x, seed=0):
    from sklearn.decomposition import TruncatedSVD
    sub = np.random.default_rng(seed).choice(x.shape[0], min(x.shape[0], 100000), replace=False)
    return TruncatedSVD(50, random_state=seed).fit(x[sub]).transform(x).astype(np.float32)


def leiden_ids(E, k=15, res=1.0, seed=0):
    """Leiden (RB configuration, resolution res) on a k-NN graph of the unit rows of E."""
    from sklearn.neighbors import NearestNeighbors
    import igraph as ig, leidenalg as la
    U = (E / np.maximum(np.linalg.norm(E, axis=1, keepdims=True), 1e-12)).astype(np.float32); n = len(U)
    nb = NearestNeighbors(n_neighbors=k).fit(U).kneighbors(U, return_distance=False)
    g = ig.Graph(n=n, edges=np.c_[np.repeat(np.arange(n), k), nb.ravel()].tolist(), directed=False).simplify()
    return np.asarray(la.find_partition(g, la.RBConfigurationVertexPartition, resolution_parameter=res, seed=seed).membership)


def id_means(x, ids, K):
    S = sp.csr_matrix((np.ones(len(ids)), (ids, np.arange(len(ids)))), shape=(K, len(ids))); cnt = np.asarray(S.sum(1)).ravel()
    return (S @ x).toarray() / np.maximum(cnt, 1)[:, None]


def lineages(mu_c, corr=LIN_CORR):
    from scipy.sparse.csgraph import connected_components
    return connected_components(sp.csr_matrix(np.corrcoef(mu_c) >= corr), directed=False)[1]


def local_foreign_mask(mu_c, ids, lin_c, nb, fov, k=15, min_expr=.02, min_genes=15):
    """K x G mask: id r keeps g if mu_c[r,g] >= the neighbour-weighted mean of its foreign-lineage neighbours' ids and >= min_expr."""
    K, G = mu_c.shape; lin = lin_c[ids]; nk = nb[:, :k]; f = (lin[nk] != lin[:, None]) & (fov[nk] == fov[:, None])
    W = np.zeros((K, K)); np.add.at(W, (np.repeat(ids, k)[f.ravel()], ids[nk].ravel()[f.ravel()]), 1.0)
    mu_for = (W / np.maximum(W.sum(1, keepdims=True), 1e-9)) @ mu_c; mask = np.zeros((K, G), bool)
    for r in range(K):
        m = (mu_c[r] >= mu_for[r]) & (mu_c[r] >= min_expr) if W[r].sum() > 0 else mu_c[r] >= min_expr
        if m.sum() < min_genes: m = np.zeros(G, bool); m[np.argsort(-((mu_c[r] - mu_for[r]) * (mu_c[r] >= .005)))[:min_genes]] = True
        mask[r] = m
    return mask, mu_for


# ------------------------------------------------------------------ states and neighbourhood
def own_states(x, ids, fov, mask, dz=DZ, min_cells=200, seed=0):
    """Z (n x dz) per id: PCA of the masked genes renormalised on them, FOV-centred, regressed on log detected genes; loadings per id."""
    from sklearn.decomposition import PCA
    Z = np.zeros((len(ids), dz), np.float32); load = {}
    for r in np.unique(ids):
        idx = np.flatnonzero(ids == r); g = np.flatnonzero(mask[r])
        if len(idx) < min_cells or len(g) <= dz: continue
        xr = renorm_log(x[idx][:, g]); X = group_center(xr.toarray(), fov[idx])
        dep = group_center(np.log1p(np.diff(xr.indptr)).astype(float)[:, None], fov[idx]); X = X - dep @ np.linalg.lstsq(dep, X, rcond=None)[0]
        X = X / np.maximum(X.std(0), 1e-6)  # genes z-scored within the id: coherent state beats per-gene Poisson noise (NV1)
        pca = PCA(dz, random_state=seed).fit(X); Z[idx] = pca.transform(X); load[r] = (g, pca.components_)
    return Z, load


def neighbourhood(pos, fov, ids, lin_c, ell=ELL, knn=KNN):
    """200-NN, distances in median nearest-neighbour spacings; kernel w over other lineages in the same FOV."""
    from scipy.spatial import cKDTree
    d, nb = cKDTree(pos).query(pos, knn + 1); d, nb = d[:, 1:], nb[:, 1:]; d = d / np.median(d[:, 0]); lin = lin_c[ids]
    w = (np.exp(-d ** 2 / (2 * ell ** 2)) * (fov[nb] == fov[:, None]) * (lin[nb] != lin[:, None])).astype(np.float32)
    return nb, w


def id_mix(nb, w, ids, K):
    T = np.zeros((len(ids), K), np.float32); np.add.at(T, (np.repeat(np.arange(len(ids)), nb.shape[1]), ids[nb].ravel()), w.ravel())
    return T / np.maximum(T.sum(1, keepdims=True), 1e-9)


def sender_means(Z, ids, nb, w, K, tau_w=TAU_W, chunk=20000):
    """n x (K*dz): per sender id the kernel mean of its state, sum w z / (sum w + tau_w)."""
    n, dz = Z.shape; Sc = np.zeros((n, K * dz), np.float32)
    for a in range(0, n, chunk):
        sl = slice(a, a + chunk); Zn = Z[nb[sl]]; cn = ids[nb[sl]]
        for c in range(K):
            wc = w[sl] * (cn == c); Sc[sl, c * dz:(c + 1) * dz] = np.einsum("nk,nkd->nd", wc, Zn) / (wc.sum(1, keepdims=True) + tau_w)
    return Sc


def kernel_mean(nb, w, F, rows=None, chunk=20000):
    if sp.issparse(nb):  # lean.ring (whole-slide arm, SF_LEAN=1): nb is a cells x cells 0/1 matrix, w None
        from .lean import kernel_mean_sparse
        return kernel_mean_sparse(nb, F, np.arange(nb.shape[0]) if rows is None else rows)
    rows = np.arange(len(nb)) if rows is None else rows; out = np.zeros((len(rows), F.shape[1]), np.float32)
    for a in range(0, len(rows), chunk):
        r = rows[a:a + chunk]; out[a:a + chunk] = np.einsum("nk,nkd->nd", w[r], F[nb[r]]) / np.maximum(w[r].sum(1, keepdims=True), 1e-9)
    return out


def ctx_bleedsafe(x, ids, fov, mask, nb, w, dz=50, seed=0):
    """fair neighbour-expression baseline: row i = kernel mean over its neighbours of the SVD-dz of the genes id(i) does not dominate."""
    from sklearn.decomposition import TruncatedSVD
    n = len(ids); out = np.zeros((n, dz), np.float32); sub = np.random.default_rng(seed).choice(n, min(n, 100000), replace=False)
    for r in np.unique(ids):
        g = np.flatnonzero(~mask[r]); rows = np.flatnonzero(ids == r); k = min(dz, len(g) - 1)  # fewer genes than dz: use what exists
        if k < 2: continue
        xg = renorm_log(x[:, g]); m1 = np.asarray(xg[sub].mean(0)).ravel(); sd = np.sqrt(np.maximum(np.asarray(xg[sub].multiply(xg[sub]).mean(0)).ravel() - m1 ** 2, 1e-6))
        xg = sp.csr_matrix(xg @ sp.diags(1 / sd)); E = np.zeros((n, dz), np.float32); E[:, :k] = group_center(TruncatedSVD(k, random_state=seed).fit(xg[sub]).transform(xg), fov)
        out[rows] = kernel_mean(nb, w, E, rows)
    return out


# ------------------------------------------------------------------ programs
def programs(Z, Sc, T, ids, fov, load, K, baselines=None, genes=None, min_cells=2000, nmax=40000, seed=0):
    """per receiver id: held-out rho of the programs (Sc) and of each baseline (n x d arrays), all beyond T; U, V program scores; read-out table."""
    dz = Z.shape[1]; rows, progs = [], []; U = np.zeros((len(Z), NAX), np.float32); V = np.zeros_like(U)
    for r in range(K):
        allr = np.flatnonzero(ids == r)
        if len(allr) < min_cells or r not in load: continue
        idx = allr if len(allr) <= nmax else np.sort(np.random.default_rng(seed).choice(allr, nmax, replace=False))
        C = np.c_[np.ones(len(idx)), T[idx]]; res = lambda M: M - C @ np.linalg.lstsq(C, M, rcond=None)[0]
        live = np.flatnonzero(Sc[idx].std(0) > 1e-6); X, Y = res(Z[idx]), res(Sc[idx][:, live]); fo = folds(fov[idx], seed=seed)
        rho = heldout_rho(X, Y, fo, NAX); row = dict(receiver=r, n=len(allr), rho1=rho[0], rho2=rho[1], rho3=rho[2])
        for name, Bm in (baselines or {}).items(): row[f"{name}_rho1"] = float(heldout_rho(X, res(Bm[idx]), fo, 1)[0])
        mx, my, A, B = cca(X, Y); Cr = np.c_[np.ones(len(allr)), T[allr]]
        U[allr] = (Z[allr] - Cr @ np.linalg.lstsq(C, Z[idx], rcond=None)[0] - mx) @ A; V[allr] = (Sc[allr][:, live] - Cr @ np.linalg.lstsq(C, Sc[idx][:, live], rcond=None)[0] - my) @ B
        g, comp = load[r]
        for k in range(NAX):
            share = {c: float(np.var(Y[:, cols] @ B[cols, k])) for c in range(K) for cols in [np.flatnonzero((live >= c * dz) & (live < (c + 1) * dz))] if len(cols)}
            tot = sum(share.values()) or 1; top = sorted(share, key=share.get, reverse=True)[:3]; gl = comp.T @ A[:, k]; o = np.argsort(-np.abs(gl))[:10]
            progs.append(dict(receiver=r, axis=k, rho=rho[k], senders=";".join(f"{c}:{share[c] / tot:.2f}" for c in top),
                              receiver_genes=";".join(f"{genes[g[i]] if genes is not None else g[i]}{gl[i]:+.2f}" for i in o)))
        rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(progs), U, V


def run(x, pos, fov, ids=None, genes=None, E=None, min_cells=2000, tau_w=TAU_W, seed=0, legacy_ctx=False):
    """end-to-end bleed-safe niche programs. Returns dict(scores, programs, U, V, ids, lineage, mask, Z, T, Sc)."""
    x = x.tocsr().astype(np.float32); fov = np.asarray(fov).astype(str)
    if ids is None: E = expr50(x, seed) if E is None else E; ids = leiden_ids(E, seed=seed)
    ids = np.unique(np.asarray(ids), return_inverse=True)[1]; K = int(ids.max()) + 1
    from scipy.spatial import cKDTree
    mu_c = id_means(x, ids, K); lin_c = lineages(mu_c); nb15 = cKDTree(pos).query(pos, 16)[1][:, 1:]
    mask, _ = local_foreign_mask(mu_c, ids, lin_c, nb15, fov)
    Z, load = own_states(x, ids, fov, mask, seed=seed); nb, w = neighbourhood(pos, fov, ids, lin_c); T = id_mix(nb, w, ids, K)
    Sc = sender_means(Z, ids, nb, w, K, tau_w); base = {"ctx_bs": ctx_bleedsafe(x, ids, fov, mask, nb, w, seed=seed)}
    if legacy_ctx: base["ctx_legacy"] = kernel_mean(nb, w, expr50(x, seed) if E is None else E)  # leak-prone (receiver bleed), diagnostic only
    scores, progs, U, V = programs(Z, Sc, T, ids, fov, load, K, base, genes, min_cells, seed=seed)
    return dict(scores=scores, programs=progs, U=U, V=V, ids=ids, lineage=lin_c[ids], mask=mask, Z=Z, T=T, Sc=Sc)


def save(out, directory, name, barcodes=None):
    from pathlib import Path
    d = Path(directory); d.mkdir(parents=True, exist_ok=True)
    out["scores"].to_csv(d / f"niche_{name}_scores.csv", index=False); out["programs"].to_csv(d / f"niche_{name}_programs.csv", index=False)
    np.savez_compressed(d / f"niche_{name}.npz", U=out["U"], V=out["V"], ids=out["ids"], lineage=out["lineage"],
                        barcodes=np.asarray(barcodes).astype(str) if barcodes is not None else np.arange(len(out["ids"])))
