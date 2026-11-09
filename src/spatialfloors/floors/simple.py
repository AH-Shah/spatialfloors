"""Floors F1 / F2 (PLAN §3): Expr PCA (scaled) and Expr PCA + all-type neighbour mean (block-scaled), k-means.
Grids (dev crop only, PLAN §6): n_pc {30, 50} x kernel {10, 25 um}. Input X = counts (csr) on the unit's panel genes."""
import numpy as np
import scipy.sparse as sp
from . import bleedsafe as bs


def lognorm(X):
    """log1p-CP10k over the given genes."""
    X = sp.csr_matrix(X, dtype=np.float32); s = np.asarray(X.sum(1)).ravel()
    Y = sp.csr_matrix(sp.diags(1e4 / np.maximum(s, 1)) @ X); Y.data = np.log1p(Y.data); return Y


def expr_pca(X, n_pc=50, seed=0):
    """F1 floor: log1p-CP10k, genes z-scored (clip 10), PCA n_pc (randomized SVD, seed)."""
    from . import lean
    if lean.on(): return lean.expr_pca(X, n_pc, seed)  # whole-slide arm only (SF_LEAN=1)
    from sklearn.decomposition import PCA
    Y = lognorm(X).toarray(); mu, sd = Y.mean(0), Y.std(0); Z = np.clip((Y - mu) / np.where(sd > 0, sd, 1), -10, 10)
    return PCA(n_pc, svd_solver="randomized", random_state=seed).fit_transform(Z).astype(np.float32)


def bleedsafe_context(X, xy_um, fov, ell_um=25., seed=0):
    """Bleed-safe neighbour mean (bleedsafe.ctx_bleedsafe): Leiden ids on Expr50 -> lineages -> local foreign mask ->
    Gaussian kernel (ell in um, 200 NN, same FOV, other lineages) mean of the SVD-50 of the genes the receiver's id does not dominate."""
    from scipy.spatial import cKDTree
    x = lognorm(X); fov = np.asarray(fov).astype(str); ids = bs.leiden_ids(bs.expr50(x, seed), seed=seed)
    ids = np.unique(ids, return_inverse=True)[1]; K = int(ids.max()) + 1
    mu_c = bs.id_means(x, ids, K); lin_c = bs.lineages(mu_c); nb15 = cKDTree(xy_um).query(xy_um, 16)[1][:, 1:]
    mask, _ = bs.local_foreign_mask(mu_c, ids, lin_c, nb15, fov)
    spacing = float(np.median(cKDTree(xy_um).query(xy_um, 2)[0][:, 1]))  # bs.neighbourhood measures ell in median NN spacings
    nb, w = bs.neighbourhood(xy_um, fov, ids, lin_c, ell=ell_um / spacing)
    return bs.ctx_bleedsafe(x, ids, fov, mask, nb, w, seed=seed), dict(n_ids=K, n_lineages=int(len(np.unique(lin_c))), spacing_um=spacing)


def shared_context(X, xy_um, fov, E, ell_um=25., seed=0):
    """Bleed-safe neighbour mean in ONE coordinate system: Gaussian kernel (ell um, 200 NN, same FOV, other lineages;
    bleedsafe.neighbourhood) mean of the crop's Expr PCA E. Same-lineage neighbours are excluded, so a cell's own type
    (and its bleed into same-type neighbours) does not enter; ctx_bleedsafe instead uses a per-receiver-id SVD, whose
    columns are not comparable across ids (P2 domain toy: k-means on it finds no domains)."""
    from scipy.spatial import cKDTree
    x = lognorm(X); fov = np.asarray(fov).astype(str); ids = bs.leiden_ids(bs.expr50(x, seed), seed=seed)
    ids = np.unique(ids, return_inverse=True)[1]; lin_c = bs.lineages(bs.id_means(x, ids, int(ids.max()) + 1))
    spacing = float(np.median(cKDTree(xy_um).query(xy_um, 2)[0][:, 1])); nb, w = bs.neighbourhood(xy_um, fov, ids, lin_c, ell=ell_um / spacing)
    return bs.kernel_mean(nb, w, E), dict(n_ids=int(ids.max()) + 1, n_lineages=int(len(np.unique(lin_c))), spacing_um=spacing)


def expr_embed(X, n_pc=100, norm="pearson", scale="z", seed=0):
    """F1 floor v2 (2026-10-29, PLAN §14, HYPOTHESES MP12): analytic Pearson residuals (theta 100, clip sqrt(n)), genes
    z-scored (clip 10), PCA n_pc 100. The v1 floor (expr_pca: log1p-CP10k, n_pc <= 50) leaked depth (CV R2 .84-.93) and
    lost rare Allen types (kNN recall .75 vs scVI .89); the landscape score rose with n_pc to 100.
    norm: pearson | logcp10k | log1p (raw counts) | sqrtcp10k; scale: z | none (centred)."""
    from sklearn.decomposition import PCA
    from . import lean
    if lean.on() and norm == "pearson": return lean.expr_embed(X, n_pc, norm, scale, seed)  # whole-slide arm only (SF_LEAN=1)
    X = sp.csr_matrix(X, dtype=np.float32); s = np.maximum(np.asarray(X.sum(1)).ravel(), 1)
    if norm == "pearson":  # in place, row blocks: the same elementwise operations as the dense form, ~1/3 of the memory (audit 2026-10-30: whole slides)
        Y = X.toarray().astype(np.float64); p = np.asarray(X.sum(0)).ravel() / s.sum(); c = np.sqrt(len(Y))
        for i in range(0, len(Y), 20000):
            mu = s[i:i + 20000, None] * p[None]; Y[i:i + 20000] -= mu; Y[i:i + 20000] /= np.sqrt(mu + mu ** 2 / 100 + 1e-12)
        np.clip(Y, -c, c, out=Y)
    else:
        Y = (sp.diags(1e4 / s) @ X if norm != "log1p" else X).toarray()
        Y = np.sqrt(Y) if norm == "sqrtcp10k" else np.log1p(Y)
    mu, sd = Y.mean(0), Y.std(0); Y -= mu
    if scale == "z": Y /= np.where(sd > 0, sd, 1); np.clip(Y, -10, 10, out=Y)
    return PCA(min(n_pc, Y.shape[1] - 1), svd_solver="randomized", random_state=seed, copy=False).fit_transform(Y).astype(np.float32)


def neighbour_mean(xy_um, fov, E, ell_um=25., knn=200):
    """BANKSY-class neighbour mean: Gaussian kernel (ell um) over the 200 NN in the same FOV, self excluded, all types."""
    from scipy.spatial import cKDTree
    from . import lean
    fov = np.asarray(fov).astype(str)
    if lean.on(): return bs.kernel_mean(*lean.knn_gauss(xy_um, ell_um, knn, [(np.unique(fov, return_inverse=True)[1].astype(np.int32), True)]), E), {}  # whole-slide arm: chunked query
    d, nb = cKDTree(xy_um).query(xy_um, knn + 1); d, nb = d[:, 1:], nb[:, 1:]
    w = (np.exp(-d ** 2 / (2 * ell_um ** 2)) * (fov[nb] == fov[:, None])).astype(np.float32); return bs.kernel_mean(nb, w, E), {}


def domain_features(X, xy_um, fov, n_pc=50, ell_um=25., seed=0, ctx="all", scale="block", lam=.5):
    """F2 floor features: [Expr PCA, context]. Locked floor (decision 2026-10-28, row 347): ctx "all", scale "block",
    lam .5. Diagnostics: ctx "bleedsafe" (bleedsafe_context, per-id spaces) or "shared"
    (shared_context) or "all" (neighbour_mean). scale "dim": every column z-scored, block weighted 1/sqrt(dim); "block": columns centred, block
    divided by its total SD (keeps the PCA variance order, so noise dimensions are not inflated). Block weights
    sqrt(1 - lam), sqrt(lam) (BANKSY's lambda: .5 = equal; BANKSY's domain-segmentation default .8)."""
    E = expr_pca(X, n_pc, seed)
    C, info = (bleedsafe_context(X, xy_um, fov, ell_um, seed) if ctx == "bleedsafe" else shared_context(X, xy_um, fov, E, ell_um, seed)
               if ctx == "shared" else neighbour_mean(xy_um, fov, E, ell_um))
    if scale == "dim": z = lambda M: (M - M.mean(0)) / np.maximum(M.std(0), 1e-6) / np.sqrt(M.shape[1])
    else: z = lambda M: (M - M.mean(0)) / max(float(np.sqrt(((M - M.mean(0)) ** 2).sum(1).mean())), 1e-6)
    return np.hstack([np.sqrt(1 - lam) * z(E), np.sqrt(lam) * z(C)]).astype(np.float32), info


def grid_smooth(xy_um, F, sigma_um, bin_um=10.):
    """Gaussian neighbour mean at a true scale (PLAN §19 operator fix): features binned on a bin_um grid of section
    coordinates, Gaussian-filtered (sigma_um) and divided by the filtered cell counts; FOV edges are not cuts (one section
    = one coordinate system). neighbour_mean instead averages the 200 nearest cells of the same FOV, so every ell above
    ~100 um (CosMx) / ~165 um (Allen) is a flat mean over that radius."""
    from scipy.ndimage import gaussian_filter
    ij = np.floor((xy_um - xy_um.min(0)) / bin_um).astype(int); sh = tuple(ij.max(0) + 1); F = np.asarray(F, np.float32)
    S = np.zeros(sh + (F.shape[1],), np.float32); C = np.zeros(sh, np.float32); np.add.at(S, (ij[:, 0], ij[:, 1]), F); np.add.at(C, (ij[:, 0], ij[:, 1]), 1.)
    s = sigma_um / bin_um; S = gaussian_filter(S, (s, s, 0), mode="constant"); C = gaussian_filter(C, s, mode="constant")
    return S[ij[:, 0], ij[:, 1]] / np.maximum(C[ij[:, 0], ij[:, 1]], 1e-9)[:, None]


def v3_config(name):
    """Parse an F2 v3 config name (scripts/v3_f2.py): 'fixed400' | 'ms<s1>_<s2>..[+c<s>[_<s>]][_lam<l>][_pc<n>]' ->
    domain_features_v3 keyword arguments."""
    import re
    if name == "fixed400": return dict(sigmas=(400.,), n_out=0)
    m = re.match(r"^ms([\d_]+?)(?:\+c([\d_]+?))?(?:_lam([.\d]+))?(?:_pc(\d+))?$", name); assert m, f"bad F2 v3 config {name}"
    kw = dict(sigmas=tuple(float(s) for s in m[1].split("_")))
    if m[2]: c = tuple(float(s) for s in m[2].split("_")); kw["comp_sigma"] = c[0] if len(c) == 1 else c
    if m[3]: kw["lam"] = float(m[3])
    if m[4]: kw["n_out"] = int(m[4])
    return kw


def domain_features_v3(X, xy_um, n_pc=50, sigmas=(25., 100., 400.), comp_sigma=None, lam=.8, n_out=60, seed=0):
    """F2 floor v3 features (PLAN §19): [Expr PCA, grid-smoothed Expr PCA at each sigma (true scale, across FOVs),
    optionally the grid-smoothed Leiden-id composition at comp_sigma (one sigma or several)]; blocks centred and divided by their total SD,
    weights sqrt(1 - lam) for expression and sqrt(lam / n_context_blocks) for each context block; PCA to n_out (GMM full on
    <= 100 dims, as the v2 lock). One sigma without comp_sigma = the v2 design with the operator fixed."""
    from sklearn.decomposition import PCA
    E = expr_pca(X, n_pc, seed); z = lambda M: (M - M.mean(0)) / max(float(np.sqrt(((M - M.mean(0)) ** 2).sum(1).mean())), 1e-6)
    ctx = [grid_smooth(xy_um, E, s) for s in sigmas]
    if comp_sigma is not None:  # one sigma or several (each a block)
        ids = bs.leiden_ids(bs.expr50(lognorm(X), seed), seed=seed); ids = np.unique(ids, return_inverse=True)[1]; OH = np.eye(int(ids.max()) + 1, dtype=np.float32)[ids]
        ctx += [grid_smooth(xy_um, OH, s) for s in np.atleast_1d(comp_sigma)]
    F = np.hstack([np.sqrt(1 - lam) * z(E)] + [np.sqrt(lam / len(ctx)) * z(c) for c in ctx]).astype(np.float32)
    return (PCA(min(n_out, F.shape[1]), random_state=seed).fit_transform(F).astype(np.float32) if n_out else F), dict(sigmas=list(sigmas), comp_sigma=comp_sigma)
