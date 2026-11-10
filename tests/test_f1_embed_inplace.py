"""Pin: the in-place, row-block expr_embed equals the original dense form bit for bit."""

import numpy as np, scipy.sparse as sp


def _dense(X, n_pc, norm, scale, seed=0):  # dense reference implementation
    from sklearn.decomposition import PCA

    X = sp.csr_matrix(X, dtype=np.float32)
    s = np.maximum(np.asarray(X.sum(1)).ravel(), 1)
    if norm == "pearson":
        D = X.toarray().astype(np.float64)
        mu = s[:, None] * (np.asarray(X.sum(0)).ravel() / s.sum())[None]
        Y = np.clip((D - mu) / np.sqrt(mu + mu**2 / 100 + 1e-12), -np.sqrt(len(D)), np.sqrt(len(D)))
    else:
        Y = (sp.diags(1e4 / s) @ X if norm != "log1p" else X).toarray()
        Y = np.sqrt(Y) if norm == "sqrtcp10k" else np.log1p(Y)
    mu, sd = Y.mean(0), Y.std(0)
    Z = np.clip((Y - mu) / np.where(sd > 0, sd, 1), -10, 10) if scale == "z" else Y - mu
    return (
        PCA(min(n_pc, Z.shape[1] - 1), svd_solver="randomized", random_state=seed).fit_transform(Z).astype(np.float32)
    )


def test_inplace_equals_dense():
    from spatialfloors.floors.simple import expr_embed

    rng = np.random.default_rng(0)
    X = sp.csr_matrix(rng.poisson(rng.gamma(0.3, 2, (1, 300)), (45000, 300)).astype(np.float32))
    for norm, scale in (("pearson", "z"), ("pearson", "none"), ("logcp10k", "z"), ("sqrtcp10k", "none")):
        assert np.array_equal(expr_embed(X, 30, norm, scale), _dense(X, 30, norm, scale)), (norm, scale)
