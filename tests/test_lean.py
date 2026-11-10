"""Whole-slide arm (floors/lean.py, SF_LEAN=1): the lean paths reproduce the registered ones on synthetic data."""

import numpy as np
import scipy.sparse as sp


def _counts(n=30000, g=300, seed=0):
    rng = np.random.default_rng(seed)
    z = rng.integers(0, 6, n)
    mu = rng.gamma(0.5, 1.0, (6, g))
    return sp.csr_matrix(rng.poisson(mu[z] * rng.uniform(0.5, 2, (n, 1))).astype(np.float32))  # > 1 lean row chunk


def test_lean_pca_equals_registered(monkeypatch):
    """The registered randomized PCA (sklearn) and its chunked replica: same scores up to summation order (float64:
    expr_embed) or float32 rounding (expr_pca), signs included."""
    from spatialfloors.floors import simple

    X = _counts()
    reg = {
        "embed_none": simple.expr_embed(X, 20, "pearson", "none"),
        "embed_z": simple.expr_embed(X, 20, "pearson", "z"),
        "pca": simple.expr_pca(X, 20),
    }
    monkeypatch.setenv("SF_LEAN", "1")
    lean = {
        "embed_none": simple.expr_embed(X, 20, "pearson", "none"),
        "embed_z": simple.expr_embed(X, 20, "pearson", "z"),
        "pca": simple.expr_pca(X, 20),
    }
    for k in ("embed_none", "embed_z"):
        err = np.abs(reg[k] - lean[k]).max() / np.abs(reg[k]).max()
        assert err < 1e-6, (k, err)
    # float32 (as registered): the 6 signal components (6 planted groups) to 1e-3; the near-degenerate noise
    # components (variance gaps < 1%) carry float32 summation-order differences, to 1e-2
    err = np.abs(reg["pca"] - lean["pca"]).max(0) / np.abs(reg["pca"]).max()
    assert err[:6].max() < 1e-3 and err.max() < 1e-2, err


def test_lean_ring_equals_dense(monkeypatch):
    from spatialfloors.floors import f3, bleedsafe as bs

    rng = np.random.default_rng(1)
    n = 4000
    xy = rng.uniform(0, 1500, (n, 2))
    D = dict(
        xy=xy,
        lin=rng.integers(0, 3, n),
        fov=(xy[:, 0] // 500).astype(int).astype(str),
        blocks=(xy[:, 1] > 750).astype(int),
    )
    F = rng.normal(size=(n, 5)).astype(np.float32)
    rows = np.sort(rng.choice(n, 1500, replace=False))
    dense = bs.kernel_mean(*f3.ring(D, 50.0, 400.0, 300), F, rows)
    monkeypatch.setenv("SF_LEAN", "1")
    W = f3.ring(D, 50.0, 400.0, 300)
    assert sp.issparse(W[0])
    assert np.allclose(bs.kernel_mean(*W, F, rows), dense, atol=1e-5)


def test_lean_knn_kernels_equal(monkeypatch):
    from spatialfloors.floors import f3, simple

    rng = np.random.default_rng(2)
    n = 45000
    xy = rng.uniform(0, 3000, (n, 2))
    E = rng.normal(size=(n, 4)).astype(np.float32)
    D = dict(
        xy=xy,
        lin=rng.integers(0, 3, n),
        fov=(xy[:, 0] // 500).astype(int).astype(str),
        blocks=(xy[:, 1] > 1500).astype(int),
    )
    reg_k, reg_m = f3.kernel(D, 25.0), simple.neighbour_mean(xy, D["fov"], E, 400.0)[0]
    monkeypatch.setenv("SF_LEAN", "1")
    D["fov"] = np.unique(D["fov"], return_inverse=True)[1].astype(np.int32)
    lk, lm = f3.kernel(D, 25.0), simple.neighbour_mean(xy, (xy[:, 0] // 500).astype(int).astype(str), E, 400.0)[0]
    assert np.array_equal(reg_k[0], lk[0]) and np.array_equal(reg_k[1], lk[1]) and np.array_equal(reg_m, lm)
