"""F1 floor v2: on Poisson cells whose depth varies 10x
independently of type, Pearson residuals encode less depth than log1p-CP10k  while keeping the types."""

import numpy as np
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import KMeans
from spatialfloors.floors.simple import expr_embed


def _cells(n=3000, g=200, k=4, seed=0):
    rng = np.random.default_rng(seed)
    prof = rng.gamma(0.5, 1, (k, g))
    prof /= prof.sum(1, keepdims=True)
    t = rng.integers(0, k, n)
    depth = np.exp(rng.uniform(np.log(20), np.log(200), n))
    return rng.poisson(prof[t] * depth[:, None]).astype(np.float32), t, np.log(depth)


def _leak(E, y):
    return 1 - np.mean((cross_val_predict(RidgeCV(), E, y, cv=5) - y) ** 2) / np.var(y)


def test_pearson_leaks_less_depth_and_keeps_types():
    X, t, ld = _cells()
    Ep, El = expr_embed(X, 20, "pearson"), expr_embed(X, 20, "logcp10k")
    assert _leak(Ep, ld) < _leak(El, ld) - 0.1, (_leak(Ep, ld), _leak(El, ld))
    assert adjusted_rand_score(t, KMeans(4, n_init=5, random_state=0).fit_predict(Ep)) > 0.8
