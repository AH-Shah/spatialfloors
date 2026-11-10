"""F3 floor on planted toys (spatialfloors.sim.niche_toy): leak tests and planted-signal checks."""

import numpy as np, pytest
from spatialfloors.sim import niche_toy
from spatialfloors.floors import f3
from spatialfloors.metrics import f3 as m3


def _D(case, seed=0):
    C, xy, fov, ids = niche_toy(case, seed=seed)
    return f3.prepare(C, xy, fov, ids=ids)


@pytest.mark.parametrize("case", ["null", "bleed"])
def test_no_signal_without_coupling(case):
    """bleed only (receiver RNA in neighbours) or nothing: the floor finds no niche signal."""
    s = f3.score_floor(_D(case), 25.0, 30)
    assert abs(f3.unit_score(s)) < 0.02 and max(abs(v[2][0]) for v in s.values()) < 0.1, s


def test_planted_coupling_found_by_context_not_mix():
    s = f3.score_floor(_D("couple"), 25.0, 30)
    t = f3.score_floor(_D("couple"), 25.0, 30, features="T")
    assert s[0][2][0] > 0.15 and abs(t[0][2][0]) < 0.05, (s[0], t[0])


def test_nonlinear_coupling_not_linear():
    """type-1 state x type-2 presence (zero-mean product over receivers): out of the ridge floor's reach (a P3 control
    for
    the learned models)."""
    s = f3.score_floor(_D("nonlinear"), 25.0, 30)
    assert abs(s[0][2][0]) < 0.05, s[0]


def test_test_targets_do_not_move_test_predictions():
    """Replacing the held-out quadrant's targets by noise leaves its predictions unchanged (fits see training rows
    only)."""
    D = _D("couple")
    nb, w = f3.kernel(D, 25.0)
    T = f3.id_mix(D, nb, w)
    rows, Y, dep = f3.target(D, 0)
    ctx = f3.context_fn(D, 0, rows, nb, w, 30)
    X = lambda q, tr: np.c_[T[rows], ctx(q, tr)]
    b = D["blocks"][rows]
    preds = []

    def fp(Xtr, Ytr, Xte):
        p = f3.ridge_fit_predict(Xtr, Ytr, Xte)
        preds.append(p)
        return p

    for Yv in (Y, np.where((b == 0)[:, None], np.random.default_rng(0).normal(size=Y.shape), Y)):
        m3.block_cv(X, Yv, D["xy"][rows], fp, f3.RF_UM, blocks=b, split=D["split"], covars=dep)
    n = len(preds) // 2
    assert n >= 3 and np.allclose(preds[0], preds[n]) and not np.allclose(preds[1], preds[n + 1])


def test_kernel_stays_in_fov_quadrant_other_lineage():
    D = _D("null")
    nb, w = f3.kernel(D, 25.0)
    i, k = np.nonzero(w > 0)
    j = nb[i, k]
    assert (
        (D["blocks"][i] == D["blocks"][j]).all()
        and (D["fov"][i] == D["fov"][j]).all()
        and (D["lin"][i] != D["lin"][j]).all()
    )
