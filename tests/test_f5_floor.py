"""F5 floors on a planted toy: 4 fine labels (2 epithelial, 1 T, 1 macrophage) from distinct profiles; the query has
other depth (x0.3) and a per-gene platform factor (log-normal sd .3). Both floors must recover the coarse classes."""

import shutil, numpy as np, pytest
from spatialfloors.floors import f5
from spatialfloors.metrics.f5 import macro_f1_gated

FINE = np.array(["EpiT", "Epi", "TCD4", "Macro"])


def _toy(seed=0, G=120, n=2000):
    rng = np.random.default_rng(seed)
    prof = rng.gamma(0.5, 1, (4, G)) + 0.05
    prof[0, :10] *= 8
    prof[1, :10] *= 6
    prof[1, 10:15] *= 4
    prof[2, 20:30] *= 8
    prof[3, 30:40] *= 8
    yr, yq = rng.integers(0, 4, n), rng.integers(0, 4, n)
    plat = np.exp(rng.normal(0, 0.3, G))
    Xr = rng.poisson(prof[yr] * 3.0)
    Xq = rng.poisson(prof[yq] * plat * 0.9)
    return Xr.astype(np.float32), FINE[yr], Xq.astype(np.float32), f5.coarse(FINE[yq]), prof


def test_coarse_map():
    assert list(f5.coarse(["EpiT", "Plasma", "Mono", "Mast", "Fibro"])) == [
        "epithelial",
        "immune_nonmyeloid",
        "myeloid",
        "other",
        "other",
    ]


@pytest.mark.parametrize("scaling", ["own", "ref"])
def test_logreg_recovers_planted_types(scaling):
    Xr, yr, Xq, cq, _ = _toy()
    p = f5.coarse(f5.logreg(Xr, yr, Xq, 1.0, scaling))
    assert macro_f1_gated(cq, p, ("epithelial", "immune_nonmyeloid", "myeloid"))[0] > 0.9


@pytest.mark.skipif(shutil.which("pixi") is None, reason="needs pixi env r")
def test_insitutype_recovers_planted_types(tmp_path):
    Xr, yr, Xq, cq, prof = _toy()
    lab, pr = f5.insitutype(Xq, prof.T, list(FINE), workdir=tmp_path)
    assert (
        macro_f1_gated(cq, f5.coarse(lab), ("epithelial", "immune_nonmyeloid", "myeloid"))[0] > 0.9
        and np.isfinite(pr).all()
    )


def test_heldout_auroc_planted():
    """markers carried by the true class: correct labels -> high AUROC; shuffled labels -> ~.5; tiny class -> .5."""
    import scipy.sparse as sp
    from spatialfloors.metrics.f5 import heldout_auroc

    rng = np.random.default_rng(0)
    n = 3000
    y = rng.integers(0, 3, n)
    cls = np.array(["epithelial", "immune_nonmyeloid", "myeloid"])
    lam = np.full((n, 30), 1.0)
    lam[np.arange(n)[:, None], y[:, None] * 10 + np.arange(10)] = 6.0
    X = sp.csr_matrix(rng.poisson(lam).astype(float))
    vn = [f"g{i}" for i in range(30)]
    mk = {c: vn[i * 10 : i * 10 + 5] for i, c in enumerate(cls)}
    assert (
        heldout_auroc(cls[y], X, vn, mk)[0] > 0.95
        and abs(heldout_auroc(cls[rng.permutation(y)], X, vn, mk)[0] - 0.5) < 0.05
    )
    p = cls[y].copy()
    p[np.flatnonzero(y == 2)[10:]] = "other"
    assert heldout_auroc(p, X, vn, mk)[1]["myeloid"]["auroc"] == 0.5
