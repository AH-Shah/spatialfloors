import numpy as np, scipy.sparse as sp, pytest
from spatialfloors.metrics import f1, f2, f3, f5
from spatialfloors import stats


def _groups(n=1500, G=200, K=5, depth=300, seed=0):
    rng = np.random.default_rng(seed)
    g = rng.integers(0, K, n)
    P = rng.dirichlet(np.full(G, 0.3), K)
    X = sp.csr_matrix(rng.poisson(P[g] * depth))
    return X, g, P


def test_count_split_sums_and_seed():
    X, *_ = _groups()
    a, b = f1.count_split(X, 0)
    assert (a + b != X).nnz == 0 and a.min() >= 0 and b.min() >= 0
    a2, _ = f1.count_split(X, 0)
    assert (a != a2).nnz == 0


def test_knn_self_excluded_even_with_duplicates():
    E = np.repeat(np.random.default_rng(1).normal(size=(50, 4)), 3, 0)
    nb = f1.knn_cosine(E, 5)
    assert not (nb == np.arange(len(E))[:, None]).any()


def test_f1_perfect_high_shuffled_zero():
    X, g, _ = _groups()
    a, b = f1.count_split(X, 0)
    genes = np.arange(X.shape[1])
    E_true = np.eye(5)[g] + 1e-3 * np.random.default_rng(0).normal(size=(len(g), 5))
    s_true = f1.knn_poisson_score(E_true, a, b, genes)
    s_shuf = [f1.knn_poisson_score(np.random.default_rng(i).normal(size=(len(g), 5)), a, b, genes) for i in (2, 3)]
    ch = f1.chance_level(E_true, a, b, genes, n_perm=2)
    # chance is not 0: k random neighbours give a noisier profile than pbar; it is stable and < 0
    assert (
        s_true > 0.3 and max(s_shuf) < 0 and abs(s_shuf[0] - s_shuf[1]) < 0.03 and abs(ch - np.mean(s_shuf)) < 0.03
    ), (s_true, s_shuf, ch)


def test_f1_perfect_tends_to_one_with_depth():
    X, g, _ = _groups(depth=5000)
    a, b = f1.count_split(X, 0)
    E = np.eye(5)[g] + 1e-3 * np.random.default_rng(0).normal(size=(len(g), 5))
    # lam = .5 caps a perfect embedding below 1 (half the profile is pbar); lam -> 0 recovers the -> 1 property
    s5, s0 = (f1.knn_poisson_score(E, a, b, np.arange(X.shape[1]), lam=l) for l in (0.5, 1e-3))
    X3, *_ = _groups()
    a3, b3 = f1.count_split(X3, 0)
    assert s0 > 0.9 and s5 > f1.knn_poisson_score(E, a3, b3, np.arange(X.shape[1])), s5


def test_f1_finite_with_zero_genes_and_cells():
    X, g, _ = _groups()
    X = X.tolil()
    X[:, :5] = 0
    X[:10] = 0
    X = X.tocsr()
    a, b = f1.count_split(X, 0)
    s = f1.knn_poisson_score(np.eye(5)[g] + 0.01, a, b, np.arange(X.shape[1]))
    assert np.isfinite(s)


def test_f2_perfect_and_random():
    r = np.repeat([0, 1, 2], 200)
    assert f2.ari_at_true_k(r, cluster_fn=lambda k, s: r)[0] == 1
    rng = np.random.default_rng(0)
    m, _ = f2.ari_at_true_k(r, cluster_fn=lambda k, s: np.random.default_rng(s).integers(0, k, len(r)))
    assert abs(m) < 0.02
    E = np.eye(3)[r] + 0.01 * rng.normal(size=(len(r), 3))
    assert f2.ari_at_true_k(r, E=E, seeds=range(2))[0] > 0.99


def test_f3_no_fit_on_test_quadrant_and_buffer():
    rng = np.random.default_rng(0)
    xy = rng.uniform(0, 100, (4000, 2))
    X = rng.normal(size=(4000, 12))
    Y = X @ rng.normal(size=(12, 20)) + 0.1 * rng.normal(size=(4000, 20))
    blocks, split = f3.quadrant_blocks(xy)
    buf = f3.edge_buffer(xy, split, 5)
    d = np.minimum(abs(xy[:, 0] - split[0]), abs(xy[:, 1] - split[1]))
    assert buf[d < 10].all() and not buf[d >= 10].any()
    log = f3.FitLog()
    seen = []

    def fp(Xtr, Ytr, Xte):
        from sklearn.linear_model import Ridge

        seen.append(len(Xtr))
        return Ridge(1.0).fit(Xtr, Ytr).predict(Xte)

    s, folds = f3.block_cv(X, Y, xy, fp, 5, fitlog=log)
    assert s > 0.9 and len(folds) == 4
    assert all(n < 3200 for n in seen)  # every fit excludes one quadrant

    # leak detector: a fit on all rows must trip the assertion
    class Leaky(f3.FitLog):
        def __call__(self, idx):
            self.rows.append(np.arange(len(xy)))
            return idx

    with pytest.raises(AssertionError):
        f3.block_cv(X, Y, xy, fp, 5, fitlog=Leaky())


def test_f3_graph_within_block():
    xy = np.array([[0, 0], [1, 1], [9, 9], [10, 10.0]])
    b, _ = f3.quadrant_blocks(xy)
    W = sp.csr_matrix(np.ones((4, 4)))
    R = f3.restrict_graph_to_blocks(W, b).toarray()
    assert (R[b[:, None] != b[None, :]] == 0).all()


def test_f5_gate_and_f1():
    rng = np.random.default_rng(0)
    n = 3000
    cls = rng.integers(0, 3, n)
    panck = np.where(cls == 0, 500, 5) * rng.lognormal(0, 0.3, n)
    cd45 = np.where(cls > 0, 500, 5) * rng.lognormal(0, 0.3, n)
    my = np.where(cls == 2, 500, 5) * rng.lognormal(0, 0.3, n)
    lab, thr = f5.if_gate(panck, cd45, my, q=(0.3, 0.3), myelo_q=(0.6, 0.7))
    truth = np.array(f5.CLASSES)[cls]
    m = lab != ""
    assert (lab[m] == truth[m]).mean() > 0.98
    assert f5.macro_f1_gated(lab, lab)[0] == 1


def test_stats_intercept_is_mean():
    D = np.array([0.1, 0.2, -0.05, 0.3, 0.0, 0.15])
    s = np.array(["a", "a", "a", "b", "b", "b"])
    est, se, df = stats.centred_intercept(D, s)
    assert np.isclose(est, D.mean()) and df == 4


def test_stats_one_stratum():
    D = np.array([0.1, 0.2, -0.05, 0.3])
    est, se, df = stats.centred_intercept(D, np.array(["a"] * 4))
    assert np.isclose(est, D.mean()) and df == 3 and np.isclose(se, D.std(ddof=1) / 2)


def test_stats_df_rule_and_decision():
    D = np.full(10, 0.001) + np.linspace(-0.002, 0.002, 10)
    r = stats.family_test(D, np.repeat(["h", "m"], 5), margin=0.03)
    assert r["df"] == 8 and r["noninferior"]
    assert list(stats.holm([0.001, 0.04, 0.03, 0.5, 0.2])) == [True, False, False, False, False]
    w = stats.family_test(
        -np.arange(1, 7) / 100, margin=0.03
    )  # floor ahead on every unit: signed ranks back non-inferiority, not superiority
    assert w["wilcoxon_noninferior_p"] == pytest.approx(1 / 64) and w["wilcoxon_superior_p"] == 1.0


def test_aggregate_units_and_margin():
    u, D, S = stats.aggregate_units(
        {"gse303070_120131_280N": 0.1, "gse303070_120131_280T": 0.3, "gse1_110": 0.0},
        {"gse303070_120131_280N": "Donor280", "gse303070_120131_280T": "Donor280", "gse1_110": "GSE1_110"},
    )
    assert dict(zip(u, D))["Donor280"] == pytest.approx(0.2) and set(S) == {"CosMx 1k", "CosMx WTx"}
    assert (
        stats.margin_from_dev("F3", [0.03, 0.05], [0, 0]) == pytest.approx(0.004)
        and stats.margin_from_dev("F1", [1.0], [0.0]) == 0.03
    )
    assert stats.margin_from_dev("F4", [0.8], [0]) == 0.03


def test_f2_scored_regions_drop_junk_and_absent():
    m, lab = f2.scored_regions(["A", "n/a", "B", "unassigned", "C", "A"], ["A", "B", "n/a"])
    assert m.tolist() == [True, False, True, False, False, True] and sorted(set(lab)) == ["A", "B"]
