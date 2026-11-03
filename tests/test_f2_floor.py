"""F2 floor on a planted toy (spatialfloors.sim.domain_toy): domains differ only in neighbourhood composition. The
other-lineage (bleed-safe) context fails it (ARI ~0: the context then depends on the receiver type); the all-type mean
recovers it."""
import numpy as np
from spatialfloors.sim import domain_toy
from spatialfloors.floors.simple import expr_pca, domain_features
from spatialfloors.metrics.f2 import ari_at_true_k, kmeans_labels, gmm_labels


def test_context_finds_planted_domains_expression_does_not():
    C, xy, fov, dom = domain_toy(); dom = dom.astype(str)
    F, _ = domain_features(C, xy, fov, 30, 50.); a_ctx = ari_at_true_k(dom, cluster_fn=lambda k, s: kmeans_labels(F, k, s), seeds=range(3))[0]
    a_exp = ari_at_true_k(dom, E=expr_pca(C, 30), seeds=range(3))[0]
    assert a_ctx > .9 and a_ctx > a_exp + .5, (a_ctx, a_exp)  # expression alone sees the type shares (~.25)


def test_bleedsafe_context_blind_to_composition_domains():
    """Pins the P2 finding: the other-lineage context cannot separate composition-only domains."""
    C, xy, fov, dom = domain_toy(); dom = dom.astype(str); F, _ = domain_features(C, xy, fov, 30, 50., ctx="bleedsafe", scale="dim")
    assert ari_at_true_k(dom, cluster_fn=lambda k, s: kmeans_labels(F, k, s), seeds=range(3))[0] < .2


def test_gmm_recovers_unequal_elongated_regions_kmeans_splits():
    """Pins the gse1_232 floor failure (PLAN §12): one large region elongated along an axis longer than the gap between
    regions; k-means at the true k cuts the large region in half, a full-covariance GMM recovers both."""
    rng = np.random.default_rng(0); n1, n2 = 8000, 2000
    A = rng.normal(size=(n1, 5)) * [6, 1, 1, 1, 1]; B = rng.normal(size=(n2, 5)) * [1, 1, 1, 1, 1] + [0, 4, 0, 0, 0]
    F = np.vstack([A, B]); dom = np.r_[np.zeros(n1), np.ones(n2)].astype(int).astype(str)
    a_km = ari_at_true_k(dom, cluster_fn=lambda k, s: kmeans_labels(F, k, s), seeds=range(3))[0]
    a_gm = ari_at_true_k(dom, cluster_fn=lambda k, s: gmm_labels(F, k, s), seeds=range(3))[0]
    assert a_gm > .8 and a_km < .5, (a_km, a_gm)
