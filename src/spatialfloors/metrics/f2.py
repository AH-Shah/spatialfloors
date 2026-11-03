"""F2 domains metric: ARI vs annotated regions at the true k, clustering seeds averaged (SPEC §F2)."""
import numpy as np
from sklearn.metrics import adjusted_rand_score

SEEDS = range(10)


def kmeans_labels(E, k, seed):
    from sklearn.cluster import KMeans
    return KMeans(k, n_init=10, random_state=seed).fit_predict(E)


def gmm_labels(E, k, seed):
    """Full-covariance GMM (floor clusterer since 2026-10-29, PLAN §12; CellCharter clusters with a GMM too). k-means
    assumes equal spherical clusters and split gse1_232's larger, elongated region along its main axis (ARI .27 vs
    GMM .93 on the same features; supervised ceiling .98)."""
    from sklearn.mixture import GaussianMixture
    return GaussianMixture(k, covariance_type="full", reg_covar=1e-4, random_state=seed).fit(E).predict(E)


def merge_labels(lab, k, how, xy=None, prob=None, knn=15):
    """Greedy merge of an over-clustered partition down to k clusters (label-free). how:
      "spatial"  merge the pair whose spatial kNN graph mixes most, edges between i and j / min(n_i, n_j): sub-clusters
                 of one heterogeneous region interleave in space, two distinct neighbouring regions only share a border
      "entropy"  merge the pair whose posteriors overlap most (largest entropy drop; Baudry et al. 2010), needs prob
      "ward"     smallest Ward cost n_i n_j / (n_i + n_j) |mu_i - mu_j|^2 on the features in prob (features x, not posteriors)"""
    lab = np.unique(np.asarray(lab), return_inverse=True)[1]; K = lab.max() + 1; grp = {i: [i] for i in range(K)}
    if how == "spatial":
        from scipy.spatial import cKDTree
        nb = cKDTree(xy).query(xy, knn + 1)[1][:, 1:]; A = np.zeros((K, K))
        np.add.at(A, (np.repeat(lab, knn), lab[nb].ravel()), 1); A = A + A.T
    while len(grp) > k:
        keys = list(grp); best, bv = None, -np.inf
        for a in range(len(keys)):
            for b in range(a + 1, len(keys)):
                i, j = keys[a], keys[b]; mi, mj = np.isin(lab, grp[i]), np.isin(lab, grp[j])
                if how == "spatial": v = A[np.ix_(grp[i], grp[j])].sum() / min(mi.sum(), mj.sum())
                elif how == "entropy":
                    ti, tj = prob[:, grp[i]].sum(1), prob[:, grp[j]].sum(1); xl = lambda t: t * np.log(np.maximum(t, 1e-300))
                    v = float((xl(ti + tj) - xl(ti) - xl(tj)).sum())
                else: ni, nj = mi.sum(), mj.sum(); v = -ni * nj / (ni + nj) * float(((prob[mi].mean(0) - prob[mj].mean(0)) ** 2).sum())
                v = -np.inf if v != v else v
                if best is None or v > bv: best, bv = (i, j), v
        i, j = best; grp[i] += grp.pop(j)
    out = np.empty_like(lab)
    for n, g in enumerate(grp.values()): out[np.isin(lab, g)] = n
    return out


def gmm_merge_labels(E, k, seed, xy, kprime, how="spatial"):
    """GMM at kprime > k components, merged to k (merge_labels)."""
    from sklearn.mixture import GaussianMixture
    g = GaussianMixture(kprime, covariance_type="full", reg_covar=1e-4, random_state=seed).fit(E); lab = g.predict(E)
    return merge_labels(lab, k, how, xy=xy, prob=g.predict_proba(E) if how == "entropy" else E)


JUNK = {"", "nan", "None", "n/a", "NA", "unassigned", "brain-unassigned"}


def scored_regions(regions, present):
    """(mask, labels) of the cells scored by F2: region in the crop's present list (>= 2% of the slide's labelled cells,
    data.crops) minus non-region labels (JUNK). Every cell is clustered; k = number of present regions (review
    2026-10-28: Allen 'n/a' / 'unassigned' and sub-2% regions were counted as regions)."""
    keep = [r for r in present if str(r) not in JUNK]; reg = np.asarray(regions, dtype=object).astype(str)
    m = np.isin(reg, keep); return m, reg[m]


def ari_bal(y, s):
    """Balanced ARI (F2 primary): ARI on equal-size resamples, 200 cells of each region with >= 200 cells, 10 draws
    (every region weighs the same, so a partition cannot score by getting the big regions and dropping a small one)."""
    rng = np.random.default_rng(0); big = [c for c in np.unique(y) if (y == c).sum() >= 200]; v = []
    for _ in range(10):
        ii = np.concatenate([rng.choice(np.flatnonzero(y == c), 200, replace=False) for c in big]); v.append(adjusted_rand_score(y[ii], s[ii]))
    return float(np.mean(v))


def ari_bal_saved(path, key, mask, regions):
    """Seed-mean balanced ARI of a saved partition (labels npz: key -> seeds x cells, or one row) on the scored cells."""
    z = np.load(path, allow_pickle=True); L = np.asarray(z[key]).astype(str); L = L[None] if L.ndim == 1 else L
    return float(np.mean([ari_bal(regions, l[mask]) for l in L]))


def ari_at_true_k(regions, cluster_fn=None, E=None, seeds=SEEDS):
    """mean ARI over seeds of cluster_fn(k, seed) vs regions, k = number of regions present. Default: k-means on E."""
    regions = np.asarray(regions); k = len(np.unique(regions))
    fn = cluster_fn or (lambda k_, s: kmeans_labels(E, k_, s))
    aris = [adjusted_rand_score(regions, fn(k, s)) for s in seeds]
    return float(np.mean(aris)), aris
