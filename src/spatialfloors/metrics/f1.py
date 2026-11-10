"""F1 identity metric: count-split kNN Poisson-deviance score (PROTOCOL §2)."""

import numpy as np
import scipy.sparse as sp

K = 15
LAM = 0.5  # prior weight on pbar (P2 dev ; was 1/(k+1), negative scores at ~32 counts/cell)


def count_split(X, seed=0):
    """X (cells x genes, integer counts, csr) -> X_a ~ Binomial(X, .5), X_b = X - X_a (same sparsity pattern)."""
    X = sp.csr_matrix(X)
    rng = np.random.default_rng(seed)
    a = X.copy()
    a.data = rng.binomial(X.data.astype(np.int64), 0.5).astype(X.dtype)
    b = X - a
    a.eliminate_zeros()
    b.eliminate_zeros()
    return a, b


def knn_cosine(E, k=K):
    """k nearest neighbours by cosine on rows of E, self excluded (explicitly, also under exact duplicates)."""
    from sklearn.neighbors import NearestNeighbors

    U = E / np.maximum(np.linalg.norm(E, axis=1, keepdims=True), 1e-12)
    nb = NearestNeighbors(n_neighbors=k + 1).fit(U).kneighbors(U, return_distance=False)
    self_ = nb == np.arange(len(U))[:, None]
    drop = np.where(self_.any(1), np.argmax(self_, 1), k)  # drop self if present, else the (k+1)-th neighbour
    keep = np.ones_like(nb, bool)
    keep[np.arange(len(U)), drop] = False
    return nb[keep].reshape(len(U), k)


def poisson_dev(x, mu):
    """Elementwise Poisson deviance 2[x log(x/mu) - (x - mu)], 0 log 0 = 0; x, mu dense float64."""
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(x > 0, x * np.log(x / mu), 0.0)
    return 2 * (t - (x - mu))


def knn_poisson_score(E, Xa, Xb, genes, k=K, chunk=20000, lam=LAM, per_gene=False):
    """mu_hat = L_b ((1 - lam) neighbour mean of X_a / L_a + lam pbar_a). 1 - sum Dev(X_b, mu_hat) / sum Dev(X_b, L_b
    pbar_a) over `genes` (HVG column indices chosen on X_a) and all cells.
    L = total counts over all columns of X_a / X_b; neighbours of each cell from E (self excluded).
    per_gene: also return the per-gene deviance sums (num, den) and the gene indices."""
    Xa, Xb = sp.csr_matrix(Xa, dtype=np.float64), sp.csr_matrix(Xb, dtype=np.float64)
    La, Lb = np.asarray(Xa.sum(1)).ravel(), np.asarray(Xb.sum(1)).ravel()
    genes = np.asarray(genes)
    genes = genes[np.asarray(Xa[:, genes].sum(0)).ravel() > 0]  # pbar_a > 0 (HVGs on X_a always are)
    A = Xa[:, genes]
    Bg = Xb[:, genes]
    pbar = np.asarray(A.sum(0)).ravel() / La.sum()
    F = sp.csr_matrix(
        sp.diags(np.where(La > 0, 1 / np.maximum(La, 1e-300), 0.0)) @ A
    )  # X_a / L_a, zero for empty cells
    nb = knn_cosine(np.asarray(E, np.float64), k)
    n = A.shape[0]
    num = den = np.zeros(len(genes))
    for s in range(0, n, chunk):
        r = np.arange(s, min(n, s + chunk))
        P = sp.csr_matrix(
            (np.full(len(r) * k, 1.0 / k), (np.repeat(np.arange(len(r)), k), nb[r].ravel())), shape=(len(r), n)
        )
        prof = (P @ F).toarray() * (1 - lam) + pbar * lam
        x = Bg[r].toarray()
        mu = Lb[r, None] * prof
        mu0 = Lb[r, None] * pbar
        z = Lb[r] == 0
        x, mu, mu0 = x[~z], mu[~z], mu0[~z]  # L_b = 0: x = mu = 0, deviance 0
        num = num + poisson_dev(x, mu).sum(0)
        den = den + poisson_dev(x, mu0).sum(0)
    return (1 - num.sum() / den.sum(), num, den, genes) if per_gene else 1 - num.sum() / den.sum()


def chance_level(E, Xa, Xb, genes, k=K, n_perm=3, seed=0, lam=LAM):
    """Score of E with rows permuted (neighbours unrelated to identity), mean over n_perm permutations. Reported next to
    every F1 score; it is < 0 because a k-cell neighbour mean is a noisier rate estimate than pbar."""
    rng = np.random.default_rng(seed)
    return float(
        np.mean(
            [
                knn_poisson_score(np.asarray(E)[rng.permutation(len(E))], Xa, Xb, genes, k, lam=lam)
                for _ in range(n_perm)
            ]
        )
    )
