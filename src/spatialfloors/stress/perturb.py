"""Input perturbations for the floor characterisation landscape (PLAN §14): count thinning, cell subsample, FOV hold-out,
gene-panel subsample. Each returns the perturbed object plus the kept indices, deterministic in seed."""
import numpy as np
import scipy.sparse as sp

NAMES = ("none", "thin50", "thin25", "sub50", "fov75", "genes75")


def thin(X, frac, seed):
    """Binomial read thinning of a count matrix (each read kept with probability frac)."""
    X = sp.csr_matrix(X, dtype=np.float32, copy=True); rng = np.random.default_rng(seed)
    X.data = rng.binomial(np.rint(X.data).astype(np.int64), frac).astype(np.float32); X.eliminate_zeros(); return X


def subsample(n, frac, seed):
    return np.sort(np.random.default_rng(seed).choice(n, int(round(frac * n)), replace=False))


def drop_fovs(fov, frac_keep, seed):
    """Cells of a random frac_keep of the FOVs (spatially contiguous hold-out)."""
    f = np.unique(np.asarray(fov).astype(str)); rng = np.random.default_rng(seed)
    keep = rng.choice(f, max(1, int(round(frac_keep * len(f)))), replace=False); return np.flatnonzero(np.isin(np.asarray(fov).astype(str), keep))


def apply(name, X, fov, seed):
    """-> (X', cell idx, gene idx) for a perturbation name in NAMES."""
    n, g = X.shape; allc, allg = np.arange(n), np.arange(g)
    if name == "none": return X, allc, allg
    if name.startswith("thin"): return thin(X, int(name[4:]) / 100, seed), allc, allg
    if name.startswith("sub"): c = subsample(n, int(name[3:]) / 100, seed); return X[c], c, allg
    if name.startswith("fov"): c = drop_fovs(fov, int(name[3:]) / 100, seed); return X[c], c, allg
    if name.startswith("genes"): gg = subsample(g, int(name[5:]) / 100, seed); return X[:, gg], allc, gg
    raise ValueError(name)
