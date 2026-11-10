"""F3 floor: bleed-safe ridge, niche -> receiver state, quadrant CV (metrics.f3.block_cv).

Shared by every F3 model (prepare): receivers = Leiden ids on Expr50 of the panel genes, lineages, local foreign mask
(floors.bleedsafe; label-free, computed once on the crop); quadrants of the whole crop; within-quadrant neighbourhood.
Target of receiver r: its masked genes renormalised on them (foreign bleed cannot enter via CP10k), FOV-centred; depth
(log detected masked genes) is a covariate removed with train-only coefficients.
Floor features: kernel id mix T and bleed-safe context (SVD of the genes r does not dominate, fitted on training rows,
FOV-centred, kernel mean). Kernel: Gaussian ell um over the 200 NN, same FOV, same quadrant, other lineages.
Fold-dependent fits (context SVD and gene scale, target covariates / z-score / PCA, feature scaler, ridge alpha) see
training quadrants only (asserted by FitLog in block_cv)."""

import numpy as np
import scipy.sparse as sp
from . import bleedsafe as bs
from .simple import lognorm
from ..metrics import f3 as m3

# receptive field shared by every F3 model: GAT 2 layers x 15-NN radius (21 um on gse1_232) ~ floor 2 x ell 25 um
RF_UM = 50.0
MIN_CELLS, KNN, ALPHAS = 1000, 200, np.logspace(-2, 4, 13)


def prepare(X, xy_um, fov, ids=None, seed=0):
    """X counts (cells x panel genes). Returns the shared F3 design dict."""
    from scipy.spatial import cKDTree

    x = lognorm(X)
    fov = np.asarray(fov).astype(str)
    xy_um = np.asarray(xy_um, float)
    from . import lean

    if lean.on():
        fov = np.unique(fov, return_inverse=True)[1].astype(np.int32)  # same partition as int codes (whole-slide arm)
    if ids is None:
        ids = bs.leiden_ids(bs.expr50(x, seed), seed=seed)
    ids = np.unique(np.asarray(ids), return_inverse=True)[1]
    K = int(ids.max()) + 1
    mu_c = bs.id_means(x, ids, K)
    lin_c = bs.lineages(mu_c)
    nb15 = cKDTree(xy_um).query(xy_um, 16)[1][:, 1:]
    mask, _ = bs.local_foreign_mask(mu_c, ids, lin_c, nb15, fov)
    blocks, split = m3.quadrant_blocks(xy_um)
    return dict(x=x, xy=xy_um, fov=fov, ids=ids, K=K, lin=lin_c[ids], mask=mask, blocks=blocks, split=split)


def kernel(D, ell_um, knn=KNN, same_lineage=False):
    """200-NN Gaussian kernel (ell in um) restricted to the same FOV and quadrant; other lineages unless
    same_lineage."""
    from scipy.spatial import cKDTree
    from . import lean

    if lean.on():
        return lean.knn_gauss(
            D["xy"],
            ell_um,
            knn,
            [(D["fov"], True), (D["blocks"], True)] + ([] if same_lineage else [(D["lin"], False)]),
        )  # whole-slide arm
    d, nb = cKDTree(D["xy"]).query(D["xy"], knn + 1)
    d, nb = d[:, 1:], nb[:, 1:]
    keep = (D["fov"][nb] == D["fov"][:, None]) & (D["blocks"][nb] == D["blocks"][:, None])
    if not same_lineage:
        keep &= D["lin"][nb] != D["lin"][:, None]
    return nb, (np.exp(-(d**2) / (2 * ell_um**2)) * keep).astype(np.float32)


def receivers(D, min_cells=MIN_CELLS):
    return [r for r in range(D["K"]) if (D["ids"] == r).sum() >= min_cells and D["mask"][r].sum() > m3.NPC]


def target(D, r):
    """(rows, Y, depth) of receiver r: masked genes renormalised on them, FOV-centred; depth = FOV-centred log
    detected."""
    rows = np.flatnonzero(D["ids"] == r)
    g = np.flatnonzero(D["mask"][r])
    xr = bs.renorm_log(D["x"][rows][:, g])
    Y = bs.group_center(xr.toarray(), D["fov"][rows])
    dep = bs.group_center(np.log1p(np.diff(xr.indptr)).astype(float)[:, None], D["fov"][rows])
    return rows, Y.astype(np.float32), dep


def id_mix(D, nb, w):
    return bs.id_mix(nb, w, D["ids"], D["K"])


def context_fn(D, r, rows, nb, w, n_pc=50, fitlog=None, seed=0):
    """Callable for block_cv: fold q -> bleed-safe context of `rows` with the SVD and gene scale fitted on cells outside
    quadrant q (fitlog records those rows, in receiver-row coordinates: the receiver's own training rows)."""
    from sklearn.decomposition import TruncatedSVD

    g = np.flatnonzero(~D["mask"][r])
    xg = bs.renorm_log(D["x"][:, g])
    k = min(n_pc, len(g) - 1)

    def f(q, tr_local):
        if fitlog is not None:
            fitlog(tr_local)
        fit = np.flatnonzero(D["blocks"] != q)
        fit = np.sort(np.random.default_rng(seed).choice(fit, min(len(fit), 100000), replace=False))
        m1 = np.asarray(xg[fit].mean(0)).ravel()
        sd = np.sqrt(np.maximum(np.asarray(xg[fit].multiply(xg[fit]).mean(0)).ravel() - m1**2, 1e-6))
        xs = sp.csr_matrix(xg @ sp.diags(1 / sd))
        E = bs.group_center(TruncatedSVD(k, random_state=seed).fit(xs[fit]).transform(xs), D["fov"]).astype(np.float32)
        return bs.kernel_mean(nb, w, E, rows)

    return f


def groups(fov, xy_um, tile_um=500.0):
    """Spatial groups for inner validation: FOVs when a unit has >= 5, else 500 um tiles (MERFISH
    sections are one FOV). Random-row inner splits leak spatial autocorrelation."""
    fov = np.asarray(fov).astype(str)
    if len(np.unique(fov)) >= 5:
        return fov
    t = np.floor(np.asarray(xy_um) / tile_um).astype(int)
    return np.char.add(np.char.add(t[:, 0].astype(str), "_"), t[:, 1].astype(str))


def ridge_fit_predict(Xtr, Ytr, Xte):
    from sklearn.linear_model import RidgeCV

    return RidgeCV(alphas=ALPHAS).fit(Xtr, Ytr).predict(Xte)


def score_floor(
    D,
    ell_um=25.0,
    n_pc=50,
    min_cells=MIN_CELLS,
    seed=0,
    features="T+ctx",
    fit_predict=ridge_fit_predict,
    full=False,
    shuffle_seed=None,
):
    """Per receiver: block_cv score of ridge on [T, ctx_bs]; features "T", "ctx" or "T+ctx". Returns {receiver: (score,
    folds, per_pc)}."""
    nb, w = kernel(D, ell_um)
    T = id_mix(D, nb, w)
    out = {}
    for r in receivers(D, min_cells):
        rows, Y, dep = target(D, r)
        ctx = context_fn(D, r, rows, nb, w, n_pc, seed=seed)
        Xf = {"T": lambda q, tr: T[rows], "ctx": ctx, "T+ctx": lambda q, tr: np.c_[T[rows], ctx(q, tr)]}[features]
        if shuffle_seed is not None:
            Xf = _shuffled(
                Xf, D["fov"][rows], D["blocks"][rows], shuffle_seed
            )  # null: features permuted within FOV x quadrant
        det = {"full": full}
        s, fo = m3.block_cv(
            Xf, Y, D["xy"][rows], fit_predict, RF_UM, blocks=D["blocks"][rows], split=D["split"], covars=dep, detail=det
        )
        out[r] = (s, fo, np.nanmean(det.get("per_pc", [[np.nan]]), 0).tolist()) + (
            (det,) if full else ()
        )  # per-PC Spearman averaged over folds (diagnostic)
    return out


def _shuffled(Xf, fov, blocks, seed):
    """Features with rows permuted within FOV x quadrant: keeps the feature distribution and the train / test split,
    breaks the cell-to-niche pairing."""
    g = np.char.add(np.asarray(fov).astype(str), np.asarray(blocks).astype(str))
    perm = np.arange(len(g))
    rng = np.random.default_rng(seed)
    for v in np.unique(g):
        i = np.flatnonzero(g == v)
        perm[i] = rng.permutation(i)
    return lambda q, tr: Xf(q, tr)[perm]


def unit_score(per_receiver):
    """Unit F3 score = mean over receivers (finite scores only)."""
    v = [r[0] for r in per_receiver.values() if np.isfinite(r[0])]
    return float(np.mean(v)) if v else float("nan")


# ---------- v3: same-lineage state rings, sender-lineage-resolved context, one ridge penalty per group ----------
RINGS_V3 = {
    "self": (50.0, 150.0, 600),
    "far": (150.0, 400.0, 1500),
    # (lo um, hi um, nearest cells of any lineage queried); far ring reaches a median 220-330 um on CosMx, 400 um on MERFISH
}
V3_ALPHAS, V3_MULTS, V3_NPC_OWN, V3_NPC_TYPED = (
    (300.0, 1000.0, 3000.0, 10000.0),
    (0.3, 1.0, 3.0, None),
    10,
    10,
)  # None = group dropped


def ring(D, lo, hi, kmax):
    """Same-lineage cells at lo <= d < hi um among the kmax nearest, same FOV and quadrant -> (nb, w) as kernel.
    Bleed-safe: no contact spill reaches beyond 50 um."""
    from . import lean

    if lean.on():
        return lean.ring(D, lo, hi, kmax)  # whole-slide arm only (SF_LEAN=1): sparse, chunked query
    from scipy.spatial import cKDTree

    n = len(D["xy"])
    dd, ii = cKDTree(D["xy"]).query(D["xy"], kmax, distance_upper_bound=hi)
    ok = (ii < n) & (dd >= lo)
    ii = np.where(ok, ii, 0)
    ok &= (
        (D["lin"][ii] == D["lin"][:, None])
        & (D["fov"][ii] == D["fov"][:, None])
        & (D["blocks"][ii] == D["blocks"][:, None])
    )
    return ii.astype(np.int32), ok.astype(np.float32)


def svd_fn(D, gsel, n_pc, seed=0):
    """fold q -> FOV-centred SVD scores of every cell on genes gsel; SVD and gene scale fitted on (<= 100k) cells
    outside
    quadrant q, as context_fn; cached per fold."""
    from sklearn.decomposition import TruncatedSVD

    xg = bs.renorm_log(D["x"][:, gsel])
    k = min(n_pc, len(gsel) - 1)
    cache = {}

    def f(q):
        if q not in cache:
            fit = np.flatnonzero(D["blocks"] != q)
            fit = np.sort(np.random.default_rng(seed).choice(fit, min(len(fit), 100000), replace=False))
            m1 = np.asarray(xg[fit].mean(0)).ravel()
            sd = np.sqrt(np.maximum(np.asarray(xg[fit].multiply(xg[fit]).mean(0)).ravel() - m1**2, 1e-6))
            xs = sp.csr_matrix(xg @ sp.diags(1 / sd))
            cache[q] = bs.group_center(TruncatedSVD(k, random_state=seed).fit(xs[fit]).transform(xs), D["fov"]).astype(
                np.float32
            )
        return cache[q]

    return f


def typed_mean(D, nb, w, E, rows, tau=None):
    """Kernel mean of E resolved by sender lineage (one block per lineage), shrunk toward 0 by tau (bs.TAU_W)."""
    tau = bs.TAU_W if tau is None else tau
    cols = []
    for l in np.unique(D["lin"]):
        wl = w[rows] * (D["lin"][nb[rows]] == l)
        cols.append(np.einsum("nk,nkd->nd", wl, E[nb[rows]]) / (wl.sum(1, keepdims=True) + tau))
    return np.hstack(cols).astype(np.float32)


class GroupFeatures:
    """Callable features for block_cv built from named groups [(name, fn(q, tr) -> array)]; remembers the group of every
    column block_cv keeps (it drops training-constant columns with the same rule)."""

    def __init__(self, parts):
        self.parts, self.names = parts, [p[0] for p in parts]

    def __call__(self, q, tr):
        blocks = [f(q, tr) for _, f in self.parts]
        cg = np.concatenate([np.full(b.shape[1], i) for i, b in enumerate(blocks)])
        X = np.hstack(blocks)
        self.live_group = cg[X[tr].std(0) > 1e-6]
        return X


class GroupRidge:
    """block_cv fit_predict (wants_rows): ridge with penalty alpha * m_g per column group; group 0
    (ctx) keeps m 1, every other group m in V3_MULTS (None = dropped); (alpha, m) chosen per fold by inner CV on the
    training rows grouped by sgroups (the receiver's quadrant ids: leave one training quadrant out, as the outer split),
    criterion = mean Pearson r over the target PCs (the metric is a rank correlation per PC). choices: per fold (alpha,
    m)."""

    wants_rows = True

    def __init__(self, feats, sgroups, alphas=V3_ALPHAS, mults=V3_MULTS):
        self.feats, self.sg, self.alphas, self.mults, self.choices = feats, np.asarray(sgroups), alphas, mults, []

    def _solve(self, X, Y, pen):
        xm, ym = X.mean(0), Y.mean(0)
        Xc = X - xm
        B = np.linalg.solve(Xc.T @ Xc + np.diag(pen), Xc.T @ (Y - ym))
        return B, xm, ym

    @staticmethod
    def _corr(P, Y):
        P, Y = P - P.mean(0), Y - Y.mean(0)
        return float(np.nanmean((P * Y).sum(0) / np.sqrt((P**2).sum(0) * (Y**2).sum(0) + 1e-12)))

    def __call__(self, Xtr, Ytr, Xte, tr, te, q):
        import itertools

        cg = self.feats.live_group
        ng = len(self.feats.names)
        g = self.sg[tr]
        combos = [(a, (1.0,) + m) for a in self.alphas for m in itertools.product(self.mults, repeat=ng - 1)]
        pens = [np.array([a * (1e12 if mm[c] is None else mm[c]) for c in cg]) for a, mm in combos]
        sc = np.zeros(len(combos))
        gu = np.unique(g)
        for v in gu if len(gu) >= 2 else []:
            itr, iva = np.flatnonzero(g != v), np.flatnonzero(g == v)
            for i, pen in enumerate(pens):
                B, xm, ym = self._solve(Xtr[itr], Ytr[itr], pen)
                sc[i] += self._corr((Xtr[iva] - xm) @ B + ym, Ytr[iva])
        best = (
            int(np.argmax(sc))
            if len(gu) >= 2
            else next(
                i for i, (a, mm) in enumerate(combos) if a == 1000.0 and mm[1] == 1.0 and all(x is None for x in mm[2:])
            )
        )
        self.choices.append(dict(q=int(q), alpha=combos[best][0], mult=dict(zip(self.feats.names, combos[best][1]))))
        B, xm, ym = self._solve(Xtr, Ytr, pens[best])
        return (Xte - xm) @ B + ym


def lock_part(rows, lock, K_lock, T_lock, Ec):
    """The v2 floor's features ('T+ctx' or 'ctx' at the lock kernel) as a block_cv callable."""
    nb, w = K_lock
    return (
        (lambda q, tr: bs.kernel_mean(nb, w, Ec(q), rows))
        if lock == "ctx"
        else (lambda q, tr: np.c_[T_lock[rows], bs.kernel_mean(nb, w, Ec(q), rows)])
    )


def v3_parts(D, r, rows, lock, K_lock, T_lock, rings, Ec, Es, shuffle=None, shuffle_typed=False):
    """Feature groups of receiver r: ctx (bleed-safe context at the lock kernel; kept, m 1), T (id mix at the lock
    kernel;
    both species, so the type mix can be weighted or dropped), self / far rings (same lineage, the receiver's own masked
    genes), typed (sender-lineage context, first V3_NPC_TYPED context components). shuffle: a permutation of the
    receiver rows applied to the ring groups (null), and to typed too with shuffle_typed (null of every v3-only group).
    `lock` names the v2 features (lock_part); v3 searches over them."""
    nb, w = K_lock
    perm = np.arange(len(rows)) if shuffle is None else shuffle
    parts = [("ctx", lambda q, tr: bs.kernel_mean(nb, w, Ec(q), rows)), ("T", lambda q, tr: T_lock[rows])]
    for name, (rn, rw) in rings.items():
        parts.append((name, (lambda rn, rw: lambda q, tr: bs.kernel_mean(rn, rw, Es(q), rows)[perm])(rn, rw)))
    pt = perm if shuffle_typed else np.arange(len(rows))
    parts.append(("typed", lambda q, tr: typed_mean(D, nb, w, Ec(q)[:, :V3_NPC_TYPED], rows)[pt]))
    return GroupFeatures(parts)


V3_LOCK = {
    "human": dict(feat="T+ctx", ell=25.0, npc=10),
    "mouse": dict(feat="ctx", ell=100.0, npc=50),
}  # v2 lock kernels


def v3_floor_parts(rows, lock, K_lock, T_lock, rings, Es, Ec):
    """F3 floor v3 as locked: the v2 lock features + the same-lineage rings
    50-150 um (600 nearest) and 150-400 um (1,500 nearest; median reach 220-330 um on CosMx) on the receiver's own
    masked genes, one ridge (alpha 1000), no penalty search."""
    lk = lock_part(rows, lock, K_lock, T_lock, Ec)
    return lambda q, tr: np.hstack([lk(q, tr)] + [bs.kernel_mean(*rings[k], Es(q), rows) for k in ("self", "far")])


def score_v3(D, species, fit_predict=None, min_cells=MIN_CELLS, alpha=1000.0):
    """Per receiver block_cv score of the F3 floor v3 (fit_predict None: v3_floor_parts + ridge, locked alpha 1000), or
    of any
    block_cv fit_predict on every v3 feature group (v3_parts: ctx, T, rings, sender-typed context; a superset of the
    floor's features, for LightGBM, the F3 test pair). Returns {receiver: (score, folds, [nan])}."""
    from sklearn.linear_model import Ridge

    L = V3_LOCK[species]
    K = kernel(D, L["ell"])
    T = id_mix(D, *K)
    R = {k: ring(D, *v) for k, v in RINGS_V3.items()}
    out = {}
    for r in receivers(D, min_cells):
        rows, Y, dep = target(D, r)
        Ec = svd_fn(D, np.flatnonzero(~D["mask"][r]), L["npc"])
        Es = svd_fn(D, np.flatnonzero(D["mask"][r]), V3_NPC_OWN)
        if fit_predict is None:
            F, fp = (
                v3_floor_parts(rows, L["feat"], K, T, R, Es, Ec),
                (lambda Xtr, Ytr, Xte: Ridge(alpha=alpha).fit(Xtr, Ytr).predict(Xte)),
            )
        else:
            F, fp = v3_parts(D, r, rows, L["feat"], K, T, R, Ec, Es), fit_predict
        s, fo = m3.block_cv(F, Y, D["xy"][rows], fp, RF_UM, blocks=D["blocks"][rows], split=D["split"], covars=dep)
        out[r] = (s, fo, [np.nan])
    return out
