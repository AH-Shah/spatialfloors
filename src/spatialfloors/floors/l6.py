"""F4 floor L6: pseudobulk two-way model on Cellina's leave-one-type-out benchmark (PROTOCOL §2), with the ladder of
simpler baselines it was chosen from. Data locations are arguments (`RAW`, `out`, `merfish_path`)."""

import time, resource
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from . import _labels_to_coarse as lm

from spatialfloors.paths import RAW as _RAW

RAW = _RAW / "GSE1"  # GSE1 CosMx WTx h5ads (cosmx_wt_{sid}.h5ad); override per call
K_NB, NPC, KMATCH, KAPPA, NDEG, PSC = 200, 50, 30, 10, 50, 0.1
CAP_H, CAP_BASE, CAP_POOL = 50000, 20000, 200000
ALPHA_REL = [1e-4, 1e-3, 1e-2, 1e-1, 1.0]
M_SECTIONS = ["C57BL6J-2.036", "C57BL6J-2.039", "C57BL6J-2.041"]
M_TYPES = ["GABAergic neuron", "astrocyte", "endothelial cell", "glutamatergic neuron", "oligodendrocyte"]
C_TYPES = ["Endothelial", "Epithelial", "Fibroblast", "Myeloid", "T_cell"]
RUNGS = [
    "baseline_repro",
    "L1_shared_shift",
    "L2_contamination",
    "L3_h_niche_ridge",
    "L4_niche_match",
    "L5_pooled_ridge",
    "L3null",
]


def log(*a):
    print(
        f"[{time.strftime('%H:%M:%S')} {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30:.1f}G]",
        *a,
        flush=True,
    )


# ---------- data (Cellina's preprocessing: NaN drop, filter 3/3, seurat_v3 HVG on counts) ----------
def prep(a, ntop):
    sc.pp.filter_cells(a, min_counts=3)
    sc.pp.filter_genes(a, min_counts=3)
    sc.pp.highly_variable_genes(
        a, flavor="seurat_v3", n_top_genes=ntop, subset=True
    )  # X holds counts (== Cellina's layer='counts')
    return (
        sp.csr_matrix(a.X, dtype=np.float64),
        np.asarray(a.obsm["spatial"], float),
        a.obs.ct.astype(str).values,
        a.obs.dom.astype(str).values,
    )


def load_crc(sid, raw=None):
    a = ad.read_h5ad(Path(raw or RAW) / f"cosmx_wt_{sid}.h5ad")
    a.obs_names_make_unique()
    ct = a.obs["ist"].astype(str).map(lm.LABEL_TO_COARSE)
    dom = a.obs["typ"].astype(str).str.extract(r"(REF|TVA|CRC)", expand=False)
    a.obsm["spatial"] = a.obs[["CenterX_global_px", "CenterY_global_px"]].values.astype(float)
    a.obs = pd.DataFrame(dict(ct=ct.values, dom=dom.values), index=a.obs_names)
    a._inplace_subset_obs((ct.notna() & dom.notna()).values)
    log(f"crc_{sid} loaded", a.shape)
    return prep(a, 2000)


def load_merfish(merfish_path):
    A = ad.read_h5ad(merfish_path)
    for sid in M_SECTIONS:
        m = (A.obs.brain_section_label == sid).values
        obs = pd.DataFrame(
            dict(ct=A.obs.cell_type.values[m], dom=A.obs.major_brain_region.values[m]), index=A.obs_names[m]
        )
        a = ad.AnnData(X=A.raw.X[np.flatnonzero(m)], obs=obs, var=pd.DataFrame(index=A.raw.var_names))
        a.obs_names_make_unique()
        a.obsm["spatial"] = np.asarray(A.obsm["X_spatial_coords"][m], float)
        a = a[(obs.ct.notna() & obs.dom.notna()).values].copy()
        log(sid, "loaded", a.shape)
        yield sid, prep(a, 1120)


# ---------- helpers ----------
def nmean(C, rows=None, w=None):
    """Mean over rows of Cellina-normalised CP10k (x / (sum + 1e-8) * 1e4), float64; optional row weights."""
    M = C if rows is None else C[rows]
    lib = np.asarray(M.sum(1)).ravel()
    f = 1e4 / (lib + 1e-8)
    if w is not None:
        f = f * w / w.sum()
        return np.asarray(M.T @ f).ravel()
    return np.asarray(M.T @ f).ravel() / M.shape[0]


def lognorm_dense(C, rows):
    M = C[rows]
    lib = np.asarray(M.sum(1)).ravel()
    return np.log1p((sp.diags(1e4 / (lib + 1e-8)) @ M).toarray()).astype(np.float32)


def niche_rows(xy, rows, src_mask, bw):
    """Row-normalised Gaussian kernel (cutoff .1, <= K_NB nearest) from `rows` to non-self cells in src_mask; |rows| x n
    csr."""
    src = np.flatnonzero(src_mask)
    tree = cKDTree(xy[src])
    kq = min(K_NB + 1, len(src))
    R, Cc, V = [], [], []
    for s0 in range(0, len(rows), 50000):
        q = rows[s0 : s0 + 50000]
        d, j = tree.query(xy[q], k=kq, workers=-1)
        j = src[j]
        w = np.exp(-(d**2) / (2 * bw**2))
        me = j == q[:, None]
        w[me] = 0
        w[~me.any(1), -1] = 0
        w[w < 0.1] = 0
        r, c = np.nonzero(w)
        R.append((r + s0).astype(np.int32))
        Cc.append(j[r, c].astype(np.int32))
        V.append(w[r, c].astype(np.float32))
    W = sp.csr_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(Cc))), shape=(len(rows), len(xy)))
    dg = np.asarray(W.sum(1)).ravel()
    return sp.diags(np.where(dg > 0, 1 / np.maximum(dg, 1e-12), 0).astype(np.float32)) @ W


def ridge(X, Y, score, rng):
    """Ridge with alpha = n * ALPHA_REL picked on a random 10% (scored on rows in `score`); returns predictor of shifts
    (B) and the alpha."""
    n = len(X)
    val = rng.random(n) < 0.1
    sv = val & score

    def fit(m, a):
        xm = X[m].mean(0)
        Xc = X[m] - xm
        e, Q = np.linalg.eigh(Xc.T @ Xc)
        XtY = np.asarray((Y[m].T @ Xc).T) if sp.issparse(Y) else Xc.T @ Y[m]
        return lambda al: (Q * (1 / (e + al))) @ (Q.T @ XtY), xm, np.asarray(Y[m].mean(0)).ravel()

    solve, xm, ym = fit(~val, None)
    Yv = Y[sv].toarray() if sp.issparse(Y) else Y[sv]
    errs = [np.mean(((X[sv] - xm) @ solve(r * (~val).sum()) + ym - Yv) ** 2) for r in ALPHA_REL]
    a = ALPHA_REL[int(np.argmin(errs))]
    solve, _, _ = fit(np.ones(n, bool), None)
    return solve(a * n), a


def group_stats(C, ct, dom, lib):
    """Per (type, domain) group: n, CP10k mean, raw gene totals, raw library total."""
    keys = sorted(set(zip(ct, dom)))
    n, M, T, L = [], [], [], []
    for t, d in keys:
        r = np.flatnonzero((ct == t) & (dom == d))
        n.append(len(r))
        M.append(nmean(C, r))
        T.append(np.asarray(C[r].sum(0)).ravel())
        L.append(lib[r].sum())
    return keys, np.array(n), np.array(M), np.array(T), np.array(L)


def wls(X, y, w, lam=1e-6):
    """Per-gene weighted least squares; X k x p (shared) or k x G x p (per gene), y/w k x G -> G x p."""
    Xg = X if X.ndim == 3 else np.broadcast_to(X[:, None, :], (X.shape[0], y.shape[1], X.shape[1]))
    A = np.einsum("kgp,kg,kgq->gpq", Xg, w, Xg) + lam * np.eye(Xg.shape[2])
    b = np.einsum("kgp,kg->gp", Xg, w * y)
    return np.linalg.solve(A, b[..., None])[..., 0]


def twoway(gs, h, hds, nu=None, rank=0, iters=20, gw=None):
    """Per-gene WLS of log(group CP10k mean + half a count) on type + domain one-hots over groups with >= 5 cells, (h,
    held-out
    domains) excluded; predicts (h, d). L6: additive. L7 (nu): + gamma_g * log(group neighbourhood CP10k mean + eps).
    L8 (rank): + rank-r type x domain interaction u_t . v_dg (ALS; u shared over genes).
    gw: per-group multipliers of the WLS weights (similarity weighting)."""
    keys, n, M, T, L = gs
    t_ = np.array([k[0] for k in keys])
    d_ = np.array([k[1] for k in keys])
    keep = (n >= 5) & ~((t_ == h) & np.isin(d_, hds))
    ts, ds = sorted(set(t_[keep])), sorted(set(d_[keep]))
    tid = np.array([ts.index(t) if t in ts else -1 for t in t_])
    did = np.array([ds.index(d) if d in ds else -1 for d in d_])
    OT, OD = np.eye(len(ts))[tid], np.eye(len(ds))[did]
    X = np.hstack([OT, OD])
    eps = (0.5 * 1e4 / L)[:, None]
    y = np.log(M + eps)
    w = 1 / (1 / (T + 0.5) + 0.05)
    if gw is not None:
        w = w * np.asarray(gw, float)[:, None]
    tgt = {hd: next(i for i, (t, d) in enumerate(keys) if t == h and d == hd) for hd in hds}
    if nu is not None:
        Xg = np.concatenate(
            [np.broadcast_to(X[:, None, :], (len(X), M.shape[1], X.shape[1])), np.log(nu + eps)[..., None]], 2
        )
        B = wls(Xg[keep], y[keep], w[keep])
        return {hd: np.exp(np.einsum("gp,gp->g", Xg[i], B)) for hd, i in tgt.items()}
    if rank == 0:
        B = wls(X[keep], y[keep], w[keep])
        return {hd: np.exp(X[i] @ B.T) for hd, i in tgt.items()}
    rng = np.random.default_rng(0)
    U = rng.normal(0, 0.1, (len(ts), rank))
    nd = len(ds)
    for _ in range(iters):
        Xi = np.hstack(
            [X, (U[tid][:, :, None] * OD[:, None, :]).reshape(len(X), -1)]
        )  # columns: types, domains, u_r x domain
        B = wls(Xi[keep], y[keep], w[keep])
        V = B[:, len(ts) + nd :].reshape(-1, rank, nd)  # G x r x D
        R = y - X @ B[:, : len(ts) + nd].T  # residual after additive part, k x G
        for j in range(len(ts)):  # u_t: WLS over (d, g) of residual on v_dg
            m = keep & (tid == j)
            F = V[:, :, did[m]].transpose(2, 0, 1).reshape(-1, rank)
            ww = w[m].ravel()
            rr = R[m].ravel()
            U[j] = np.linalg.solve((F * ww[:, None]).T @ F + np.eye(rank), (F * ww[:, None]).T @ rr)
    Xi = np.hstack([X, (U[tid][:, :, None] * OD[:, None, :]).reshape(len(X), -1)])
    B = wls(Xi[keep], y[keep], w[keep])
    return {hd: np.exp(Xi[i] @ B.T) for hd, i in tgt.items()}


def group_nbr(xy, src_mask, gid, ng, bw, CP):
    """Mean neighbourhood CP10k per group (all cells; neighbours from src_mask), ng x G."""
    acc = np.zeros((ng, len(xy)), np.float32)
    cnt = np.bincount(gid, minlength=ng)
    for s0 in range(0, len(xy), 100000):
        r = np.arange(s0, min(len(xy), s0 + 100000))
        Wc = niche_rows(xy, r, src_mask, bw)
        acc += (
            sp.csr_matrix((np.ones(len(r), np.float32), (gid[r], np.arange(len(r)))), shape=(ng, len(r))) @ Wc
        ).toarray()
    return np.asarray(CP.T @ acc.T).T / np.maximum(cnt, 1)[:, None]


def cf_mean(base_ln, shift):
    """Per-cell counterfactual: log1p CP10k + shift -> expm1, clip, renormalise to 1e4 (Cellina saves normalised cf);
    pseudobulk mean."""
    x = np.clip(np.expm1(base_ln + shift), 0, None)
    x = x / (x.sum(1, keepdims=True) + 1e-8) * 1e4
    return x.mean(0)


def score(mc, mt, mf, bleed):
    gt = np.log2((mt + 1e-6) / (mc + 1e-6))
    cf = np.log2((mf + 1e-6) / (mc + 1e-6))
    deg = np.argsort(-np.abs(gt))[:NDEG]
    inter = np.intersect1d(deg, np.argsort(-np.abs(cf))[:NDEG])
    nb = deg[~bleed[deg]]
    g2 = np.log2((mt + PSC) / (mc + PSC))
    c2 = np.log2((mf + PSC) / (mc + PSC))
    d2 = np.argsort(-np.abs(g2))[:NDEG]
    pr = lambda a, b: pearsonr(a, b)[0] if np.std(b) > 0 else np.nan
    return dict(
        pearson=pr(gt[deg], cf[deg]),
        spearman=spearmanr(gt[deg], cf[deg])[0],
        precision=len(inter) / NDEG,
        direction_match_k=float((np.sign(gt[inter]) == np.sign(cf[inter])).sum()) / NDEG,
        mse_lfc=float(np.mean((gt[deg] - cf[deg]) ** 2)),
        pearson_clean=pr(g2[d2], c2[d2]),
        pearson_nonbleed=pr(gt[nb], cf[nb]) if len(nb) >= 10 else np.nan,
        n_nonbleed=len(nb),
        n_zero_ctrl=int((mc[deg] == 0).sum()),
    )


# ---------- one slide ----------
def run_slide(sid, data, bw, ctrl_doms, hds, htypes, out):
    out = Path(out)
    C, xy, ct, dom = data
    n, G = C.shape
    rng = np.random.default_rng(0)
    rows, lfc = [], {}
    lib = np.asarray(C.sum(1)).ravel()
    CP = sp.csr_matrix(sp.diags(1e4 / (lib + 1e-8)) @ C)
    LN = CP.copy()
    LN.data = np.log1p(LN.data)  # CP10k and log1p CP10k, sparse
    types = sorted(set(ct))
    doms = sorted(set(dom))
    tmeans = {t: nmean(C, np.flatnonzero(ct == t)) for t in types}
    gs = group_stats(C, ct, dom, lib)
    gkey = gs[0]
    for h in htypes:
        t0 = time.time()
        is_h = ct == h
        test = is_h & np.isin(dom, hds)
        train = ~test
        ctrl = np.flatnonzero(is_h & np.isin(dom, ctrl_doms))
        hT_all = np.flatnonzero(is_h & train)
        hT = rng.choice(hT_all, min(CAP_H, len(hT_all)), replace=False)
        base = rng.choice(ctrl, min(CAP_BASE, len(ctrl)), replace=False)
        other = np.flatnonzero(train & ~is_h)
        pool = np.union1d(rng.choice(other, min(CAP_POOL, len(other)), replace=False), hT)
        # PCA of log1p CP10k on training cells; scores for every cell (neighbours can be any non-h cell)
        fitr = rng.choice(np.flatnonzero(train), min(50000, train.sum()), replace=False)
        pca = PCA(NPC, svd_solver="randomized", random_state=0).fit(LN[fitr].toarray())
        Z = (np.asarray(LN @ pca.components_.T) - pca.mean_ @ pca.components_.T).astype(np.float32)
        bleed = np.max([tmeans[t] for t in types if t != h], 0) > 2 * np.maximum(tmeans[h], 0.1)
        mc = nmean(C, ctrl)
        mh = nmean(C, hT_all)
        mc_eb = (len(ctrl) * mc + KAPPA * mh) / (len(ctrl) + KAPPA)
        anch = {hd: np.flatnonzero(is_h & (dom == hd)) for hd in hds}
        anc = {hd: rng.choice(v, min(CAP_BASE, len(v)), replace=False) for hd, v in anch.items()}
        R = np.unique(np.concatenate([hT, base, pool, ctrl] + list(anc.values())))
        ix = lambda a: np.searchsorted(R, a)
        W = niche_rows(xy, R, ~is_h, bw)
        NZ = np.asarray(W @ Z)
        # L2 beta: h's own CP10k on its neighbourhood CP10k (per gene OLS slope, clipped [0, 1])
        Wh = W[ix(hT)]
        Y2 = CP[hT].toarray()
        X2 = np.vstack([(Wh[s0 : s0 + 5000] @ CP).toarray() for s0 in range(0, len(hT), 5000)])
        Xc = X2 - X2.mean(0)
        vx = (Xc**2).sum(0)
        beta = np.clip(np.where(vx > 0, (Xc * (Y2 - Y2.mean(0))).sum(0) / np.maximum(vx, 1e-12), 0), 0, 1)
        del Y2, X2, Xc
        nu_c = np.asarray(W[ix(ctrl)].mean(0) @ CP).ravel()
        # L3: h-niche ridge on h's training cells
        mu3, sd3 = NZ[ix(hT)].mean(0), NZ[ix(hT)].std(0) + 1e-6
        S3 = (NZ - mu3) / sd3
        B3, a3 = ridge(S3[ix(hT)], LN[hT], np.ones(len(hT), bool), rng)
        # L5: pooled ridge on all training cells: type + domain one-hots, Nz, h x Nz
        mu5, sd5 = NZ[ix(pool)].mean(0), NZ[ix(pool)].std(0) + 1e-6
        S5 = (NZ - mu5) / sd5

        def feat5(r, d, s):
            T = (ct[r][:, None] == np.array(types)).astype(np.float32)
            D = (np.asarray(d)[:, None] == np.array(doms)).astype(np.float32)
            hh = (ct[r] == h)[:, None]
            return np.hstack([T, D, s, s * hh])

        B5, a5 = ridge(feat5(pool, dom[pool], S5[ix(pool)]), LN[pool], ct[pool] == h, rng)
        base_ln = lognorm_dense(C, base)
        Sb3, Sb5 = S3[ix(base)], S5[ix(base)]
        F_own = feat5(base, dom[base], Sb5)
        tree4 = cKDTree(S3[ix(hT)])
        p6 = twoway(gs, h, hds)
        p8 = twoway(gs, h, hds, rank=1)
        gid = np.array([gkey.index(k) for k in zip(ct, dom)])
        p7 = twoway(gs, h, hds, nu=group_nbr(xy, ~is_h, gid, len(gkey), bw, CP))
        for hd in hds:
            task = f"{h}_{hd}"
            tg = anch[hd]
            mt = nmean(C, tg)
            preds = {}
            # published baseline (eval_loo.py: training-set shift of the held-out domain applied to raw control counts)
            rest, tgt = np.flatnonzero(train & (dom != hd)), np.flatnonzero(train & (dom == hd))
            delta = np.log2((np.asarray(C[tgt].mean(0)).ravel() + 1e-8) / (np.asarray(C[rest].mean(0)).ravel() + 1e-8))
            acc = np.zeros(G)
            for s0 in range(0, len(ctrl), 20000):
                X = np.clip((C[ctrl[s0 : s0 + 20000]].toarray() + 1) * 2**delta - 1, 0, None)
                acc += (X / (X.sum(1, keepdims=True) + 1e-8) * 1e4).sum(0)
            preds["baseline_repro"] = acc / len(ctrl)
            # L1 shared shift of the other types
            oth = [
                t
                for t in types
                if t != h and ((ct == t) & np.isin(dom, ctrl_doms)).sum() >= 5 and ((ct == t) & (dom == hd)).sum() >= 5
            ]
            wt = np.array([((ct == t) & (dom == hd)).sum() for t in oth], float)
            sh = np.average(
                [
                    np.log1p(nmean(C, np.flatnonzero((ct == t) & (dom == hd))))
                    - np.log1p(nmean(C, np.flatnonzero((ct == t) & np.isin(dom, ctrl_doms))))
                    for t in oth
                ],
                0,
                wt,
            )
            preds["L1_shared_shift"] = np.clip(np.expm1(np.log1p(mc_eb) + sh), 0, None)
            preds["L6_twoway_pseudobulk"] = p6[hd]
            preds["L7_twoway_contam"] = p7[hd]
            preds["L8_twoway_rank1"] = p8[hd]
            # L2 contamination
            nu_t = np.asarray(W[ix(anc[hd])].mean(0) @ CP).ravel()
            preds["L2_contamination"] = np.clip(mc_eb + beta * (nu_t - nu_c), 0, None)
            # L3 / L3null / L5: per control cell, a random anchor's niche vs its own (Cellina's anchor pairing)
            pa = rng.choice(anc[hd], len(base))
            preds["L3_h_niche_ridge"] = cf_mean(base_ln, (S3[ix(pa)] - Sb3) @ B3)
            pn = rng.choice(base, len(base))
            preds["L3null"] = cf_mean(base_ln, (S3[ix(pn)] - Sb3) @ B3)
            F_cf = feat5(pa, np.full(len(pa), hd), S5[ix(pa)])
            F_cf[:, : len(types)] = F_own[:, : len(types)]  # keep type h
            preds["L5_pooled_ridge"] = cf_mean(base_ln, (F_cf - F_own) @ B5)
            # L4 niche matching
            _, j = tree4.query(S3[ix(anc[hd])], k=KMATCH, workers=-1)
            cnt = np.bincount(j.ravel(), minlength=len(hT)).astype(float)
            preds["L4_niche_match"] = nmean(C, hT[cnt > 0], cnt[cnt > 0])
            for m, mf in preds.items():
                rows.append(
                    dict(
                        sid=sid,
                        holdout_celltype=task,
                        model_name=m,
                        n_ctrl=len(ctrl),
                        n_target=len(tg),
                        alpha3=a3,
                        alpha5=a5,
                        **score(mc, mt, mf, bleed),
                    )
                )
            lfc[task] = dict(mc=mc, mt=mt, bleed=bleed, **preds)
        log(
            sid,
            h,
            f"{time.time() - t0:.0f}s",
            " ".join(f"{r['model_name'][:6]} {r['pearson']:.3f}" for r in rows[-len(preds) * len(hds) :][: len(preds)]),
        )
    pd.DataFrame(rows).to_csv(out / f"{sid}.csv", index=False)
    np.savez_compressed(out / f"{sid}_lfc.npz", **{f"{t}|{k}": v for t, d in lfc.items() for k, v in d.items()})


# ---------- v3: two-way floor variants only (contamination bandwidth, similarity-weighted shift, EB groups) ----------
V3_BW_UM, V3_SIM_POW, V3_KAPPA = (100.0, 30.0, 15.0), 2.0, 10.0


def type_similarity(gs, h, hds, ctrl_doms, power=V3_SIM_POW, floor=0.05, centred=True):
    """Per-group multiplier: (correlation of the type's control-domain log CP10k profile with h's, after removing the
    mean
    over types (centred) or not, clipped at 0) ** power, floored; h's own groups 1. Training groups only ((h, held-out)
    excluded). Centred correlations of distinct GSE1 types are all <= 0 (crc_232), so centred weights are uniform."""
    keys, n, M, T, L = gs
    t_ = np.array([k[0] for k in keys])
    d_ = np.array([k[1] for k in keys])
    types = sorted(set(t_))
    prof = {
        t: np.log1p(np.average(M[(t_ == t) & np.isin(d_, ctrl_doms)], 0, n[(t_ == t) & np.isin(d_, ctrl_doms)]))
        for t in types
        if ((t_ == t) & np.isin(d_, ctrl_doms)).any()
    }
    mu = np.mean(list(prof.values()), 0) if centred else 0.0
    s = {t: max(np.corrcoef(prof[t] - mu, prof[h] - mu)[0, 1], 0.0) ** power if t != h else 1.0 for t in prof}
    return np.array([max(s.get(t, floor), floor) for t in t_])


def eb_groups(gs, h, hds, kappa=V3_KAPPA):
    """Group CP10k means shrunk toward their type's mean over its training groups: (n M + kappa M_type) / (n +
    kappa)."""
    keys, n, M, T, L = gs
    t_ = np.array([k[0] for k in keys])
    d_ = np.array([k[1] for k in keys])
    tr = ~((t_ == h) & np.isin(d_, hds))
    M2 = M.copy()
    for t in set(t_):
        m = (t_ == t) & tr
        Mt = np.average(M[m], 0, n[m])
        i = t_ == t
        M2[i] = (n[i, None] * M[i] + kappa * Mt) / (n[i, None] + kappa)
    return keys, n, M2, T, L


def run_twoway_v3(sid, data, px_per_um, ctrl_doms, hds, htypes, out):
    """F4 v3 candidates on one slide: L6 (twoway), L7 at each V3_BW_UM (neighbourhood contamination term),
    L9 (L6 with similarity-weighted WLS), L9 + L7 at each bandwidth, L6_eb (EB-shrunk group means), L8 rank 1.
    Scored as run_slide (score: official / spike-free / non-bleed pearson) -> <out>/<sid>_v3.csv."""
    out = Path(out)
    C, xy, ct, dom = data
    n, G = C.shape
    rows = []
    lib = np.asarray(C.sum(1)).ravel()
    CP = sp.csr_matrix(sp.diags(1e4 / (lib + 1e-8)) @ C)
    types = sorted(set(ct))
    tmeans = {t: nmean(C, np.flatnonzero(ct == t)) for t in types}
    gs = group_stats(C, ct, dom, lib)
    gkey = gs[0]
    gid = np.array([gkey.index(k) for k in zip(ct, dom)])
    for h in htypes:
        t0 = time.time()
        is_h = ct == h
        ctrl = np.flatnonzero(is_h & np.isin(dom, ctrl_doms))
        bleed = np.max([tmeans[t] for t in types if t != h], 0) > 2 * np.maximum(tmeans[h], 0.1)
        mc = nmean(C, ctrl)
        nus = {bw: group_nbr(xy, ~is_h, gid, len(gkey), bw * px_per_um, CP) for bw in V3_BW_UM}
        sim = type_similarity(gs, h, hds, ctrl_doms)
        P = {
            "L6_twoway_pseudobulk": twoway(gs, h, hds),
            "L9_simweighted": twoway(gs, h, hds, gw=sim),
            "L6_eb": twoway(eb_groups(gs, h, hds), h, hds),
            "L8_twoway_rank1": twoway(gs, h, hds, rank=1),
        }
        for bw, nu in nus.items():
            P[f"L7_bw{bw:g}"] = twoway(gs, h, hds, nu=nu)
            P[f"L9+L7_bw{bw:g}"] = twoway(gs, h, hds, nu=nu, gw=sim)
        simu = type_similarity(
            gs, h, hds, ctrl_doms, power=8.0, centred=False
        )  # round 2: uncentred profiles (power 8 for contrast)
        P["L9u_simweighted"] = twoway(gs, h, hds, gw=simu)
        P["L9u+L7_bw15"] = twoway(gs, h, hds, nu=nus[15.0], gw=simu)
        for hd in hds:
            mt = nmean(C, np.flatnonzero(is_h & (dom == hd)))
            for m, p in P.items():
                rows.append(
                    dict(
                        sid=sid,
                        holdout_celltype=f"{h}_{hd}",
                        model_name=m,
                        n_ctrl=len(ctrl),
                        **score(mc, mt, p[hd], bleed),
                    )
                )
        log(
            sid,
            h,
            f"{time.time() - t0:.0f}s",
            " ".join(
                f"{r['model_name']} {r['pearson']:.3f}/{r['pearson_clean']:.3f}" for r in rows[-len(P) * len(hds) :]
            ),
        )
    pd.DataFrame(rows).to_csv(out / f"{sid}_v3.csv", index=False)
    return pd.DataFrame(rows)
