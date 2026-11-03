"""F3 niche -> receiver state: quadrant CV with edge buffer; train-only PCA / scalers asserted by index (SPEC §F3)."""
import os
import numpy as np
from scipy.stats import spearmanr

NPC = 10
CLIP_X = os.environ.get("SF_F3_CLIP", "0") == "1"  # dev A/B (PLAN §12): clip scaled test features to the training range


def quadrant_blocks(xy):
    """4 contiguous blocks split at median x and median y -> ids 0..3, plus the split point."""
    mx, my = np.median(xy[:, 0]), np.median(xy[:, 1])
    return (xy[:, 0] > mx).astype(int) + 2 * (xy[:, 1] > my).astype(int), (mx, my)


def edge_buffer(xy, split, receptive_field):
    """True for cells within 2 x receptive_field of a block edge (the median lines): not scored."""
    mx, my = split; d = np.minimum(np.abs(xy[:, 0] - mx), np.abs(xy[:, 1] - my)); return d < 2 * receptive_field


def restrict_graph_to_blocks(W, blocks):
    """Keep only edges i-j with blocks[i] == blocks[j] (sparse adjacency W)."""
    W = W.tocoo(); m = blocks[W.row] == blocks[W.col]
    import scipy.sparse as sp
    return sp.csr_matrix((W.data[m], (W.row[m], W.col[m])), shape=W.shape)


class FitLog:
    """Records the row indices every fit sees; block_cv asserts none is a test row."""
    def __init__(self): self.rows = []
    def __call__(self, idx): self.rows.append(np.asarray(idx)); return idx


def block_cv(X, Y, xy, fit_predict, receptive_field, npc=NPC, fitlog=None, blocks=None, split=None, covars=None, detail=None):
    """X features (cells x p, already built within blocks) or a callable X(q, train_rows) -> features of all rows (for
    fold-dependent feature fits, which must call fitlog on the rows they fit on); Y receiver state (cells x q, bleed-safe).
    blocks / split: the crop's quadrants (default: from xy itself; pass the crop's when the rows are one receiver id).
    Per held-out quadrant, on training rows only: Y regressed on covars (+ intercept) then z-scored, PCA top npc, scaler on
    X; fit_predict(Xtr, Ytr_pcs, Xte) -> predicted PCs. Score = mean over folds of the mean Spearman over PCs on scored test
    rows (not in the edge buffer). detail (dict, optional) receives per-fold per-PC Spearman under 'per_pc'."""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    if blocks is None: blocks, split = quadrant_blocks(xy)
    buf = edge_buffer(xy, split, receptive_field); log = fitlog or FitLog(); out = []
    for q in range(4):
        tr = np.flatnonzero(blocks != q); te = np.flatnonzero((blocks == q) & ~buf)
        if len(te) < 20 or len(tr) < 20: continue
        Xq = X(q, tr) if callable(X) else X; Yq = Y
        live = Xq[tr].std(0) > 1e-6; Xq = Xq[:, live]  # near-constant training columns (e.g. an id only seen far away) would scale to inf
        if covars is not None:
            C = np.c_[np.ones(len(Y)), covars]; Yq = Y - C @ np.linalg.lstsq(C[log(tr)], Y[tr], rcond=None)[0]
        sy = StandardScaler().fit(Yq[log(tr)]); Yq = sy.transform(Yq)
        sx = StandardScaler().fit(Xq[log(tr)]); pca = PCA(min(npc, Y.shape[1]), random_state=0).fit(Yq[log(tr)])
        Ptr, Pte = pca.transform(Yq[tr]), pca.transform(Yq[te])
        Xtr_, Xte_ = sx.transform(Xq[tr]), sx.transform(Xq[te])
        if CLIP_X: Xte_ = np.clip(Xte_, Xtr_.min(0), Xtr_.max(0))  # a rare id seen in a few training cells scales to |z| ~ 1e3 in test
        pred = (fit_predict(Xtr_, Ptr, Xte_, tr, te, q) if getattr(fit_predict, "wants_rows", False) else fit_predict(Xtr_, Ptr, Xte_))  # graph models also get row ids and the fold
        rho = [spearmanr(pred[:, j], Pte[:, j])[0] for j in range(Pte.shape[1])]; out.append(np.nanmean(rho))
        if detail is not None: detail.setdefault("per_pc", []).append(rho); detail.setdefault("evr", []).append(pca.explained_variance_ratio_.tolist())
        if detail is not None and detail.get("full"):  # stress panel (PLAN §13): fold predictions in PC and z-scored gene space
            detail.setdefault("folds", []).append(dict(te=te, pred=np.asarray(pred), true=Pte, Yz=Yq[te], comp=pca.components_))
        for r in log.rows: assert not np.isin(r, te).any() and not np.isin(r, np.flatnonzero(blocks == q)).any(), "fit saw held-out quadrant"
        log.rows.clear()
    return float(np.mean(out)) if out else float("nan"), out
