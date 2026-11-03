"""F3 stress panel (PLAN §13 / §14): profiles any niche -> receiver-state model from its block_cv fold predictions.
  agreement   top-10 PC mean Spearman (primary), PC1, top-3, median per-gene Spearman (predictions mapped back to z-scored
              genes; N2 row MP3: the 10-PC mean dilutes one-axis signals), PC-space R2
  detail      worst receiver; receivers above the shuffled-feature null + .01 (null run at the anchors)
  orthogonal  protein (IF) shift-null check, N2 libella.niche.validate (vendored floors.validate): per receiver,
              |corr| of u with the kernel mean of the neighbours' residual IF beyond the id mix, vs toroidal FOV shifts
              (a within-FOV permutation null gives false positives, N2 §55). u_true = the receiver's own state PC1 (is there
              niche coupling in an RNA-free modality; model-free), u_pred = the model's out-of-fold PC1 prediction (does
              what the model predicts from neighbour RNA carry protein-confirmed neighbourhood signal, not bleed or noise)."""
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, rankdata
from spatialfloors.floors import validate as V

MS = ("top10", "pc1", "top3", "gene_median", "r2")


def colspearman(A, B):
    """Spearman of each column of A with the same column of B (NaN where either is constant)."""
    ra, rb = rankdata(A, axis=0), rankdata(B, axis=0); ra -= ra.mean(0); rb -= rb.mean(0)
    with np.errstate(invalid="ignore", divide="ignore"): return (ra * rb).sum(0) / np.sqrt((ra ** 2).sum(0) * (rb ** 2).sum(0))


def receiver_panel(s, det):
    pc1, top3, gene, r2 = [], [], [], []
    for f in det.get("folds", []):
        P, T = f["pred"], f["true"]; rho = [spearmanr(P[:, j], T[:, j])[0] for j in range(T.shape[1])]
        pc1.append(rho[0]); top3.append(np.nanmean(rho[:3])); r2.append(1 - ((P - T) ** 2).sum() / ((T - T.mean(0)) ** 2).sum())
        gene.append(float(np.nanmedian(colspearman(P @ f["comp"], f["Yz"]))))
    return dict(top10=float(s), pc1=float(np.nanmean(pc1)), top3=float(np.nanmean(top3)), gene_median=float(np.nanmean(gene)), r2=float(np.nanmean(r2)))


def panel(out):
    """out: {receiver: (score, folds, per_pc, detail)} (score_floor / gat.score with full=True) -> per-receiver panel."""
    return {int(r): receiver_panel(v[0], v[3]) for r, v in out.items()}


def summary(per, null=None):
    s = {k: float(np.nanmean([x[k] for x in per.values()])) for k in MS}; s["worst_top10"] = float(np.nanmin([x["top10"] for x in per.values()]))
    if null: s["above_null"] = float(np.mean([per[r]["top10"] > null[r]["top10"] + .01 for r in per if r in null]))
    return s


def oof(D, out, r):
    """(cells, u_true, u_pred) of receiver r: scored test cells over all folds, PC1 truth and prediction."""
    rows = np.flatnonzero(D["ids"] == r); f = out[r][3]["folds"]
    return rows[np.concatenate([x["te"] for x in f])], np.concatenate([x["true"][:, 0] for x in f]), np.concatenate([x["pred"][:, 0] for x in f])


class Protein:
    """IF field of one unit: residual IF (FOV, log library regressed out), observed field and nrep toroidal-shift fields."""
    def __init__(self, D, P, lib, nrep=50, seed=0):
        self.D = D; self.Pr = V.residual_protein(P, lib, D["fov"]); pos, fov = D["xy"], D["fov"]; sp_ = float(np.median(cKDTree(pos).query(pos, 2)[0][:, 1]))
        f = lambda q: V.field(pos, fov, D["lin"], D["ids"], self.Pr, D["K"], q, sp_); rng = np.random.default_rng(seed)
        self.obs = f(pos); self.null = [f(V.shift(pos, fov, sp_, rng)) for _ in range(nrep)]

    def __call__(self, out):
        """Per receiver: p_shift of u_true and u_pred; summary = share of receivers with p < .05 and the mean max|r|."""
        res = {}
        for r in out:
            c, ut, up = oof(self.D, out, r); Pn0, T0, w0 = self.obs; ok = w0[c] >= V.WMIN; c, ut, up = c[ok], ut[ok], up[ok]; C = np.c_[np.ones(len(c)), T0[c]]
            for nm, u in (("true", ut), ("pred", up)):
                o = V._stat(u, C, Pn0, T0, c)[0]; nl = []
                for Pn, Tq, w in self.null:
                    v = w[c] >= V.WMIN; nl.append(V._stat(u[v], C[v], Pn, Tq, c[v])[0])
                res.setdefault(int(r), {})[nm] = dict(max_abs=float(o), p=float((1 + (np.array(nl) >= o).sum()) / (1 + len(nl))))
        return dict(per_receiver=res, prot_true_sig=float(np.mean([v["true"]["p"] < .05 for v in res.values()])),
                    prot_pred_sig=float(np.mean([v["pred"]["p"] < .05 for v in res.values()])), prot_pred_r=float(np.mean([v["pred"]["max_abs"] for v in res.values()])))
