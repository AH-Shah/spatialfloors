# Vendored VERBATIM from libella/src/libella/niche/validate.py (N2 trusted baselines, CALCULUS §55); only change: the
# ELL / KNN import points at floors.bleedsafe (identical constants).
"""Protein validation of niche programs with a spatial-shift null (CALCULUS §55; independent IF modality, no RNA bleed).

Statistic per receiver id: max over IF channels of |corr(res_T(u), res_T(Pn))|, u = program axis 1, Pn = kernel mean (l 5 spacings,
200 NN, same FOV, other lineages) of FOV/depth-residualised IF, both residualised on [1, T] at their own location.
Null: each FOV's protein field is shifted toroidally (wrap in the FOV box, |s| >= 15 spacings) and Pn, T are re-evaluated at the shifted
points; both fields keep their spatial autocorrelation (a within-FOV permutation null gives false positives, §55 NS1)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .bleedsafe import ELL, KNN

SMIN, WMIN = 15.0, .5


def residual_protein(P, lib, fov):
    """log1p IF channels with FOV and log library regressed out."""
    X = np.c_[np.ones(len(P)), np.log1p(lib), pd.get_dummies(fov).values.astype(float)]; P = np.log1p(np.clip(P, 0, None)); return P - X @ np.linalg.lstsq(X, P, rcond=None)[0]


def field(pos, fov, lin, ids, Pr, K, q, sp):
    """kernel mean of Pr and id mix T at query points q (one per cell, in its own FOV), excluding cells of the cell's own lineage."""
    n = len(q); Pn = np.zeros((n, Pr.shape[1])); T = np.zeros((n, K)); ws = np.zeros(n)
    for f in np.unique(fov):
        m = np.flatnonzero(fov == f); k = min(KNN, len(m)); d, nb = cKDTree(pos[m]).query(q[m], k); d, nb = d.reshape(len(m), k) / sp, m[nb.reshape(len(m), k)]
        w = np.exp(-d ** 2 / (2 * ELL ** 2)) * (lin[nb] != lin[m, None]); s = w.sum(1); ws[m] = s; s = np.maximum(s, 1e-9)[:, None]
        Pn[m] = np.einsum("nk,nkd->nd", w, Pr[nb]) / s; t = np.zeros((len(m), K)); np.add.at(t, (np.repeat(np.arange(len(m)), k), ids[nb].ravel()), w.ravel()); T[m] = t / s
    return Pn, T, ws


def shift(pos, fov, sp, rng):
    q = pos.copy()
    for f in np.unique(fov):
        m = fov == f; lo = pos[m].min(0); L = pos[m].max(0) + 1e-6 - lo
        if np.hypot(*(L / 2)) < SMIN * sp: s = L / 2  # crash fix 2026-11-01: FOV sliver at a crop edge (gse1_231 FOV 9: 10 cells,
        else:                                         # bruker_breast 314: 18) can never reach SMIN; the loop spun forever. Largest shift instead.
          while True:
            s = rng.uniform(0, L)
            if np.hypot(*np.minimum(s, L - s)) >= SMIN * sp: break
        q[m] = lo + np.mod(pos[m] - lo + s, L)
    return q


def _stat(u, C, Pn, T, i):
    ur = u - C @ np.linalg.lstsq(C, u, rcond=None)[0]; Cq = np.c_[np.ones(len(i)), T[i]]; pn = Pn[i] - Cq @ np.linalg.lstsq(Cq, Pn[i], rcond=None)[0]
    r = np.array([np.corrcoef(ur, pn[:, j])[0, 1] for j in range(pn.shape[1])]); return np.nanmax(np.abs(r)), r


def protein_shift(U, ids, lineage, pos, fov, Pr, receivers=None, nrep=100, seed=0):
    """U: n x 3 program scores (0 for non-receivers); Pr: residualised IF (residual_protein). Returns per-id max|r|, shift95, p."""
    ids = np.asarray(ids); K = int(ids.max()) + 1; fov = np.asarray(fov).astype(str); sp = float(np.median(cKDTree(pos).query(pos, 2)[0][:, 1]))
    receivers = receivers if receivers is not None else [r for r in range(K) if np.abs(U[ids == r, 0]).sum() > 0]; rng = np.random.default_rng(seed)
    Pn0, T0, w0 = field(pos, fov, lineage, ids, Pr, K, pos, sp); ok0 = w0 >= WMIN; obs, null, Cs = {}, {r: [] for r in receivers}, {}
    for r in receivers:
        i = np.flatnonzero((ids == r) & ok0); Cs[r] = (i, np.c_[np.ones(len(i)), T0[i]]); obs[r] = _stat(U[i, 0], Cs[r][1], Pn0, T0, i)
    for _ in range(nrep):
        Pn, Tq, w = field(pos, fov, lineage, ids, Pr, K, shift(pos, fov, sp, rng), sp); ok = w >= WMIN
        for r in receivers:
            i, C = Cs[r]; v = ok[i]; null[r].append(_stat(U[i[v], 0], C[v], Pn, Tq, i[v])[0])
    rows = []
    for r in receivers:
        s = np.array(null[r]); rows.append(dict(receiver=r, n=int((ids == r).sum()), max_abs=obs[r][0], shift95=np.percentile(s, 95), p_shift=(1 + (s >= obs[r][0]).sum()) / (1 + len(s)),
                                                **{f"r_ch{j}": v for j, v in enumerate(obs[r][1])}))
    return pd.DataFrame(rows)
