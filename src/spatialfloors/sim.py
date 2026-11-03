"""Planted-signal toys for P2 leak tests and planted-signal checks.
niche_toy: libella/niche/checks.py `synth` (CALCULUS §60 NV1), same generator with FOVs tiled 3 x 2 so the four
quadrants are populated, returning counts and um coordinates
(grid spacing 1 -> SPACING_UM). 'bleed': each type has its own smooth state field, no cross-type coupling, counts +=
bleed x mean of the 8 nearest cells; 'couple': type-0 state follows its type-1 neighbours' kernel-mean state (beta);
'null': independent per-cell states, no bleed, no coupling. 'nonlinear': type-0 state follows type-1 state x type-2
presence (zero-mean product: no linear signal in either alone). 'quadratic': type-0 state follows |type-1 neighbours'
state| (no linear signal)."""
import numpy as np
import scipy.sparse as sp

SPACING_UM = 10.


def niche_toy(case, seed=0, nf=6, side=70, bleed=.3, beta=.8):
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(seed); pos, fov = [], []
    for f in range(nf):
        g = np.stack(np.meshgrid(np.arange(side), np.arange(side)), -1).reshape(-1, 2).astype(float) + rng.normal(0, .25, (side * side, 2)); pos.append(g + side * np.array([f % 3, f // 3])); fov += [f] * len(g)  # FOVs tile 3 x 2
    pos, fov = np.concatenate(pos), np.array(fov).astype(str); n = len(pos); d, nb = cKDTree(pos).query(pos, 40)
    w = np.exp(-d ** 2 / (2 * 3. ** 2)); ids = np.einsum("nk,nkc->nc", w, rng.normal(size=(n, 3))[nb]).argmax(1)
    if case == "bleed":
        F = np.einsum("nk,nkc->nc", w, rng.normal(size=(n, 3))[nb]); s = F[np.arange(n), ids]
    elif case == "null":
        s = rng.normal(size=n)
    else:
        s = rng.normal(size=n); m1 = (ids[nb] == 1) * w; s1 = (m1 * s[nb]).sum(1) / np.maximum(m1.sum(1), 1e-9)
        if case == "compo": s1 = m1.sum(1) / w.sum(1)  # receiver state linear in the local type-1 share: control for type-label models (NCEM)
        if case == "nonlinear":
            p2 = ((ids[nb] == 2) * w).sum(1) / w.sum(1); r0 = ids == 0  # standardised over receivers: no linear part in either factor
            z1 = (s1 - s1[r0].mean()) / s1[r0].std(); z2 = (p2 - p2[r0].mean()) / p2[r0].std(); s1 = z1 * z2
        if case == "quadratic":  # symmetric in the sender state: zero linear correlation
            r0 = ids == 0; z1 = (s1 - s1[r0].mean()) / s1[r0].std(); s1 = np.abs(z1)
        s = np.where(ids == 0, beta * (s1 - s1[ids == 0].mean()) / s1[ids == 0].std() + .6 * rng.normal(size=n), s)
    s = (s - s.mean()) / s.std(); lam = np.full((n, 90), .3); lam[:, 60:] = 1.0
    for t in range(3):
        m = ids == t; lam[m, t * 20:(t + 1) * 20] = 4.0; lam[m, t * 20:t * 20 + 10] *= np.exp(.5 * s[m])[:, None]
    if case == "bleed": lam = lam + bleed * lam[nb[:, 1:9]].mean(1)
    C = sp.csr_matrix(rng.poisson(lam).astype(np.float32))
    return C, pos * SPACING_UM, fov, ids


def domain_toy(seed=0, side=100, nf=4):
    """F2 toy: 2 domains (left / right half of each FOV) differing only in cell-type composition (types 0/1/2 at
    .6/.3/.1 vs .1/.3/.6), each type with fixed marker profile (same in both domains): domains are visible only through
    the neighbourhood. Returns counts, xy um, fov, domain labels."""
    rng = np.random.default_rng(seed); pos, fov, dom = [], [], []
    for f in range(nf):
        g = np.stack(np.meshgrid(np.arange(side), np.arange(side)), -1).reshape(-1, 2).astype(float) + rng.normal(0, .25, (side * side, 2))
        pos.append(g + side * np.array([f % 2, f // 2])); fov += [f] * len(g); dom.append((g[:, 0] >= side / 2).astype(int))
    pos, fov, dom = np.concatenate(pos), np.array(fov).astype(str), np.concatenate(dom)
    pr = np.where(dom[:, None] == 0, [.6, .3, .1], [.1, .3, .6]); t = (rng.random(len(dom))[:, None] > np.cumsum(pr, 1)).sum(1)
    lam = np.full((len(dom), 90), .3); lam[:, 60:] = 1.
    for k in range(3): lam[t == k, k * 20:(k + 1) * 20] = 4.
    return sp.csr_matrix(rng.poisson(lam).astype(np.float32)), pos * SPACING_UM, fov, dom
