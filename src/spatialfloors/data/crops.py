"""Contiguous ~50k-cell crops chosen by code from coordinates + region labels only (PROTOCOL §3)."""
import numpy as np

N_TARGET, MIN_SHARE, PRESENT_MIN, GRID = 50000, .20, .02, 40


def square_window(xy, c, n):
    """The n cells nearest to centre c in L-inf distance: a contiguous axis-aligned square; returns indices,
    half-side."""
    d = np.max(np.abs(xy - c), 1); idx = np.argpartition(d, n - 1)[:n]; return np.sort(idx), float(d[idx].max())


def choose_crop(xy, regions=None, n_target=N_TARGET, grid=GRID, min_share=MIN_SHARE, present_min=PRESENT_MIN):
    """Slides with <= 1.1 n_target cells are kept whole. Else candidate centres on a grid x grid lattice over the 2-98%
    coordinate range. With regions: present = regions with >= present_min of the slide's labelled cells; pick the window
    maximising the minimum present-region share (ties: smaller half-side), and record whether it reaches min_share.
    Without regions: the smallest square holding n_target cells (densest contiguous tissue)."""
    xy = np.asarray(xy, float); n = len(xy)
    rec = dict(n_slide=n, n_target=n_target, rule=f"grid {grid}x{grid}, L-inf square, present >= {present_min}, min share {min_share}")
    lab = None if regions is None else np.asarray(regions, dtype=object)
    if lab is not None:
        ok = np.array([isinstance(v, str) and v != "" and v != "nan" for v in lab]); u, c = np.unique(lab[ok], return_counts=True)
        present = [r for r, k in zip(u, c) if k / ok.sum() >= present_min]; rec["slide_shares"] = {r: round(k / ok.sum(), 4) for r, k in zip(u, c)}
        rec["present"] = present
    if n <= 1.1 * n_target:
        idx = np.arange(n); rec.update(whole_slide=True, centre=None, half_side=None)
    else:
        qx, qy = np.quantile(xy[:, 0], [.02, .98]), np.quantile(xy[:, 1], [.02, .98]); best = None
        for cx in np.linspace(*qx, grid):
            for cy in np.linspace(*qy, grid):
                idx, h = square_window(xy, np.array([cx, cy]), n_target)
                if lab is not None:
                    li = lab[idx]; m = ok[idx]; sh = min(((li[m] == r).sum() / max(m.sum(), 1) for r in present), default=0.0); key = (-sh, h)
                else:
                    key = (h,)
                if best is None or key < best[0]: best = (key, idx, (cx, cy), h)
        _, idx, cen, h = best; rec.update(whole_slide=False, centre=[float(cen[0]), float(cen[1])], half_side=h)
    if lab is not None:
        li = lab[idx]; m = ok[idx]; sh = {r: round(float((li[m] == r).sum() / max(m.sum(), 1)), 4) for r in rec["present"]}
        rec.update(crop_shares=sh, min_share=min(sh.values()) if sh else None, meets_min_share=bool(sh) and min(sh.values()) >= min_share,
                   n_unlabelled_in_crop=int((~m).sum()))
    rec["n_crop"] = int(len(idx)); return idx, rec
