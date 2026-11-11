"""Figure: the 12 registered pairs of round 1 on one axis. x = D / margin of the pair's family (D = learned - floor per
unit; estimate and one-sided 95% bounds from the tier-1 report), so the non-inferiority limit is +1 for every pair.
Reads only results/reports/P5_corrected.json.
Usage (env floors): python scripts/fig_forest.py <out.png>   (also writes the same name with .pdf)"""

import json, sys
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

PUB = Path(__file__).resolve().parents[1]
P = json.loads((PUB / "results/reports/P5_corrected.json").read_text())["pairs"]
ROWS = [  # family, pair key, label
    ("F1 cell identity", "F1|scvi", "scVI"),
    ("F1 cell identity", "F1|nicheformer", "Nicheformer (zero-shot)"),
    ("F2 spatial domains", "F2|cellcharter", "CellCharter"),
    ("F2 spatial domains", "F2|novae", "Novae (mouse, zero-shot)"),
    ("F3 context to receiver state", "F3|lightgbm", "LightGBM"),
    ("F3 context to receiver state", "F3|ncem", "NCEM"),
    ("F4 counterfactual", "F4_pearson|cellina", "Cellina, official"),
    ("F4 counterfactual", "F4_pearson|cpa", "CPA, official"),
    ("F4 counterfactual", "F4_pearson_clean|cellina", "Cellina, spike-free"),
    ("F4 counterfactual", "F4_pearson_clean|cpa", "CPA, spike-free"),
    ("F5 label transfer", "F5|scanvi", "scANVI"),
    ("F5 label transfer", "F5|tangram", "Tangram"),
]
INK, INK2, MUTED, GRID, AXIS, SURF, SERIES = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb", "#2a78d6"
plt.rcParams.update(
    {"font.family": "sans-serif", "font.sans-serif": ["Helvetica Neue", "Arial", "DejaVu Sans"], "font.size": 9}
)

fig, ax = plt.subplots(figsize=(7.2, 4.9), dpi=300)
fig.patch.set_facecolor(SURF)
ax.set_facecolor(SURF)
y, ys, fam_y, prev = 0.0, [], [], None
for fam, key, lab in ROWS:
    if fam != prev:
        y += 0.9 if prev else 0.0
        fam_y.append((y, fam))
        y += 0.8
        prev = fam
    ys.append(y)
    y += 1.0
XMIN, XMAX = -10.5, 2.0
TR = ax.get_yaxis_transform()  # x in axes fraction, y in data units
COLS = (1.19, 1.39, 1.555, 1.66)  # right edges of the table columns
for x in range(-10, 1, 2):
    ax.axvline(x, color=GRID, lw=0.6, zorder=0)
ax.axvline(0, color=AXIS, lw=1.0, zorder=1)
ax.axvline(1, color=INK2, lw=1.0, zorder=1)
for (fam, key, lab), yy in zip(ROWS, ys):
    v = P[key]
    m = v["margin"]
    lo, est, hi = v["lower"] / m, v["est"] / m, v["upper"] / m
    ax.plot([max(lo, XMIN), hi], [yy, yy], color=SERIES, lw=2, solid_capstyle="round", zorder=2)
    ax.plot([est], [yy], "o", ms=7.5, mfc=SERIES, mec=SURF, mew=1.6, zorder=3)
    ax.text(XMIN - 0.25, yy, lab, ha="right", va="center", color=INK)
    cells = (f"{v['est']:+.4f}", f"{v['upper']:+.4f}", f"{m:.4f}".lstrip("0"), str(v["n_units"]))
    for cx, txt in zip(COLS, cells):
        ax.text(cx, yy, txt, ha="right", va="center", color=INK2, family="monospace", fontsize=7.5, transform=TR)
for yy, fam in fam_y:
    ax.text(XMIN - 0.25, yy, fam, ha="right", va="center", color=MUTED, fontsize=8)
top = ys[0] - 1.9
for cx, t in zip(COLS, ("D", "upper", "margin", "units")):
    ax.text(cx, top + 0.8, t, ha="right", va="center", color=MUTED, fontsize=7.5, transform=TR)
ax.text(1.15, top, "non-inferiority limit", ha="left", va="center", color=INK2, fontsize=7.5)
ax.text(-0.15, top, "no difference", ha="right", va="center", color=MUTED, fontsize=7.5)
ax.set_ylim(ys[-1] + 0.8, top - 0.9)
ax.set_xlim(XMIN, XMAX)
ax.set_xticks(list(range(-10, 1, 2)) + [1])
ax.set_xticklabels([str(x) for x in range(-10, 1, 2)] + ["+1"])
ax.set_yticks([])
ax.tick_params(axis="x", colors=MUTED, length=0, labelsize=8)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_xlabel("Learned model minus floor, in units of the pre-registered margin", color=INK2, fontsize=8.5, labelpad=6)
fig.text(
    0.012,
    0.975,
    "Round 1: every upper bound is left of the non-inferiority limit",
    ha="left",
    va="top",
    color=INK,
    fontsize=10.5,
    fontweight="bold",
)
fig.text(
    0.012,
    0.925,
    "Twelve pre-registered pairs. Dot: mean difference over units. Line: one-sided 95% bounds. Left of 0: the floor scores higher.",
    ha="left",
    va="top",
    color=INK2,
    fontsize=7.8,
)
fig.subplots_adjust(left=0.25, right=0.69, top=0.86, bottom=0.11)
out = Path(sys.argv[1])
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, facecolor=SURF)
fig.savefig(out.with_suffix(".pdf"), facecolor=SURF)
print("wrote", out, out.with_suffix(".pdf"))
