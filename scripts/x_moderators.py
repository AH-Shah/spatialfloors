"""Exploratory (PROTOCOL §6): unit size as a moderator of the tier-1 differences. Per test pair of F1, F2, F3, F5: Spearman
correlation between the unit's D (tier-1 report, per_unit) and log10 of the cells it was scored on (sum over the unit's crops
that have an output for the family; registration/units/units.tsv). F4 units are whole slides / sections and are not included.
No decision weight; written after the tier-1 report.
Usage (env floors): python scripts/x_moderators.py <runs dir> [tier-1 report name, default P5_corrected] -> <runs dir>/X_moderators.{json,md}"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr

R = Path(sys.argv[1]); NAME = sys.argv[2] if len(sys.argv) > 2 else "P5_corrected"
P = json.loads((R / f"{NAME}.json").read_text()); T = pd.read_csv("registration/units/units.tsv", sep="\t"); T = T[T.role == "confirmation"]
FILE = {"F1": "f1_{}_expr_embed_pearson_none_pc50.json", "F2": "labels/f2_{}_floor.npz", "F3": "f3_{}_ridge.json", "F5": "f5_{}_floors.json"}
out = dict(tier1=NAME, pairs={})
for k, v in P["pairs"].items():
    fam = k.split("|")[0]
    if fam not in FILE: continue
    scored = T[[(R / FILE[fam].format(c)).exists() for c in T.crop]]  # tier-1 unit key: GSE303070 donor, otherwise the crop
    size = scored.assign(key=[u if c.startswith("gse303070") else c for c, u in zip(scored.crop, scored.unit)]).groupby("key").n_cells.sum()
    units = [u for u in v["per_unit"] if u in size.index]; D = np.array([v["per_unit"][u] for u in units]); n = np.log10(size[units].values.astype(float))
    rho, p = spearmanr(D, n) if len(units) > 3 else (np.nan, np.nan)
    out["pairs"][k] = dict(n_units=len(units), spearman_rho=float(rho), p_two_sided=float(p), cells_min=int(size[units].min()), cells_max=int(size[units].max()))
(R / "X_moderators.json").write_text(json.dumps(out, indent=1))
md = [f"Unit size as a moderator (exploratory). Spearman of per-unit D ({NAME}) vs log10 scored cells.", "", "| pair | units | cells (min-max) | rho | p (two-sided) |", "|---|---|---|---|---|"]
md += [f"| {k} | {v['n_units']} | {v['cells_min']:,}-{v['cells_max']:,} | {v['spearman_rho']:+.2f} | {v['p_two_sided']:.3f} |" for k, v in out["pairs"].items()]
(R / "X_moderators.md").write_text("\n".join(md) + "\n"); print(len(out["pairs"]), "pairs ->", R / "X_moderators.md")
