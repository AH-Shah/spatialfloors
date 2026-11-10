"""Summary of the F5 whole-slide diagnostic (exploratory; 2026-11-10): scANVI trained on the reference + 50,000 random
cells of each whole slide and applied to every cell (B_ratio, scripts/x_f5_scale.py) vs the floor on the same slide;
compared with the registered whole-slide run and the crop. Same test as the arm (one-sided 95% bounds, margin .0255).
Usage (env floors): python scripts/x_f5_scale_summary.py -> runs/xw/diag/X_f5_scale.{json,md}"""
import json
from pathlib import Path
import numpy as np
from spatialfloors import stats

R = Path("runs/xw/diag"); W = json.loads(Path("runs/xw/X_wholeslide.json").read_text())["pairs"]["F5|scanvi"]
S = ("110", "120", "210", "221", "222", "242"); rows = {}
for s in S:
    d = json.loads((R / f"X_f5_scale_{s}.json").read_text()); b = d["variants"]["B_ratio"]
    rows[s] = dict(n_cells=d["n_cells"], n_ref=d["n_ref"], ratio=d["n_cells"] / d["n_ref"], scanvi_registered=d["registered_whole"], scanvi_subsampled=b["heldout_auroc"],
                   floor=d["floor_whole"], D_registered=d["registered_whole"] - d["floor_whole"], D_subsampled=b["heldout_auroc"] - d["floor_whole"], D_crop=W["D_crop"][s],
                   A_epochs=d["variants"].get("A_epochs", {}).get("heldout_auroc"))
out = dict(per_slide=rows, tests={})
for k in ("D_registered", "D_subsampled", "D_crop"):
    D = np.array([rows[s][k] for s in S]); t = stats.family_test(D, None, margin=.0255)
    out["tests"][k] = dict(mean=float(D.mean()), lower=t["lower"], upper=t["upper"], noninferior=t["noninferior"])
md = ["# X_f5_scale: scANVI on whole slides with a crop-size query (exploratory diagnostic)", "",
      "scANVI (locked nl30_nb) trained on the reference + 50,000 random cells of the slide, then applied to every cell; floor = InSituType + spatial prior on the whole slide. On GSE1 221 a run with crop-equivalent epochs on all cells (69 / 10) scored .588, as the registered run.", "",
      "| Slide | Cells | Query : reference | scANVI registered | scANVI subsampled | Floor | D registered | D subsampled | D crop |", "|---|---|---|---|---|---|---|---|---|"]
md += [f"| {s} | {r['n_cells']:,} | {r['ratio']:.1f} | {r['scanvi_registered']:.3f} | {r['scanvi_subsampled']:.3f} | {r['floor']:.3f} | {r['D_registered']:+.3f} | {r['D_subsampled']:+.3f} | {r['D_crop']:+.3f} |" for s, r in rows.items()]
md += ["", "| D (6 slides) | Mean | One-sided 95% bounds | Non-inferior (margin .0255) |", "|---|---|---|---|"]
md += [f"| {k} | {v['mean']:+.4f} | {v['lower']:+.4f}, {v['upper']:+.4f} | {v['noninferior']} |" for k, v in out["tests"].items()]
(R / "X_f5_scale.json").write_text(json.dumps(out, indent=1, default=float)); (R / "X_f5_scale.md").write_text("\n".join(md) + "\n"); print("\n".join(md))
