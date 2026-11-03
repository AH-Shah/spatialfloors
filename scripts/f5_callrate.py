"""F5 call-rate check (PLAN §14 (d)): does held-out-marker AUROC reward abstention? Per unit: InSituType (P3 profile / NegPrb)
at posterior thresholds 0-.95, logistic regression, scANVI and Tangram with their call shares; over the floor grid, Spearman
of AUROC with the 'other' (uncalled) share, and the AUROC of the floor at no threshold vs the learned models.
Usage: python scripts/f5_callrate.py [units]"""
import glob, json, sys
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
U = sys.argv[1:] or [Path(f).stem[3:] for f in sorted(glob.glob("runs/landscape/f5_*.jsonl"))]
M = ("heldout_auroc", "worst_class", "share_other", "gated_macro_f1", "excess_over_threshold", "depth_strat_mean")
f = lambda v: "-" if v is None else f"{v:.3f}"
for u in U:
    R = [json.loads(l) for l in open(f"runs/landscape/f5_{u}.jsonl")]; print(f"## {u}\n| row | " + " | ".join(M) + " |")
    sel = [r for r in R if r["id"].startswith("ist|all|") and r["id"].endswith("|none")] + [r for r in R if r.get("part") == "learned" and r.get("perturb") == "none"]
    sel += [max((r for r in R if r.get("part") == "logreg"), key=lambda r: r["heldout_auroc"])]
    for r in sel: print(f"| {r['id'][:40]} | " + " | ".join(f(r.get(m)) for m in M) + " |")
    G = [r for r in R if r.get("part") in ("ist", "logreg")]
    print(f"grid rho(AUROC, share_other) {spearmanr([r['heldout_auroc'] for r in G], [r['share_other'] for r in G])[0]:+.2f}; ist t0 only: best "
          f"{max(r['heldout_auroc'] for r in G if r.get('threshold') == 0):.3f}\n")
