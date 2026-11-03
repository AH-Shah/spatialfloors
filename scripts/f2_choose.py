"""F2 floor config choice (PLAN §14 choice rule, F2 decision row 359): over the grid, each config's balanced ARI as a share of
the unit's best, worst over development units; also its seed stability and min region recall. Prints the top configs.
Usage: python scripts/f2_choose.py <unit> <unit> ..."""
import json, sys
import numpy as np
U = sys.argv[1:]; R = {u: {r["id"]: r for r in map(json.loads, open(f"runs/landscape/f2_{u}.jsonl")) if r.get("part") == "grid"} for u in U}
ids = set.intersection(*(set(v) for v in R.values())); best = {u: max(r["ari_bal"] for r in R[u].values()) for u in U}
rows = sorted(ids, key=lambda i: -min(R[u][i]["ari_bal"] / best[u] for u in U))
for i in rows[:12]: print(i, f"worst share {min(R[u][i]['ari_bal'] / best[u] for u in U):.3f} |", " | ".join(f"{u[:8]} bal {R[u][i]['ari_bal']:.3f} ari {R[u][i]['ari']:.3f} stab {R[u][i]['seed_stability']:.2f} minrec {R[u][i]['min_region_recall']:.2f}" for u in U))
