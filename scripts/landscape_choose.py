"""Cross-unit floor config choice (PLAN §14 choice rule): per grid config, the primary as a share of each unit's best and
the worst share over development units; plateau count = units where the config is within EPS of the unit's best.
Generalises f2_choose.py. Usage: python scripts/landscape_choose.py <fam> <primary> [n] -> prints the top n configs."""
import glob, json, sys
fam, prim, n = sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 10; EPS = .02
R = {}
for f in sorted(glob.glob(f"runs/landscape/{fam}_*.jsonl")):
    rows = [r for r in map(json.loads, open(f)) if r.get("part") == "grid" and r.get(prim) == r.get(prim)]
    if rows: R[rows[0]["unit"]] = {r["id"]: r[prim] for r in rows}
ids = set.intersection(*(set(v) for v in R.values())); best = {u: max(v.values()) for u, v in R.items()}
key = lambda i: min(R[u][i] / best[u] for u in R)
for i in sorted(ids, key=lambda i: -key(i))[:n]:
    print(f"{i} worst share {key(i):.3f} plateau {sum(best[u] - R[u][i] <= EPS for u in R)}/{len(R)} |", " ".join(f"{u[-5:]} {R[u][i]:.3f}" for u in R))
