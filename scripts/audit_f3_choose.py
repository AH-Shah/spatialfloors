"""Audit 2026-11-02: PLAN §14 choice rule for the F3 floor per species with the plateau tolerance at the F3 margin
(.0043) instead of .02 (5x the margin, which put 85-133 of 144 configs on the 'plateau' of each human unit).
Usage: python scripts/audit_f3_choose.py <eps> <unit ...> -> top configs by worst share, plateau count at eps."""
import json, sys
eps, units = float(sys.argv[1]), sys.argv[2:]; R = {}
for u in units:
    rows = [r for r in map(json.loads, open(f"runs/landscape/f3_{u}.jsonl")) if r.get("part") == "grid" and r.get("top10") == r.get("top10")]
    R[u] = {r["id"]: r["top10"] for r in rows}
ids = set.intersection(*(set(v) for v in R.values())); best = {u: max(v.values()) for u, v in R.items()}
lock = "grid|T+ctx|ell25|pc50|acv"; print("best", {u[-6:]: round(b, 4) for u, b in best.items()}, "| locked", {u[-6:]: round(R[u].get(lock, float("nan")), 4) for u in R})
key = lambda i: min(R[u][i] - best[u] for u in R)
for i in sorted(ids, key=lambda i: -key(i))[:8]:
    print(f"{i:32s} worst gap {key(i):+.4f} plateau {sum(best[u] - R[u][i] <= eps for u in R)}/{len(R)} |", " ".join(f"{u[-5:]} {R[u][i]:.4f}" for u in R))
