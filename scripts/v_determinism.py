"""Validation gate G4 (PLAN §17): a floor re-run must reproduce its first run bit-identically on CPU. Compares the numeric
leaves of two result JSONs (provenance / timing fields ignored). Usage: python scripts/v_determinism.py <a.json> <b.json>"""
import json, sys
SKIP = {"s", "fit_s", "spatialfloors_commit", "time_s"}


def leaves(x, p=""):
    if isinstance(x, dict): return {k: v for kk, vv in x.items() if kk not in SKIP for k, v in leaves(vv, f"{p}/{kk}").items()}
    if isinstance(x, list): return {k: v for i, vv in enumerate(x) for k, v in leaves(vv, f"{p}/{i}").items()}
    return {p: x} if isinstance(x, (int, float)) and not isinstance(x, bool) else {}


a, b = (leaves(json.load(open(f))) for f in sys.argv[1:3]); diff = {k: (a[k], b.get(k)) for k in a if a[k] != b.get(k)}
print("G4", sys.argv[1], "identical" if not diff and a.keys() == b.keys() else f"DIFFERS {len(diff)} / {len(a)}: {list(diff.items())[:5]}"); sys.exit(1 if diff else 0)
