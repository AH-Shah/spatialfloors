"""NCEM test configuration (PROTOCOL §3): the lock rule (highest development mean, ties -> grid order) applied to the
authors' NCEM grid {ED-NCEM, linear interactions} x radius {19, 46 um} on the development units (runs/p3v3, rows 379).
manifests/locked_v3.json lists both models at the selected radius; this file names the one that forms the test pair.
Usage: python scripts/lock_ncem.py -> manifests/locked_v3_ncem.json"""
import json
from pathlib import Path
import numpy as np

DEV = {"human": ["gse1_232", "gse0_S0", "gse303070_148191A_118N", "gse303070_148191A_278N", "gse303070_148191A_278T", "gse303070_148191A_46T"],
       "mouse": ["allen_animal4_C57BL6J-4.003"]}
GRID = ["edA_e1_r19", "linA_l20_r19", "edA_e1_r46", "linA_l20_r46"]
out = dict(rule=__doc__.split("\n")[0], source="runs/p3v3/f3_<unit>_ncemA.json", locked={})
for s, us in DEV.items():
    g = {u: json.loads(Path(f"runs/p3v3/f3_{u}_ncemA.json").read_text())["grid"] for u in us}
    m = {c: float(np.mean([g[u][c]["unit"] for u in us])) for c in GRID}; best = max(GRID, key=lambda c: (m[c], -GRID.index(c)))
    out["locked"][s] = dict(config=best, dev_mean=m[best], grid_means=m, units=us)
Path("manifests/locked_v3_ncem.json").write_text(json.dumps(out, indent=1) + "\n")
print({s: (v["config"], round(v["dev_mean"], 4)) for s, v in out["locked"].items()})
