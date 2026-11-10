"""Exploratory: configurations explored per model, floors vs learned models.
(a) development landscape (runs/landscape/f*_*.jsonl, perturb == none): distinct configurations per family for the
    floor side (every non-learned part: grid, ref, ist, logreg, anchor, ...) and per learned model, with the units run;
    a configuration = the row id without the split seed and perturbation. GAT (F3) was dropped before the lock.
(b) tier-2 grid per model (manifests/tier2.json "grids"); (c) floor lock files (locked -> v2 -> v3).
Counts only; nothing is re-scored. The landscape covers the development and validation units; later floor revisions
 ran under runs/g7v3 and are counted as lock versions, not configurations.
Usage: python scripts/x_effort.py -> runs/p4/X_effort.{json,md}"""

import json, glob, re
from collections import defaultdict
from pathlib import Path

cfg, units = defaultdict(set), defaultdict(set)
for f in glob.glob("runs/landscape/f*_*.jsonl"):
    fam = Path(f).name.split("_")[0].upper()
    for line in open(f):
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            continue
        if r.get("perturb", "none") != "none":
            continue
        side = r.get("model") if r.get("part") == "learned" else f"floor: {r.get('part', '?')}"
        c = re.sub(r"\|split\d+$", "", r["id"]).replace("|none", "")
        cfg[(fam, side)].add(c)
        units[(fam, side)].add(r.get("unit"))
G = json.loads(Path("manifests/tier2.json").read_text())["grids"]
locks = [
    p for p in ("manifests/locked.json", "manifests/locked_v2.json", "manifests/locked_v3.json") if Path(p).exists()
]
rows = [
    dict(
        family=f,
        side=s,
        dev_configs=len(cfg[(f, s)]),
        dev_units=sorted(units[(f, s)]),
        tier2_configs=len(G.get(f, {}).get("floor" if s.startswith("floor") else s, {}).get("configs", [])) or None,
    )
    for f, s in sorted(cfg, key=lambda k: (k[0], not k[1].startswith("floor"), k[1]))
]
Path("runs/p4/X_effort.json").write_text(json.dumps(dict(rows=rows, floor_lock_files=locks), indent=1))
md = [
    "# X_effort: configurations explored per model (exploratory; counts only)",
    "",
    f"Floor lock versions: {len(locks)} ({', '.join(Path(p).name for p in locks)}). Development configurations from runs/landscape (development + validation units, no perturbation).",
    "",
    "| Family | Model | Development configurations | Development units | Tier-2 grid |",
    "|---|---|---|---|---|",
]
md += [
    f"| {r['family']} | {r['side']} | {r['dev_configs']} | {len(r['dev_units'])} | {r['tier2_configs'] or '–'} |"
    for r in rows
]
Path("runs/p4/X_effort.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
