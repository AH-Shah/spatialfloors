"""Exploratory: F4 tier-1 D split by cancer (CRC = GSE1
slides, mouse = Allen sections) and by provenance (prospective = run here; retrospective = the Cellina authors'
published
runs, PROTOCOL §3). Reads the registered tier-1 report only; same test as tier 1 (one-sided 95% bounds, exact t,
margin .03), no strata (each subset is one platform). The subsets were chosen after the tier-1 report was read.
Usage (env floors): python scripts/x_f4_split.py runs/p4/P5_corrected.json -> runs/p4/X_f4_split.{json,md}"""

import json, sys
from pathlib import Path
import numpy as np
from spatialfloors import stats
from spatialfloors.provenance import commit

T1 = json.loads(Path(sys.argv[1]).read_text())["pairs"]
OUT = Path(sys.argv[1]).parent
PROSPECTIVE = {"gse1_110", "gse1_222", "allen_animal1", "allen_animal3", "allen_yao638850"}
SUBSETS = {
    "CRC (GSE1)": lambda u: u.startswith("gse1_"),
    "mouse (Allen)": lambda u: u.startswith("allen_"),
    "prospective": lambda u: u in PROSPECTIVE,
    "retrospective": lambda u: u not in PROSPECTIVE,
    "CRC prospective": lambda u: u.startswith("gse1_") and u in PROSPECTIVE,
    "CRC retrospective": lambda u: u.startswith("gse1_") and u not in PROSPECTIVE,
}
out = dict(commit=commit(), source=sys.argv[1], pairs={})
md = [
    "# X_f4_split: F4 by cancer and by provenance (exploratory; subsets chosen after tier 1)",
    "",
    "| Pair | Subset | n | Mean D | One-sided 95% bounds | Above .03 (units) | Ahead (units) |",
    "|---|---|---|---|---|---|---|",
]
for k in ("F4_pearson|cellina", "F4_pearson|cpa", "F4_pearson_clean|cellina", "F4_pearson_clean|cpa"):
    pu = T1[k]["per_unit"]
    out["pairs"][k] = dict(per_unit=pu, subsets={})
    for name, f in SUBSETS.items():
        us = sorted(u for u in pu if f(u))
        D = np.array([pu[u] for u in us])
        if len(us) < 2:
            continue
        t = stats.family_test(D, None, margin=0.03)
        o = dict(
            n=len(us),
            units=us,
            mean=float(D.mean()),
            lower=t["lower"],
            upper=t["upper"],
            p_superior_uncorrected=t["p_superior"],
            n_above_margin=int((D > 0.03).sum()),
            n_ahead=int((D > 0).sum()),
        )
        out["pairs"][k]["subsets"][name] = o
        md.append(
            f"| {k} | {name} | {o['n']} | {o['mean']:+.4f} | {o['lower']:+.4f}, {o['upper']:+.4f} | {o['n_above_margin']} | {o['n_ahead']} |"
        )
md += [
    "",
    "Per unit (D = learned − floor): "
    + "; ".join(
        f"{k}: " + ", ".join(f"{u} {v:+.4f}" for u, v in sorted(T1[k]["per_unit"].items())) for k in out["pairs"]
    ),
]
(OUT / "X_f4_split.json").write_text(json.dumps(out, indent=1, default=float))
(OUT / "X_f4_split.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
