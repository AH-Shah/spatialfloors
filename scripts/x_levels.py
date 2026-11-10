"""Exploratory: tier-1 D re-aggregated at three replication
levels. unit = registered unit (patient, donor or animal); slide = units on one physical slide averaged (GSE303070's 11
donors sit on 5 slides; elsewhere unit = slide); study = units of one data release averaged. Same test as tier 1
(stats.family_test, centred platform strata when the strata leave df >= 1, else none). The donor level must reproduce
the registered estimate (checked).
Usage (env floors): python scripts/x_levels.py runs/p4/P5_corrected.json registration/units/units.tsv ->
runs/p4/X_levels.{json,md}"""

import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from spatialfloors import stats
from spatialfloors.provenance import commit

T1 = json.loads(Path(sys.argv[1]).read_text())["pairs"]
OUT = Path(sys.argv[1]).parent
U = pd.read_csv(sys.argv[2], sep="\t")
STUDY = (
    ("gse1_", "GSE1 CRC (Crowell et al.)"),
    ("Donor", "GSE303070"),
    ("spatch_", "SPATCH"),
    ("allen_yao", "Allen Yao"),
    ("allen_animal", "Allen Zhuang"),
    ("bruker_kidney", "Bruker kidney"),
    ("bruker_skin", "Bruker skin"),
    ("bruker_mobreast", "Bruker breast multiomic"),
    ("bruker_breast", "Bruker breast"),
)
study = lambda u: next(s for p, s in STUDY if u.startswith(p))
slide_of, stratum_of = {}, {}
for _, r in U.iterrows():
    stratum_of[r.unit] = r.stratum
    slide_of[r.unit] = (
        "gse303070_" + "_".join(r.crop.split("_")[1:-1]) if str(r.crop).startswith("gse303070_") else r.unit
    )


def test(groups, mg):
    keys = sorted(groups)
    D = np.array([np.mean([g[0] for g in groups[k]]) for k in keys])
    S = [groups[k][0][1] for k in keys]
    use = S if 1 < len(set(S)) < len(D) - 1 else None
    t = stats.family_test(D, use, margin=mg)
    return dict(
        n=len(keys),
        est=t["est"],
        lower=t["lower"],
        upper=t["upper"],
        noninferior=t["noninferior"],
        strata=use is not None,
        groups=keys,
    )


out = dict(commit=commit(), pairs={})
md = [
    "# X_levels: tier-1 D at donor, slide and study level (exploratory)",
    "",
    "| Pair | Level | n | D | One-sided 95% bounds | Non-inferior |",
    "|---|---|---|---|---|---|",
]
for k, v in T1.items():
    if "|" not in k or "per_unit" not in v:
        continue
    mg = v["margin"]
    pu = v["per_unit"]
    o = {}
    for lev, key in (("donor", lambda u: u), ("slide", lambda u: slide_of.get(u, u)), ("study", study)):
        g = {}
        for u, d in pu.items():
            g.setdefault(key(u), []).append((d, stratum_of.get(u, "?")))
        if len(g) < 2:
            o[lev] = dict(n=len(g), note="one group")
            continue
        o[lev] = test(g, mg)
    assert abs(o["donor"]["est"] - v["est"]) < 1e-9, (k, o["donor"]["est"], v["est"])
    out["pairs"][k] = o
    for lev, r in o.items():
        md.append(
            f"| {k} | {lev} | {r['n']} | "
            + (
                f"{r['est']:+.4f} | {r['lower']:+.4f}, {r['upper']:+.4f} | {'yes' if r['noninferior'] else 'NO'} |"
                if "est" in r
                else "– | – | one group |"
            )
        )
(OUT / "X_levels.json").write_text(json.dumps(out, indent=1, default=float))
(OUT / "X_levels.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
