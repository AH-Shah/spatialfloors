"""Leave-one-cancer-out (PROTOCOL §4, sensitivity) with the SPATCH Xenium 5K sample in the CRC group.
The registered analysis code (p5_analyze.py) labels SPATCH "tissue unverified", so its leave-one-cancer-out leaves
SPATCH out as a group of its own. The sample is colon adenocarcinoma (DATA.md); here it counts as CRC. Recomputed from
the per-unit D of the registered tier-1 report with the registered test (stats.family_test, platform strata, margin).
Check: with the registered labels every registered leave-one-cancer-out value is reproduced (|diff| < 1e-12).
Usage (env floors): python scripts/p5_loco_crc.py results/reports/P5_corrected.json -> P5_loco_crc.{json,md} beside it"""

import hashlib, json, sys
from pathlib import Path
import numpy as np
from spatialfloors import stats
from spatialfloors.provenance import commit

SRC = Path(sys.argv[1])
P = json.loads(SRC.read_text())["pairs"]
REGISTERED = {  # p5_analyze.py CANCER, as registered
    "gse1": "CRC",
    "gse0": "CRC",
    "gse303070": "CRC",
    "spatch": "SPATCH (tissue unverified)",
    "allen": "mouse brain",
    "bruker_breast": "breast",
    "bruker_mobreast": "breast",
    "bruker_skin": "skin SCC",
    "bruker_kidney": "RCC",
}
CORRECTED = REGISTERED | {"spatch": "CRC"}  # colon adenocarcinoma (DATA.md)
as_crop = lambda u: "gse303070_" if u.startswith("Donor") else u  # GSE303070 units are donors
label = lambda m, u: next(v for k, v in sorted(m.items(), key=lambda kv: -len(kv[0])) if as_crop(u).startswith(k))


def loco(v, labels):
    """As p5_analyze.py: drop one cancer at a time, test the remaining units (needs > 2 units and > 1 cancer)."""
    units = sorted(v["per_unit"])
    canc = {u: label(labels, u) for u in units}
    if len(set(canc.values())) < 2 or len(units) < 3:
        return None
    out = {}
    for k in sorted(set(canc.values())):
        keep = [u for u in units if canc[u] != k]
        if len(keep) > 1:
            D = np.array([v["per_unit"][u] for u in keep])
            t = stats.family_test(D, [stats.stratum(as_crop(u)) for u in keep], margin=v["margin"])
            out[k] = {kk: t[kk] for kk in ("est", "upper", "noninferior")} | {"n_units": len(keep)}
    return out


out = dict(
    source=SRC.name,
    source_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),
    commit=commit(),
    relabelled={u: "CRC" for v in P.values() for u in v.get("per_unit", {}) if u.startswith("spatch")},
    pairs={},
)
for k, v in P.items():
    if "|" not in k or "per_unit" not in v:
        continue
    reg, cor = loco(v, REGISTERED), loco(v, CORRECTED)
    ref = v.get("leave_one_cancer_out")
    assert (reg is None) == (ref is None), k
    for g in ref or {}:  # the registered labels reproduce the registered report
        assert set(reg) == set(ref), (k, sorted(reg), sorted(ref))
        for f in ("est", "upper"):
            assert abs(reg[g][f] - ref[g][f]) < 1e-12, (k, g, f, reg[g][f], ref[g][f])
        assert reg[g]["noninferior"] == ref[g]["noninferior"], (k, g)
    if cor is None:
        continue
    out["pairs"][k] = dict(
        changed=reg != cor,
        registered=reg,
        corrected=cor,
        noninferior_every_split=dict(
            registered=all(s["noninferior"] for s in reg.values()),
            corrected=all(s["noninferior"] for s in cor.values()),
        ),
    )
OUT = SRC.parent / "P5_loco_crc"
OUT.with_suffix(".json").write_text(json.dumps(out, indent=1, default=float))

md = [
    "# Leave-one-cancer-out with SPATCH Xenium 5K counted as CRC",
    "",
    f"Source: `{SRC.name}` (sha256 {out['source_sha256'][:12]}), per-unit D of the registered tier-1 report. "
    "Registered labels reproduce every registered value; relabelled: "
    + ", ".join(f"`{u}`" for u in out["relabelled"])
    + " (SPATCH, tissue unverified -> CRC).",
    "",
    "| Pair | Left out | Registered: units, est, upper, NI | Corrected: units, est, upper, NI |",
    "|---|---|---|---|",
]
fmt = lambda s: (
    f"{s['n_units']}, {s['est']:+.4f}, {s['upper']:+.4f}, {'yes' if s['noninferior'] else 'no'}" if s else "-"
)
yn = lambda b: "yes" if b else "no"
nm = lambda k: k.replace("|", " ")  # pair keys hold "|", the markdown cell separator
for k, r in out["pairs"].items():
    if not r["changed"]:
        continue
    for g in sorted(set(r["registered"]) | set(r["corrected"])):
        md.append(f"| {nm(k)} | {g} | {fmt(r['registered'].get(g))} | {fmt(r['corrected'].get(g))} |")
md += [
    "",
    "Pairs whose splits do not change (SPATCH not among their units): "
    + ", ".join(nm(k) for k, r in out["pairs"].items() if not r["changed"])
    + ".",
    "",
    "| Pair | Non-inferior in every split: registered | corrected |",
    "|---|---|---|",
]
md += [
    f"| {nm(k)} | {yn(r['noninferior_every_split']['registered'])} | {yn(r['noninferior_every_split']['corrected'])} |"
    for k, r in out["pairs"].items()
]
OUT.with_suffix(".md").write_text("\n".join(md) + "\n")
print("\n".join(md))
