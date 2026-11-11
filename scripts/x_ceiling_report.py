"""Exploratory: X_ceiling report from results/p4/x_ceiling/<family>_<crop>.json (scripts/x_ceiling.py).
Per crop: a labelling exceeds its ceiling when the mean over draws of (its half-B score - the oracle at its call rates)
minus 2 standard errors is > 0. Per pair: crops -> units (registration/units/units.tsv, mean over a unit's crops);
means over units of floor, learned and ceiling (oracle at the learned model's call rates) on half B;
headroom share = (learned - floor) / (ceiling at the floor's call rates - floor) on those means. Full-data scores are checked against the
registered per-crop scores (<round-1 runs>/f5_<crop>_<model>.json, locked configuration).
Usage (env floors): python scripts/x_ceiling_report.py <round-1 runs dir> -> results/reports/X_ceiling.{json,md}"""

import csv, glob, json, sys
from pathlib import Path
import numpy as np

P = Path(__file__).resolve().parents[1]
R = Path(sys.argv[1])
UNIT = {r["crop"]: r["unit"] for r in csv.DictReader(open(P / "registration/units/units.tsv"), delimiter="\t")}
LOCK = json.loads((P / "manifests/locked_v3.json").read_text())["locked"]["F5"]
STEM = dict(floor="floors", scanvi="scanvi", tangram="tangram")
out, md = dict(F5=dict(crops={}, pairs={})), []
crops = {}
for f in sorted(glob.glob(str(P / "results/p4/x_ceiling/f5_*.json"))):
    r = json.loads(Path(f).read_text())
    u, dr = r["unit"], r["half_b_draws"]
    exceed = {}
    for n in r["full"]:
        d = np.array(dr[n]) - np.array(dr[f"oracle_{n}"])
        exceed[n] = bool(d.mean() - 2 * d.std(ddof=1) / np.sqrt(len(d)) > 0)
    reg = {
        n: json.loads((R / f"f5_{u}_{STEM[n]}.json").read_text())["grid"][LOCK[STEM[n]]["human"]["config"]]["heldout_auroc"]
        for n in r["full"]
    }
    crops[u] = dict(half_b=r["half_b_mean"], exceeds_ceiling=exceed, full_equals_registered=all(r["full"][n] == reg[n] for n in reg))
out["F5"]["crops"] = crops
units = sorted({UNIT[c] for c in crops})
um = lambda key: np.array([np.mean([crops[c]["half_b"][key] for c in crops if UNIT[c] == x]) for x in units])
for m in ("scanvi", "tangram"):
    fl, le, ce, cf = um("floor"), um(m), um(f"oracle_{m}"), um("oracle_floor")
    out["F5"]["pairs"][f"F5|{m}"] = dict(
        n_units=len(units),
        floor_half_b=float(fl.mean()),
        learned_half_b=float(le.mean()),
        ceiling_half_b=float(ce.mean()),
        ceiling_at_floor_rate_half_b=float(cf.mean()),
        headroom_share=float((le.mean() - fl.mean()) / (cf.mean() - fl.mean())),
        crops_learned_above_ceiling=int(sum(v["exceeds_ceiling"][m] for v in crops.values())),
        crops_floor_above_ceiling=int(sum(v["exceeds_ceiling"]["floor"] for v in crops.values())),
        n_crops=len(crops),
    )
out["F5"]["full_scores_equal_registered"] = all(v["full_equals_registered"] for v in crops.values())
# F2: one crop per unit; floor checked against the registered per-crop floor (P5_corrected f2_blob_reference)
P5 = json.loads((P / "results/reports/P5_corrected.json").read_text())["f2_blob_reference"]
c2 = {r["unit"]: r for r in (json.loads(Path(f).read_text()) for f in sorted(glob.glob(str(P / "results/p4/x_ceiling/f2_*.json"))))}
out["F2"] = dict(crops={u: {k: r[k] for k in ("ceiling", "ceiling_quadrant_folds", "scores")} for u, r in c2.items()}, pairs={})
out["F2"]["floor_equals_registered"] = all(r["scores"]["floor"] == P5[u]["floor"] for u, r in c2.items())
for m in ("cellcharter", "novae"):
    cs = [r for r in c2.values() if m in r["scores"]]
    fl, le, ce = (np.array([r["scores"]["floor"] for r in cs]), np.array([r["scores"][m] for r in cs]), np.array([r["ceiling"] for r in cs]))
    out["F2"]["pairs"][f"F2|{m}"] = dict(
        n_units=len(cs),
        floor=float(fl.mean()),
        learned=float(le.mean()),
        ceiling=float(ce.mean()),
        ceiling_quadrant_folds=float(np.mean([r["ceiling_quadrant_folds"] for r in cs])),
        headroom_share=float((le.mean() - fl.mean()) / (ce.mean() - fl.mean())),
        units_learned_above_ceiling=int((le > ce).sum()),
        units_floor_above_ceiling=int((fl > ce).sum()),
    )
(P / "results/reports/X_ceiling.json").write_text(json.dumps(out, indent=1, default=float))
md = [
    "# X_ceiling (exploratory)",
    "",
    "F5: held-out-marker counts split by binomial thinning into halves A and B (5 draws x 2 directions). Ceiling = per",
    "class, the cells with the highest half-A marker score called at the learned model's call rates, scored on half B.",
    "Floor and learned labellings scored on half B. Means over units; headroom share = (learned - floor) / (ceiling at",
    "the floor's call rates - floor). Above ceiling: mean over draws of (score - own-rate oracle) - 2 SE > 0.",
    "",
    "| pair | units | floor | learned | ceiling | ceiling at floor rate | headroom share | crops learned > ceiling | crops floor > ceiling |",
    "|---|---|---|---|---|---|---|---|---|",
]
for k, v in out["F5"]["pairs"].items():
    md.append(
        f"| {k} | {v['n_units']} | {v['floor_half_b']:.4f} | {v['learned_half_b']:.4f} | {v['ceiling_half_b']:.4f} | "
        f"{v['ceiling_at_floor_rate_half_b']:.4f} | {v['headroom_share']:+.3f} | {v['crops_learned_above_ceiling']}/{v['n_crops']} | "
        f"{v['crops_floor_above_ceiling']}/{v['n_crops']} |"
    )
md += ["", f"Full-data scores equal the registered per-crop scores: {out['F5']['full_scores_equal_registered']}", ""]
md += [
    "F2: ceiling = out-of-fold logistic regression from the locked floor's features (ell 400 µm, λ .8) to the annotated",
    "regions, 5 random cell folds, balanced ARI on the scored cells; quadrant folds (4 spatial blocks) reported beside it.",
    "Floor and learned: balanced ARI of the saved partitions (seed mean). One crop per unit; means over units.",
    "",
    "| pair | units | floor | learned | ceiling | ceiling, quadrant folds | headroom share | units learned > ceiling | units floor > ceiling |",
    "|---|---|---|---|---|---|---|---|---|",
]
for k, v in out["F2"]["pairs"].items():
    md.append(
        f"| {k} | {v['n_units']} | {v['floor']:.4f} | {v['learned']:.4f} | {v['ceiling']:.4f} | {v['ceiling_quadrant_folds']:.4f} | "
        f"{v['headroom_share']:+.3f} | {v['units_learned_above_ceiling']}/{v['n_units']} | {v['units_floor_above_ceiling']}/{v['n_units']} |"
    )
md += ["", f"Floor equals the registered per-crop floor: {out['F2']['floor_equals_registered']}", ""]
# F3: crops -> units (mean over a unit's crops); scores = the registered per-crop unit scores
c3 = {r["unit"]: r for r in (json.loads(Path(f).read_text()) for f in sorted(glob.glob(str(P / "results/p4/x_ceiling/f3_*.json"))))}
if c3:
    out["F3"] = dict(crops={u: {k: r[k] for k in ("ceiling", "scores")} for u, r in c3.items()}, pairs={})
    for m in ("lightgbm", "ncem"):
        cs = [c for c, r in c3.items() if m in r["scores"] and "floor" in r["scores"] and r["ceiling"] is not None]
        us = sorted({UNIT[c] for c in cs})
        mu = lambda f: np.array([np.mean([f(c3[c]) for c in cs if UNIT[c] == x]) for x in us])
        fl, le, ce = mu(lambda r: r["scores"]["floor"]), mu(lambda r: r["scores"][m]), mu(lambda r: r["ceiling"])
        out["F3"]["pairs"][f"F3|{m}"] = dict(
            n_units=len(us), floor=float(fl.mean()), learned=float(le.mean()), ceiling=float(ce.mean()),
            headroom_share=float((le.mean() - fl.mean()) / (ce.mean() - fl.mean())),
            units_learned_above_ceiling=int((le > ce).sum()), units_floor_above_ceiling=int((fl > ce).sum()),
        )
    md += [
        "F3: ceiling = per receiver, masked-gene counts split by binomial thinning (2 draws x 2 directions); per held-out",
        "quadrant the full-data training fit (depth covariate, z-score, PCA top 10) applied to each half's test rows;",
        "per PC sqrt of the Spearman-Brown corrected half-A / half-B Spearman rho; mean over PCs, folds, draws, receivers.",
        "Floor and learned: registered per-crop unit scores. Crops -> units as registered; means over units.",
        "",
        "| pair | units | floor | learned | ceiling | headroom share | units learned > ceiling | units floor > ceiling |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for k, v in out["F3"]["pairs"].items():
        md.append(
            f"| {k} | {v['n_units']} | {v['floor']:.4f} | {v['learned']:.4f} | {v['ceiling']:.4f} | {v['headroom_share']:+.4f} | "
            f"{v['units_learned_above_ceiling']}/{v['n_units']} | {v['units_floor_above_ceiling']}/{v['n_units']} |"
        )
    md.append("")
# F1: oracle reference (not a bound)
c1 = {r["unit"]: r for r in (json.loads(Path(f).read_text()) for f in sorted(glob.glob(str(P / "results/p4/x_ceiling/f1_*.json"))))}
if c1:
    out["F1"] = dict(crops={u: {k: r[k] for k in ("ceiling", "scores")} for u, r in c1.items()}, pairs={})
    for m in ("scvi", "nicheformer"):
        cs = [c for c, r in c1.items() if m in r["scores"] and "floor" in r["scores"]]
        us = sorted({UNIT[c] for c in cs})
        mu = lambda f: np.array([np.mean([f(c1[c]) for c in cs if UNIT[c] == x]) for x in us])
        fl, le, ce = mu(lambda r: r["scores"]["floor"]), mu(lambda r: r["scores"][m]), mu(lambda r: r["ceiling"])
        out["F1"]["pairs"][f"F1|{m}"] = dict(
            n_units=len(us), floor=float(fl.mean()), learned=float(le.mean()), reference=float(ce.mean()),
            share_of_reference_gap=float((le.mean() - fl.mean()) / (ce.mean() - fl.mean())) if ce.mean() - fl.mean() >= 0.01 else None,
            units_learned_above_reference=int((le > ce).sum()), units_floor_above_reference=int((fl > ce).sum()),
        )
    md += [
        "F1: oracle reference (not a bound) = neighbours from the expression PCA (50 PCs) of the held-out count half",
        "itself, kNN-Poisson score as registered (count split seed 0). Floor and learned: registered rescored scores.",
        "Share of reference gap = (learned - floor) / (reference - floor); n/a when reference - floor < .01.",
        "",
        "| pair | units | floor | learned | reference | share of reference gap | units learned > reference | units floor > reference |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for k, v in out["F1"]["pairs"].items():
        sh = "n/a" if v["share_of_reference_gap"] is None else format(v["share_of_reference_gap"], "+.3f")
        md.append(
            f"| {k} | {v['n_units']} | {v['floor']:.4f} | {v['learned']:.4f} | {v['reference']:.4f} | {sh} | "
            f"{v['units_learned_above_reference']}/{v['n_units']} | {v['units_floor_above_reference']}/{v['n_units']} |"
        )
    md.append("")
CK = P / "results/p4/x_ceiling/checks"
ck = {f.stem: json.loads(f.read_text()) for f in sorted(CK.glob("*.json"))}
out["checks"] = ck
md += ["Checks before the round-1 units (results/p4/x_ceiling/checks): F5 on development unit gse1_232 (floor = v2",
       "configuration), F2 on sim.domain_toy and on gse1_232, F3 on gse1_232 (development grid has no v3 floor; LightGBM",
       "shown).", ""]
for n, r in ck.items():
    if n.startswith("f5"):
        h = r["half_b_mean"]
        md.append(f"- {n}: half B " + ", ".join(f"{k} {v:.4f}" for k, v in h.items()))
    else:
        q = f", quadrant folds {r['ceiling_quadrant_folds']:.4f}" if "ceiling_quadrant_folds" in r else ""
        md.append(f"- {n}: ceiling {r['ceiling']:.4f}{q}, " + ", ".join(f"{k} {v:.4f}" for k, v in r["scores"].items()))
md.append("")
(P / "results/reports/X_ceiling.json").write_text(json.dumps(out, indent=1, default=float))
(P / "results/reports/X_ceiling.md").write_text("\n".join(md))
print("\n".join(md))
