"""Exploratory: specification curve of the round-1 comparison. For each registered pair, every estimate of
D = learned - floor that was already computed under another analysis choice, collected from the saved reports
(no model is run, no cell is re-scored). Choices: units (leave one cancer out, one platform, subsets), aggregation
(slide, study), tuning (both sides at their best confirmation-unit configuration), metric (distance, call rate, depth,
scoring region, protein-gated F1), floor (earlier or stronger floor versions), model (other configurations of the
learned model) and scale (whole slides). Each row keeps the estimate, the one-sided 95% bounds and the pair's locked
margin. Bounds that a report gives only as an upper bound are completed by symmetry (lower = 2 est - upper); platform
subsets and the interior re-score are tested here from the saved per-unit D with the registered test. Two routes to
the same set of units count as one specification.
The rows of one pair share units and are not independent: the curve is descriptive, no test is made across rows.
Usage (env floors): python scripts/x_speccurve.py results/reports registration/units/units.tsv
  -> results/reports/X_speccurve.{json,md}"""

import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from spatialfloors import stats
from spatialfloors.provenance import commit

R = Path(sys.argv[1])
J = lambda f: json.loads((R / f).read_text())
U = pd.read_csv(sys.argv[2], sep="\t")
T1 = {k: v for k, v in J("P5_corrected.json")["pairs"].items() if "|" in k and "per_unit" in v}
as_crop = lambda u: "gse303070_" if u.startswith("Donor") else u  # GSE303070 units are donors
strata_of = lambda units: [stats.stratum(as_crop(u)) for u in units]
rows = []


def add(pair, choice, label, source, est, lower=None, upper=None, n=None, same_scale=True, bounds="report"):
    """One specification. same_scale = the score is on the registered metric's scale, so the margin applies."""
    if lower is None and upper is not None:
        lower, bounds = 2 * est - upper, "symmetry"
    mg = T1[pair]["margin"]
    rows.append(
        dict(
            pair=pair,
            choice=choice,
            label=label,
            source=source,
            n=n,
            est=est,
            lower=lower,
            upper=upper,
            margin=mg if same_scale else None,
            est_over_margin=est / mg if same_scale else None,
            noninferior=bool(upper < mg) if same_scale and upper is not None else None,
            floor_ahead=bool(upper < 0) if upper is not None else None,
            learned_ahead=bool(lower > 0) if lower is not None else None,
            bounds=bounds,
        )
    )


def test(units, D, mg):
    """The registered test on a set of units (platform strata centred when they leave a residual df)."""
    S = strata_of(units)
    t = stats.family_test(np.asarray(D, float), S if 1 < len(set(S)) < len(D) - 1 else None, margin=mg)
    return dict(est=t["est"], lower=t["lower"], upper=t["upper"], n=len(D), bounds="computed here")


# registered estimate, leave one cancer out (SPATCH counted as CRC), one platform at a time
LOCO = J("P5_loco_crc.json")["pairs"]
for k, v in T1.items():
    add(k, "registered", "registered", "P5_corrected", v["est"], v["lower"], v["upper"], v["n_units"])
    for c, t in (LOCO.get(k, {}).get("corrected") or {}).items():
        add(k, "units", f"without {c}", "P5_loco_crc", t["est"], None, t["upper"], t["n_units"])
    if "without_pretraining_overlap" in v:
        t = v["without_pretraining_overlap"]
        add(k, "units", "without units in the model's pretraining data", "P5_corrected", t["est"], None, t["upper"])
    units = sorted(v["per_unit"])
    S = strata_of(units)
    if len(set(S)) > 1:
        for s in sorted(set(S)):
            sub = [u for u, x in zip(units, S) if x == s]
            if len(sub) > 1:
                t = test(sub, [v["per_unit"][u] for u in sub], v["margin"])
                assert abs(t["est"] - v["by_stratum"][s]) < 1e-9, (k, s)  # the reported platform mean
                add(k, "units", f"{s} only", "P5_corrected", **t)

# slide and study level
for k, v in J("X_levels.json")["pairs"].items():
    for lev in ("slide", "study"):
        t = v[lev]
        if "est" in t and t["n"] != v["donor"]["n"]:
            add(k, "aggregation", f"one value per {lev}", "X_levels", t["est"], t["lower"], t["upper"], t["n"])

# both sides at their best configuration on the confirmation units
for k, t in J("T2_corrected.json")["capacity"].items():
    add(k, "tuning", "both sides at their best configuration (tier 2)", "T2_corrected", t["est"], t["lower"], t["upper"], t["n_units"])  # fmt: skip

# F1 with Euclidean neighbours
E = J("X_f1_euclid.json")["tests"]
for m in ("scvi", "nicheformer"):
    for lab, key in (
        ("Euclidean neighbours on both sides", f"{m} euclid − floor euclid"),
        ("Euclidean for the model, cosine for the floor", f"{m} euclid − floor cosine"),
        ("the model's better distance, cosine for the floor", f"{m} best − floor cosine"),
    ):
        t = E[key]
        add(f"F1|{m}", "metric", lab, "X_f1_euclid", t["est"], t["lower"], t["upper"], t["n_units"])

# F5 at the comparator's call rate, and with depth held fixed
M = J("X_f5_matched.json")["tests"]
DP = J("X_f5_depth.json")["pairs"]
for m in ("scanvi", "tangram"):
    for lab, key in (
        ("floor thresholded at its own call rate", "auroc_floor_ownrate"),
        ("floor at the model's call rate", "auroc_floor_matched"),
    ):
        t = M[f"{m} - {key} [mean]"]
        add(f"F5|{m}", "metric", lab, "X_f5_matched", t["est"], t["lower"], t["upper"], t["n_units"])
    t = DP[m]["depth_strat"]
    assert abs(DP[m]["heldout_auroc"]["est"] - T1[f"F5|{m}"]["est"]) < 1e-9
    add(f"F5|{m}", "metric", "AUROC within depth quintiles", "X_f5_depth", t["est"], t["lower"], t["upper"], t["n_units"])  # fmt: skip

# scoring restricted to the interior of the crop
for k, v in J("X_interior.json")["pairs"].items():
    r0 = v["by_radius"]["0"]
    t0 = test(r0["units"], r0["D"], T1[k]["margin"])
    assert max(abs(t0[x] - T1[k][x]) for x in ("est", "lower", "upper")) < 1e-9, k  # R = 0 is the registered test
    for r in ("100", "250", "500"):
        b = v["by_radius"][r]
        add(
            k,
            "metric",
            f"scored {r} µm inside the crop edge",
            "X_interior",
            **test(b["units"], b["D"], T1[k]["margin"]),
        )

# other floors, other configurations of the learned model, a protein-gated metric
S2 = J("P5_secondary.json")["comparisons"]
SEC = (
    ("F2|cellcharter", "floor", "multi-scale floor", "F2 CellCharter - multi-scale floor", True),
    ("F2|novae", "floor", "multi-scale floor", "F2 Novae - multi-scale floor (mouse)", True),
    ("F2|cellcharter", "model", "CellCharter on PCA instead of scVI", "F2 CellCharter on PCA - locked floor", True),
    ("F3|lightgbm", "floor", "floor without ring features (v2)", "F3 LightGBM - v2 floor", True),
    ("F3|ncem", "floor", "floor without ring features (v2)", "F3 NCEM - v2 floor", True),
    ("F3|ncem", "model", "encoder-decoder NCEM", "F3 ED-NCEM - floor", True),
    ("F5|scanvi", "floor", "floor without the spatial prior", "F5 scANVI - InSituType without prior (v2 floor)", True),
    (
        "F5|tangram",
        "floor",
        "floor without the spatial prior",
        "F5 Tangram - InSituType without prior (v2 floor)",
        True,
    ),
    ("F5|scanvi", "metric", "protein-gated macro-F1 (GSE1)", "F5 IF-gated macro-F1 (GSE1): scanvi - floor", False),
    ("F5|tangram", "metric", "protein-gated macro-F1 (GSE1)", "F5 IF-gated macro-F1 (GSE1): tangram - floor", False),
)
for pair, choice, lab, key, same in SEC:
    t = S2[key]
    add(pair, choice, lab, "P5_secondary", t["est"], t["lower"], t["upper"], t["n_units"], same_scale=same)
t = J("P5_corrected.json")["secondary"]["F3|ncem_vs_input_ridge"]
add("F3|ncem", "floor", "ridge on NCEM's own inputs", "P5_corrected", t["est"], t["lower"], t["upper"], t["n_units"])

# LightGBM with early stopping, Tangram on normalised counts
ES = J("X_f3_lgbm_es.json")["pairs"]
assert abs(ES["registered"]["est"] - T1["F3|lightgbm"]["est"]) < 1e-9
for lab, key in (
    ("early stopping, all features", "all_feats"),
    ("early stopping, the floor's features", "floor_feats"),
):
    t = ES[key]
    add("F3|lightgbm", "model", lab, "X_f3_lgbm_es", t["est"], t["lower"], t["upper"], t["n_units"])
TN = J("X_f5_tangram_norm.json")["pairs"]
assert abs(TN["raw"]["est"] - T1["F5|tangram"]["est"]) < 1e-9
for lab, key in (("input as counts per 10k", "cp10k"), ("input as log counts per 10k", "logcp10k")):
    t = TN[key]
    add("F5|tangram", "model", lab, "X_f5_tangram_norm", t["est"], t["lower"], t["upper"], t["n_units"])

# F4 by cancer and by provenance
for k, v in J("X_f4_split.json")["pairs"].items():
    for lab, t in v["subsets"].items():
        add(k, "units", f"{lab} only", "X_f4_split", t["mean"], t["lower"], t["upper"], t["n"])

# whole slides (six GSE1 slides)
W = J("X_wholeslide.json")["pairs"]
for k, lab in (
    ("F1|scvi", "whole slides (6 GSE1)"),
    ("F2|cellcharter", "whole slides (6 GSE1)"),
    ("F3|lightgbm", "whole slides (6 GSE1)"),
    ("F5|scanvi", "whole slides (6 GSE1)"),
):
    t = W[k]
    add(k, "scale", lab, "X_wholeslide", t["est"], t["lower"], t["upper"], t["n"])
for k, key, lab in (
    ("F3|lightgbm", "F3 LightGBM - v2 floor", "whole slides, floor without ring features (v2)"),
    ("F5|scanvi", "F5 scANVI - floor at scANVI's call rate", "whole slides, floor at the model's call rate"),
):
    t = W[key]
    add(k, "scale", lab, "X_wholeslide", t["est"], t["lower"], t["upper"], t["n"])

# the same set of units reached by two routes (e.g. "without CRC" = "MERFISH mouse only") is one specification
merged = []
for x in rows:
    same = [
        y
        for y in merged
        if y["pair"] == x["pair"]
        and abs(y["est"] - x["est"]) < 1e-9
        and None not in (x["upper"], y["upper"])
        and abs(y["upper"] - x["upper"]) < 1e-9
    ]
    if same and x["choice"] == same[0]["choice"] == "units":
        same[0]["label"] += " = " + x["label"]
        same[0]["source"] += ", " + x["source"] if x["source"] not in same[0]["source"] else ""
    else:
        merged.append(x)
rows = merged

# summary per pair (rows on the registered scale)
summary = {}
for k in T1:
    r = [x for x in rows if x["pair"] == k and x["margin"] is not None]
    e = np.array([x["est"] for x in r])
    summary[k] = dict(
        n_specs=len(r),
        est_min=float(e.min()),
        est_median=float(np.median(e)),
        est_max=float(e.max()),
        registered=T1[k]["est"],
        margin=T1[k]["margin"],
        n_est_negative=int((e < 0).sum()),
        n_noninferior=sum(x["noninferior"] is True for x in r),
        n_floor_ahead=sum(x["floor_ahead"] is True for x in r),
        n_learned_ahead=sum(x["learned_ahead"] is True for x in r),
        n_not_noninferior=sum(x["noninferior"] is False for x in r),
        other_scale=len([x for x in rows if x["pair"] == k and x["margin"] is None]),
    )
tot = [x for x in rows if x["margin"] is not None]
overall = dict(
    n_specs=len(tot),
    n_est_negative=sum(x["est"] < 0 for x in tot),
    n_noninferior=sum(x["noninferior"] is True for x in tot),
    n_floor_ahead=sum(x["floor_ahead"] is True for x in tot),
    n_learned_ahead=sum(x["learned_ahead"] is True for x in tot),
    n_other_scale=len(rows) - len(tot),
)
out = dict(commit=commit(), overall=overall, summary=summary, specifications=rows)

f = lambda x: "–" if x is None else f"{x:+.4f}"
yn = lambda x: "–" if x is None else ("yes" if x else "no")
md = [
    "# X_speccurve: specification curve of the round-1 comparison (exploratory)",
    "",
    "D = learned − floor under every choice already computed in the saved reports; one-sided 95% bounds; δ = the "
    "pair's locked margin. Rows of one pair share units and are not independent: descriptive only.",
    "",
    f"On the registered scale: {overall['n_specs']} specifications; D < 0 in {overall['n_est_negative']}; floor "
    f"non-inferior in {overall['n_noninferior']}; floor ahead (upper bound < 0) in {overall['n_floor_ahead']}; learned "
    f"model ahead (lower bound > 0) in {overall['n_learned_ahead']}. {overall['n_other_scale']} further specifications "
    "use another score scale (no margin).",
    "",
    "| Pair | Specifications | Registered D | Min | Median | Max | δ | D < 0 | Non-inferior | Floor ahead | Learned ahead |",
    "|---|---|---|---|---|---|---|---|---|---|---|",
]
for k, s in summary.items():
    md.append(
        f"| {k} | {s['n_specs']} | {s['registered']:+.4f} | {s['est_min']:+.4f} | {s['est_median']:+.4f} | "
        f"{s['est_max']:+.4f} | {s['margin']:.4f} | {s['n_est_negative']} | {s['n_noninferior']} | "
        f"{s['n_floor_ahead']} | {s['n_learned_ahead']} |"
    )
md += [
    "",
    "Bounds: report = as in the source report; symmetry = lower bound completed as 2 est − upper; computed here = the "
    "registered test on the saved per-unit D.",
]
for k in T1:
    md += [
        "",
        f"## {k}",
        "",
        "| Choice | Specification | Source | n | D | Bounds | D / δ | Non-inferior | Bounds from |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for x in sorted([x for x in rows if x["pair"] == k], key=lambda x: (x["margin"] is None, x["est"])):
        md.append(
            f"| {x['choice']} | {x['label']} | {x['source']} | {x['n'] if x['n'] is not None else '–'} | {f(x['est'])} | "
            f"{f(x['lower'])}, {f(x['upper'])} | "
            f"{'other scale' if x['margin'] is None else format(x['est_over_margin'], '+.2f')} | "
            f"{yn(x['noninferior'])} | {x['bounds']} |"
        )
(R / "X_speccurve.json").write_text(json.dumps(out, indent=1, default=float))
(R / "X_speccurve.md").write_text("\n".join(md) + "\n")
print("\n".join(md[:22]))
print("X_SPECCURVE_DONE", len(rows))
