"""Replication round: predicted power per test pair from the round-1 per-unit differences (pilot data only).
Per pair: sd of D = learned - floor over the round-1 confirmation units (residual sd of the registered model, platform
strata centred), N = the planned round-2 units of the family (registration_v2/SELECTION.md), df = N - strata, margin as
locked in round 1. Floor non-inferior if the one-sided 95% upper bound of mean D < margin: power at the round-1 mean D,
at D = 0 and at D = +margin/2. Learned superior (Holm over the round-2 pairs, first step alpha / pairs): power at
D = +margin and the smallest D detected with 80% power. Sensitivity: the same with sd x 1.5 (new tissues, platforms).
Joint = product over pairs. N is the plan before the metric-blind preparation and is updated with units.tsv.
Usage (env floors): python scripts/v2_power.py results/reports/P5_corrected.json
  -> registration_v2/configs/power_v2.json, registration_v2/POWER.md"""

import json, sys
from pathlib import Path
import numpy as np, scipy.stats as st
from scipy.optimize import brentq
from spatialfloors.provenance import commit

PUB = Path(__file__).resolve().parents[1]
T1 = json.loads(Path(sys.argv[1]).read_text())["pairs"]
# planned round-2 units and platform strata (Xenium Prime 5K, CosMx 960-plex); F2 = the five 10x sections with H&E
N = dict(F1=(15, 2), F2=(5, 1), F3=(15, 2), F5=(13, 2))
PAIRS = ["F1|scvi", "F1|nicheformer", "F2|cellcharter", "F3|lightgbm", "F3|ncem", "F5|scanvi", "F5|tangram"]
ALPHA = 0.05
SD = (("power", 1.0), ("power_sd_x1.5", 1.5))


def power(mean, sd, m, n, strata):
    df, se = n - strata, sd / np.sqrt(n)
    tq, th = st.t.ppf(1 - ALPHA, df), st.t.ppf(1 - ALPHA / len(PAIRS), df)
    ni = lambda mu: float(st.nct.cdf(-tq, df, (mu - m) / se))
    sup = lambda mu: float(st.nct.sf(th, df, mu / se))
    return dict(
        df=df,
        se=float(se),
        noninferior=dict(round1_mean=ni(mean), zero=ni(0.0), half_margin=ni(m / 2)),
        superior_at_margin=sup(m),
        superior_d80=float(brentq(lambda mu: sup(mu) - 0.8, 0, 50 * se)),
    )


out = dict(commit=commit(), source=Path(sys.argv[1]).name, alpha=ALPHA, n_planned=N, pairs={})
joint = {lab: dict(round1_mean=1.0, zero=1.0, half_margin=1.0) for lab, _ in SD}
for k in PAIRS:
    v = T1[k]
    n, strata = N[k[:2]]
    sd = float(v["se"] * np.sqrt(v["n_units"]))  # residual sd of the registered model
    o = dict(
        round1_units=v["n_units"],
        round1_mean=v["est"],
        sd=sd,
        sd_raw=float(np.std(list(v["per_unit"].values()), ddof=1)),
        margin=v["margin"],
        N=n,
        strata=strata,
    )
    for lab, f in SD:
        o[lab] = power(v["est"], sd * f, v["margin"], n, strata)
        for x in joint[lab]:
            joint[lab][x] *= o[lab]["noninferior"][x]
    out["pairs"][k] = o
out["joint_noninferior"] = joint
cfg = PUB / "registration_v2/configs"
cfg.mkdir(parents=True, exist_ok=True)
(cfg / "power_v2.json").write_text(json.dumps(out, indent=1) + "\n")

two = lambda o, f: f"{f(o['power']):.2f} ({f(o['power_sd_x1.5']):.2f})"
md = [
    "# Replication round: predicted power (from round-1 pilot data)",
    "",
    "sd = residual sd of the per-unit D in round 1; N = planned round-2 units; δ = locked margin. Non-inferiority: "
    "one-sided 95% upper bound of mean D < δ. Superiority: Holm first step (α / 7). In brackets: sd × 1.5.",
    "",
    "| Pair | Round-1 mean D | sd | δ | N (df) | NI power at round-1 mean | at D = 0 | at D = +δ/2 | Superiority power "
    "at D = +δ | D for 80% superiority power |",
    "|---|---|---|---|---|---|---|---|---|---|",
]
for k, o in out["pairs"].items():
    md.append(
        f"| {k} | {o['round1_mean']:+.4f} | {o['sd']:.4f} | {o['margin']:.4f} | {o['N']} ({o['power']['df']}) | "
        f"{two(o, lambda p: p['noninferior']['round1_mean'])} | {two(o, lambda p: p['noninferior']['zero'])} | "
        f"{two(o, lambda p: p['noninferior']['half_margin'])} | {two(o, lambda p: p['superior_at_margin'])} | "
        f"{o['power']['superior_d80']:+.4f} ({o['power_sd_x1.5']['superior_d80']:+.4f}) |"
    )
md += [
    "",
    "Joint non-inferiority power (all seven pairs): "
    + "; ".join(
        f"{name} {joint['power'][x]:.2f} ({joint['power_sd_x1.5'][x]:.2f})"
        for x, name in (("round1_mean", "at the round-1 means"), ("zero", "at D = 0"), ("half_margin", "at D = +δ/2"))
    )
    + ".",
]
(PUB / "registration_v2/POWER.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
