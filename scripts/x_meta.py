"""Exploratory: the round-1 differences on one scale, and a random-effects meta-analysis over studies.
Common scale: D = learned - floor as a percentage of the floor's above-chance score A (mean over the pair's units of
floor - chance; chance: F1 the permutation chance of each crop, F2 0, F3 0, F4 0, F5 0.5). The locked margins were set
on the same idea (10% of the development floor's above-chance score, capped), so D and the margin are comparable across
families in these units.
Meta-analysis: units are nested in studies (data releases; the study map of x_levels.py). One-stage random-intercept
model D_us = mu + b_s + e_us, b_s ~ N(0, tau^2), e_us ~ N(0, sigma^2), fitted by REML; one within-study variance for
all studies, because five studies have one unit. Reported: mu with one-sided 95% bounds (t, studies - 1 df), tau, and
the 95% prediction interval for the true difference in a new study, mu +/- t(studies - 2) sqrt(tau^2 + se^2)
(needs three studies). Check without the variance model: mean and prediction interval of the unweighted study means.
Usage (env floors): python scripts/x_meta.py results/reports results/p4/f1_rescored.csv
  -> results/reports/X_meta.{json,md}"""

import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import optimize, stats as st
from spatialfloors.provenance import commit

PUB = Path(__file__).resolve().parents[1]
R = Path(sys.argv[1])
T1 = {k: v for k, v in json.loads((R / "P5_corrected.json").read_text())["pairs"].items() if "per_unit" in v}
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
CHANCE = dict(F2=0.0, F3=0.0, F4=0.0, F5=0.5)

# F1: floor score and permutation chance per unit (crops of one donor averaged), locked floor configuration
L1 = json.loads((PUB / "manifests/locked_v3.json").read_text())["locked"]["F1"]["expr_embed"]
assert L1["human"]["config"] == L1["mouse"]["config"]
donors = {
    f"gse303070_{d['slide'].replace('Tumor', 'Tumour')}_{d['pref']}": d["donor_cells"]
    for d in json.loads((PUB / "manifests/gse303070_donors.json").read_text())["check_vs_fov_csv"]
}
f1 = pd.read_csv(sys.argv[2]).drop_duplicates(["unit", "model", "config"], keep="last")
f1 = f1[(f1.model == "expr_embed") & (f1.config == L1["human"]["config"])].assign(u=lambda d: d.unit.map(lambda c: donors.get(c, c)))  # fmt: skip
f1 = f1.groupby("u")[["score", "chance"]].mean()


def above_chance(k, v):
    """A = mean over the pair's units of the floor's score minus chance."""
    if k.startswith("F1"):
        u = f1.loc[sorted(v["per_unit"])]
        assert abs(u.score.mean() - v["floor_mean"]) < 1e-9, (k, u.score.mean(), v["floor_mean"])
        return float((u.score - u.chance).mean()), float(u.chance.mean())
    c = CHANCE[k[:2]]
    return v["floor_mean"] - c, c


def reml(y, g):
    """Random-intercept model by REML. y = per-unit D, g = study of each unit."""
    y, g = np.asarray(y, float), np.asarray(g)
    S = sorted(set(g))
    n = np.array([(g == s).sum() for s in S], float)
    m = np.array([y[g == s].mean() for s in S])
    ssw = float(sum(((y[g == s] - y[g == s].mean()) ** 2).sum() for s in S))
    dfw = int(n.sum() - len(S))
    assert dfw > 0, "no study with two units"
    sc = y.std(ddof=1) ** 2  # the fit is done on y / sd

    def nll(p):
        s2, t2 = p
        w = n / (s2 + n * t2)
        mu = (w * m).sum() / w.sum()
        return 0.5 * (dfw * np.log(s2) + np.log(s2 + n * t2).sum() + np.log(w.sum()) + ssw / sc / s2 + (w * (m - mu) ** 2).sum() / sc)  # fmt: skip

    best = min(
        (optimize.minimize(nll, x0, method="L-BFGS-B", bounds=[(1e-8, None), (0, None)]) for x0 in ((0.5, 0.5), (1, 0), (0.1, 1))),
        key=lambda r: r.fun,
    )  # fmt: skip
    s2, t2 = best.x * sc
    w = n / (s2 + n * t2)
    mu, se, k = float((w * m).sum() / w.sum()), float(1 / np.sqrt(w.sum())), len(S)
    tq = st.t.ppf(0.95, k - 1)
    o = dict(k=k, n=int(n.sum()), mu=mu, se=se, lower=mu - tq * se, upper=mu + tq * se, tau=float(np.sqrt(t2)), sigma=float(np.sqrt(s2)))  # fmt: skip
    o["pi"] = [mu - (h := st.t.ppf(0.975, k - 2) * np.sqrt(t2 + se**2)), mu + h] if k > 2 else None
    # check: unweighted study means
    sd = m.std(ddof=1)
    o["study_means"] = dict(
        mean=float(m.mean()),
        lower=float(m.mean() - tq * sd / np.sqrt(k)),
        upper=float(m.mean() + tq * sd / np.sqrt(k)),
        pi=[float(m.mean() - (h := st.t.ppf(0.975, k - 1) * sd * np.sqrt(1 + 1 / k))), float(m.mean() + h)],
    )
    o["studies"] = {s: dict(n=int(a), mean=float(b)) for s, a, b in zip(S, n, m)}
    return o


out = dict(commit=commit(), chance="F1 permutation chance per crop; F2, F3, F4 0; F5 0.5", pairs={})
pc = lambda x, A: 100 * x / A
for k, v in T1.items():
    A, ch = above_chance(k, v)
    units = sorted(v["per_unit"])
    r = reml([v["per_unit"][u] for u in units], [study(u) for u in units])
    out["pairs"][k] = dict(
        n_units=v["n_units"],
        floor_mean=v["floor_mean"],
        chance_mean=ch,
        above_chance=A,
        registered=dict(est=v["est"], lower=v["lower"], upper=v["upper"], margin=v["margin"]),
        registered_pct={x: pc(v[x], A) for x in ("est", "lower", "upper", "margin")},
        meta=r,
        meta_pct=dict(
            mu=pc(r["mu"], A),
            lower=pc(r["lower"], A),
            upper=pc(r["upper"], A),
            tau=pc(r["tau"], A),
            pi=[pc(x, A) for x in r["pi"]] if r["pi"] else None,
            study_means={x: pc(r["study_means"][x], A) for x in ("mean", "lower", "upper")}
            | dict(pi=[pc(x, A) for x in r["study_means"]["pi"]]),
        ),
        noninferior_meta=bool(r["upper"] < v["margin"]),
        pi_below_margin=bool(r["pi"][1] < v["margin"]) if r["pi"] else None,
        check_pi_below_margin=bool(r["study_means"]["pi"][1] < v["margin"]),
    )

md = [
    "# X_meta: common scale and random-effects meta-analysis over studies (exploratory)",
    "",
    "All differences as % of the floor's above-chance score A (D = learned − floor; negative = floor ahead). "
    "δ = the locked margin on the same scale.",
    "",
    "| Pair | Units | Floor | Chance | A | D (registered) | D % of A | One-sided 95% bounds, % | δ, % |",
    "|---|---|---|---|---|---|---|---|---|",
]
for k, p in out["pairs"].items():
    q = p["registered_pct"]
    md.append(
        f"| {k} | {p['n_units']} | {p['floor_mean']:.3f} | {p['chance_mean']:+.3f} | {p['above_chance']:.3f} | "
        f"{p['registered']['est']:+.4f} | {q['est']:+.1f} | {q['lower']:+.1f}, {q['upper']:+.1f} | {q['margin']:.1f} |"
    )
md += [
    "",
    "Random-effects model over studies (REML, one within-study variance). μ = mean difference over studies; τ = "
    "between-study SD; PI = 95% prediction interval for the true difference in a new study (three studies needed). "
    "Check = unweighted study means, mean and 95% prediction interval for a new study's mean. All in % of A.",
    "",
    "| Pair | Studies | Units | μ | One-sided 95% bounds | τ | 95% PI | Non-inferior (μ) | PI below δ | Check: mean | Check: PI | Check: PI below δ |",
    "|---|---|---|---|---|---|---|---|---|---|---|---|",
]
yn = lambda x: "–" if x is None else ("yes" if x else "no")
for k, p in out["pairs"].items():
    q, r = p["meta_pct"], p["meta"]
    md.append(
        f"| {k} | {r['k']} | {r['n']} | {q['mu']:+.1f} | {q['lower']:+.1f}, {q['upper']:+.1f} | {q['tau']:.1f} | "
        + (f"{q['pi'][0]:+.1f}, {q['pi'][1]:+.1f}" if q["pi"] else "–")
        + f" | {yn(p['noninferior_meta'])} | {yn(p['pi_below_margin'])} | {q['study_means']['mean']:+.1f} | "
        f"{q['study_means']['pi'][0]:+.1f}, {q['study_means']['pi'][1]:+.1f} | {yn(p['check_pi_below_margin'])} |"
    )
md += ["", "Study means of D (units):", ""]
for k, p in out["pairs"].items():
    md.append(f"- {k}: " + "; ".join(f"{s} {x['mean']:+.4f} ({x['n']})" for s, x in p["meta"]["studies"].items()))
(R / "X_meta.json").write_text(json.dumps(out, indent=1, default=float))
(R / "X_meta.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
print("X_META_DONE")
