"""Audit 2026-11-02: confirmation power per (family, model) pair for floor non-inferiority (one-sided 95% upper bound <
margin, df = N - 3 strata), with sd from development + validation D (231w excluded, as v_gate_report G6), at three true
mean D: the observed mean, 0 (learned = floor) and +margin/2 (learned slightly better). Joint = product over pairs (the
intersection-union claim needs every pair). Usage: python scripts/audit_power.py -> prints"""
import json
from pathlib import Path
import numpy as np, scipy.stats as st

P = json.loads(Path("runs/p3/P5_gate.json").read_text())["pairs"]; P.update(json.loads(Path("runs/v_retro/P5_retro.json").read_text())["pairs"])
NC = dict(F1=25, F2=10, F3=24, F5=20, F4=10); joint = {"obs": 1., "zero": 1., "half": 1.}
for k, v in P.items():
    f = k.split("|")[0].split("_")[0]; D = np.array([d for u, d in v["per_unit"].items() if not u.endswith("w")])
    if len(D) < 3: print(k, "n<3, skipped"); continue
    n, sd, m = NC[f], D.std(ddof=1), v["margin"]; df = n - 3; se = sd / np.sqrt(n); tq = st.t.ppf(.95, df)
    pw = {lab: float(st.nct.cdf(-tq, df, (mu - m) / se)) for lab, mu in (("obs", D.mean()), ("zero", 0.), ("half", m / 2))}
    for lab in joint: joint[lab] *= pw[lab]
    print(f"{k:22s} N {n:2d} sd {sd:.4f} margin {m:.4f} | power at mean D {D.mean():+.4f}: {pw['obs']:.2f}  at D=0: {pw['zero']:.2f}  at D=+m/2: {pw['half']:.2f}")
print("joint (all pairs pass):", {k: round(x, 2) for k, x in joint.items()})
