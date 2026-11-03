"""Predicted confirmation power per test pair at the locked v3 configurations (PROTOCOL §4). Per-unit D = learned - floor
on development + validation units (whole slide 231w excluded; F4: the five retrospective units), sd of D, N = the
family's confirmation units, df = N - strata; floor non-inferior if the one-sided 95% upper bound of mean D < margin.
Power at the observed mean D, at D = 0 and at D = +margin/2; joint = product over pairs (every pair must pass).
Usage: python scripts/power_v3.py -> stdout + manifests/power_v3.json"""
import json
from pathlib import Path
import numpy as np, scipy.stats as st

J = lambda p: json.loads(Path(p).read_text()); L = J("manifests/locked_v3.json")["locked"]; NC = J("manifests/locked_v3_ncem.json")["locked"]
N = dict(F1=(25, 4), F2=(10, 2), F3=(24, 3), F4=(10, 2), F5=(20, 2))  # confirmation units, strata
DEV_H = ["gse1_232", "gse0_S0", "gse303070_148191A_118N", "gse303070_148191A_278N", "gse303070_148191A_278T", "gse303070_148191A_46T"]
VAL, A4 = ["gse1_231", "bruker_breast"], "allen_animal4_C57BL6J-4.003"
donor = lambda u: {"gse303070_148191A_278N": "Donor278", "gse303070_148191A_278T": "Donor278", "gse303070_148191A_118N": "Donor118", "gse303070_148191A_46T": "Donor46"}.get(u, u)
sp = lambda u: "mouse" if u.startswith("allen") else "human"


def per_unit(d):
    out = {}
    for u, v in d.items(): out.setdefault(donor(u), []).append(v)
    return {k: float(np.mean(v)) for k, v in out.items()}


D = {}
for src in ("runs/p3v2/P5_gate.json", "runs/v_retro/P5_retro.json"):  # F1, F2 (configs unchanged in v3), F4 retrospective
    for k, v in J(src)["pairs"].items():
        if k.split("|")[0] in ("F1", "F2", "F4_pearson"): D[k] = {u: d for u, d in v["per_unit"].items() if not u.endswith("w")}
f3fl = {u: J(f"runs/v3/f3alt_{u}.json")["unit_mean"]["fixed_sf"] for u in DEV_H + [A4]} | {u: J(f"runs/v3/val/f3_{u}_ridge.json")["grid"]["v3"]["unit"] for u in VAL}
lg = {u: J(f"runs/v3/lgbm/f3_{u}_lightgbm.json")["grid"]["leaves15_lr0.05"]["unit"] for u in f3fl}
nc = {u: J(f"runs/p3v3/f3_{u}_ncemA.json")["grid"][NC[sp(u)]["config"]]["unit"] for u in f3fl if Path(f"runs/p3v3/f3_{u}_ncemA.json").exists()}
D["F3|lightgbm"] = per_unit({u: lg[u] - f3fl[u] for u in lg}); D["F3|ncem"] = per_unit({u: nc[u] - f3fl[u] for u in nc})
f5fl = {u: J(f"runs/v3/f5_{u}.json")["summaries"]["v2+spatial"]["heldout_auroc"] for u in DEV_H if u != "gse0_S0"} | {u: J(f"runs/v3/val/f5_{u}_floors.json")["grid"]["insitutype_spatial"]["heldout_auroc"] for u in VAL}
for m in ("scanvi", "tangram"):
    le = {u: J(f"runs/p3v2/f5_{u}_{m}.json")["grid"][L["F5"][m]["human"]["config"]]["heldout_auroc"] for u in f5fl if Path(f"runs/p3v2/f5_{u}_{m}.json").exists()}
    D[f"F5|{m}"] = per_unit({u: le[u] - f5fl[u] for u in le})
out, joint = {}, {"obs": 1., "zero": 1., "half": 1.}
for k in sorted(D):
    f = k.split("|")[0].split("_")[0]; x = np.array(list(D[k].values())); m = L[f]["margin"]; n, s = N[f]
    if len(x) < 3: out[k] = dict(n_dev=len(x), note="fewer than 3 development / validation units: not estimated"); print(k, "n < 3"); continue
    df = n - s; se = x.std(ddof=1) / np.sqrt(n); tq = st.t.ppf(.95, df)
    pw = {lab: float(st.nct.cdf(-tq, df, (mu - m) / se)) for lab, mu in (("obs", x.mean()), ("zero", 0.), ("half", m / 2))}
    for lab in joint: joint[lab] *= pw[lab]
    out[k] = dict(n_dev=len(x), mean_D=float(x.mean()), sd_D=float(x.std(ddof=1)), margin=m, N=n, df=df, power=pw, per_unit=D[k])
    print(f"{k:20s} dev/val n {len(x)} mean D {x.mean():+.4f} sd {x.std(ddof=1):.4f} | N {n} margin {m:.4f} | power obs {pw['obs']:.2f} D=0 {pw['zero']:.2f} D=+m/2 {pw['half']:.2f}")
out["joint"] = joint; print("joint:", {k: round(v, 2) for k, v in joint.items()})
Path("manifests/power_v3.json").write_text(json.dumps(out, indent=1) + "\n")
