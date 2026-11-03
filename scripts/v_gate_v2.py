"""Gate re-run for lock v2 (PLAN §18 procedure step 4): G2 floor landscape at the v2 floors (tolerance = the family
margin, §18) on development + validation units; G6 power and G7 crop check from <R>/P5_gate.json (v2 configs).
Same rules as scripts/v_gate_report.py otherwise. Usage (env floors): python scripts/v_gate_v2.py [runs/p3v2] -> stdout"""
import json, sys
from pathlib import Path
import numpy as np
from scipy import stats as st

R = Path(sys.argv[1] if sys.argv[1:] else "runs/p3v2"); L = json.loads(Path("manifests/locked_v2.json").read_text())["locked"]
U = ["gse1_232", "gse0_S0"] + [f"gse303070_148191A_{s}" for s in ("118N", "278N", "278T", "46T")] + ["allen_animal4_C57BL6J-4.003", "gse1_231", "gse1_231w", "bruker_breast"]
VAL = ("gse1_231", "gse1_231w", "bruker_breast")
sp = lambda u: "mouse" if u.startswith("allen") else "human"
F3ID = {"human": "grid|T+ctx|ell25|pc10|a1000.0", "mouse": "grid|ctx|ell100|pc50|a1000.0"}
FAM = dict(f1=("excess_k15", lambda s: "grid|pearson|none|pc50|none|split0", L["F1"]["margin"]),
           f2=("ari_bal", lambda s: "grid|gmm_full|ell400|lam0.8|pc50", L["F2"]["margin"]),
           f3=("top10", lambda s: F3ID[s], L["F3"]["margin"]))
print("## G2 v2 (locked floor vs best grid config; tolerance = margin)")
for f, (p, lid, eps) in FAM.items():
    for u in U:
        q = Path(f"runs/landscape/{f}_{u}.jsonl")
        if not q.exists(): continue
        Rw = {r["id"]: r for r in map(json.loads, open(q)) if r.get("part") == "grid" and r.get(p) == r.get(p)}
        if not Rw: continue
        b = max(Rw.values(), key=lambda r: r[p]); i = lid(sp(u)); lv = Rw[i][p] if i in Rw else np.nan
        print(f"{f} {u[-14:]:>14} {'VAL' if u in VAL else 'dev'} locked {lv:.4f} best {b[p]:.4f} ({b['id'][5:]}) gap {b[p] - lv:.4f} {'ok' if b[p] - lv <= eps else 'OFF'}")
p5 = R / "P5_gate.json"
if p5.exists():
    P = json.loads(p5.read_text())["pairs"]; NC = dict(F1=25, F2=10, F3=24, F5=20)
    print("## G6 power v2 (dev + validation D, 231w excluded) / G7 crop check (231 crop vs 231w whole)")
    for k, v in P.items():
        f = k.split("|")[0].split("_")[0]; D = np.array([d for u, d in v["per_unit"].items() if not u.endswith("w")]); n = NC.get(f)
        if n is None or len(D) < 3: continue
        sd = D.std(ddof=1); se = sd / np.sqrt(n); df = n - 3; tq = st.t.ppf(.95, df)
        pw = lambda mu: float(st.nct.cdf(-tq, df, (mu - v["margin"]) / se))
        a, w = v["per_unit"].get("gse1_231"), v["per_unit"].get("gse1_231w")
        g7 = "" if a is None or w is None else f" | G7 crop {a:+.4f} whole {w:+.4f} {'pass' if abs(w - a) <= v['margin'] else 'FAIL'}"
        print(f"{k}: n {len(D)} mean D {D.mean():+.4f} sd {sd:.4f} margin {v['margin']:.4f} power@obs {pw(D.mean()):.2f} @0 {pw(0.):.2f}{g7}")
