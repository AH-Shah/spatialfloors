"""Planted-signal controls for the authors'-setup NCEM port (models.ncem.grid, HYPOTHESES row 375) before any slide:
compo (receiver state linear in the local type-1 share: what type-label NCEM can see), couple (linear in the sender
state: only type+ctx labels can see it), bleed (no signal, neighbour marker bleed). Each model is checked against ridge
on exactly its own inputs (models.ncem.input_ridge) and the toy floor (T+ctx ridge, 30 PCs). Pass rules: compo r0 PC1
>= .15 and >= its input ridge - .03 (type-label models); couple r0 PC1 >= .15 (type+ctx); bleed |unit| < .02 and
|r0 PC1| < .10 (all). Env external/ncem_env, PYTHONPATH=src. Usage: external/ncem_env/bin/python
scripts/p3_posctl_ncem_auth.py [radius_um] -> runs/p3/posctl_ncem_auth_r<R>.json"""
import json, sys, time
from pathlib import Path
from spatialfloors.sim import niche_toy
from spatialfloors.floors import f3
from spatialfloors.models import ncem
from spatialfloors.provenance import stamp

R = float(sys.argv[1]) if len(sys.argv) > 1 else 46.; G = ncem.grid((R,)); tag = f"r{R:g}"
M = {"edA": G[f"edA_e1_{tag}"], "linA": G[f"linA_l20_{tag}"], "edT_type": G[f"edT_type_{tag}"], "edT_ctx": dict(G[f"edT_type+ctx_{tag}"], n_pc=30)}
INP = {"edA": "mean", "linA": "presence", "edT_type": "mean"}
rep = {}
for case in ("compo", "couple", "bleed"):
    C, xy, fov, ids = niche_toy(case); D = f3.prepare(C, xy, fov, ids=ids); rep[case] = {}
    runs = [("floor", lambda: (f3.score_floor(D, 25., 30), {})), ("ridge_mean", lambda: (ncem.input_ridge(D, R, "mean"), {})), ("ridge_presence", lambda: (ncem.input_ridge(D, R, "presence"), {}))]
    runs += [(k, (lambda kw=kw: ncem.score(D, 25., **kw))) for k, kw in M.items()]
    for name, fn in runs:
        t = time.time(); s, ep = fn(); rep[case][name] = dict(unit=f3.unit_score(s), r0_pc1=s[0][2][0], r0=s[0][0], s=round(time.time() - t, 1), train={int(k): v for k, v in ep.items()})
        print(case, name, {k: round(v, 3) for k, v in rep[case][name].items() if k != "train"}, {k: (v["epochs"], [[round(x, 2) for x in l] for l in v["val"]]) for k, v in ep.items()} if ep else "", flush=True)
X = rep; chk = {}
for m in M:
    c = dict(bleed=abs(X["bleed"][m]["unit"]) < .02 and abs(X["bleed"][m]["r0_pc1"]) < .10)
    if m in INP: c["compo"] = X["compo"][m]["r0_pc1"] >= .15 and X["compo"][m]["r0_pc1"] >= X["compo"]["ridge_" + INP[m]]["r0_pc1"] - .03
    else: c["couple"] = X["couple"][m]["r0_pc1"] >= .15
    c["pass"] = all(c.values()); chk[m] = c
rep["checks"] = chk; print(json.dumps(chk)); Path(f"runs/p3/posctl_ncem_auth_{tag}.json").write_text(json.dumps(stamp(rep), indent=1, default=float)); print("POSCTL_DONE")
