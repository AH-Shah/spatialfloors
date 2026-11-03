"""P3 positive / negative controls for F3 (decision 2026-10-28: planted toys, option b). sim.niche_toy cases:
  quadratic  (|type-1 neighbours' state|): learned model must beat the ridge floor on receiver 0, PC1 rho >= floor + .08
             (floor blind by construction; LightGBM ceiling probe .12, scratch 2026-10-28)
  nonlinear  (sender state x third-type presence): reported only; not recoverable even with an oracle product feature
             (ridge on [features, z1 x p2] .03), so no pass criterion
  couple     (linear): learned PC1 rho >= .15 (finds what the floor finds)
  bleed      (receiver RNA in neighbours only): learned |unit score| < .02 and |PC1| < .10 (no leak)
Models: ridge floor, LightGBM (leaves 15, lr .05), GAT (2 layers, hidden 128; device from argv, default cpu).
Usage (env scvi): python scripts/p3_posctl_f3.py [cpu|mps] -> runs/p3/posctl_f3.json"""
import json, sys, time
from pathlib import Path
import numpy as np
from spatialfloors.sim import niche_toy
from spatialfloors.floors import f3
from spatialfloors.models import gat
from spatialfloors.provenance import stamp

dev = sys.argv[1] if len(sys.argv) > 1 else "cpu"; rep = {}


def lgbm(Xtr, Ytr, Xte):
    from lightgbm import LGBMRegressor
    return np.column_stack([LGBMRegressor(num_leaves=15, learning_rate=.05, n_estimators=100, n_jobs=4, random_state=0, verbose=-1).fit(Xtr, Ytr[:, j]).predict(Xte) for j in range(Ytr.shape[1])])


for case in ("quadratic", "nonlinear", "couple", "bleed"):
    C, xy, fov, ids = niche_toy(case); D = f3.prepare(C, xy, fov, ids=ids); rep[case] = {}
    for name, fn in (("ridge", lambda: f3.score_floor(D, 25., 30)), ("lightgbm", lambda: f3.score_floor(D, 25., 30, fit_predict=lgbm)),
                     ("gat", lambda: gat.score(D, 25., 2, 128, 30, device=dev)[0])):
        t = time.time(); s = fn(); rep[case][name] = dict(unit=f3.unit_score(s), r0_pc1=s[0][2][0], r0=s[0][0], s=round(time.time() - t, 1))
        print(case, name, {k: round(v, 3) for k, v in rep[case][name].items()}, flush=True)
R = rep; chk = {}
for m in ("lightgbm", "gat"):
    chk[m] = dict(quadratic=R["quadratic"][m]["r0_pc1"] >= R["quadratic"]["ridge"]["r0_pc1"] + .08, couple=R["couple"][m]["r0_pc1"] >= .15,
                  bleed=abs(R["bleed"][m]["unit"]) < .02 and abs(R["bleed"][m]["r0_pc1"]) < .10)
    chk[m]["pass"] = all(chk[m].values())
rep["checks"] = chk; print(json.dumps(chk)); Path("runs/p3/posctl_f3.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
