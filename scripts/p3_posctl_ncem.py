"""F3 planted-signal controls for NCEM (PLAN §18; same toys and pass rules as scripts/p3_posctl_f3.py, row 349):
quadratic r0 PC1 >= ridge + .08; couple r0 PC1 >= .15; bleed |unit| < .02 and |r0 PC1| < .10; nonlinear reported only.
NCEM at lr 5e-3 with neighbour labels 'type' and 'type+ctx' (30 context PCs as the toy floor). Env external/ncem_env, PYTHONPATH=src.
Usage: external/ncem_env/bin/python scripts/p3_posctl_ncem.py -> runs/p3/posctl_ncem.json"""
import json, time
from pathlib import Path
from spatialfloors.sim import niche_toy
from spatialfloors.floors import f3
from spatialfloors.models import ncem
from spatialfloors.provenance import stamp

rep = {}
for case in ("quadratic", "nonlinear", "couple", "bleed"):
    C, xy, fov, ids = niche_toy(case); D = f3.prepare(C, xy, fov, ids=ids); rep[case] = {}
    runs = [("ridge", lambda: f3.score_floor(D, 25., 30))] + [(f"ncem_{h}", (lambda h=h: ncem.score(D, 25., latent_dim=30, cond_dim=64, lr=5e-3, batch_size=4, h0=h, n_pc=30)[0])) for h in ("type", "type+ctx")]  # grid's high-capacity config
    for name, fn in runs:
        t = time.time(); s = fn(); rep[case][name] = dict(unit=f3.unit_score(s), r0_pc1=s[0][2][0], r0=s[0][0], s=round(time.time() - t, 1))
        print(case, name, {k: round(v, 3) for k, v in rep[case][name].items()}, flush=True)
R = rep; chk = {}
for m in ("ncem_type", "ncem_type+ctx"):
    chk[m] = dict(quadratic=R["quadratic"][m]["r0_pc1"] >= R["quadratic"]["ridge"]["r0_pc1"] + .08, couple=R["couple"][m]["r0_pc1"] >= .15,
                  bleed=abs(R["bleed"][m]["unit"]) < .02 and abs(R["bleed"][m]["r0_pc1"]) < .10)
    chk[m]["pass"] = all(chk[m].values())
rep["checks"] = chk; print(json.dumps(chk)); Path("runs/p3/posctl_ncem.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
