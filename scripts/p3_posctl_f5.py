"""Positive control for F5 learned models: the planted-type toy of
tests/test_f5_floor.py (4 fine labels, query depth x0.3 and per-gene platform factor). Pass: coarse macro-F1 >= .90
for scANVI (n_latent 10, nb) and Tangram (clusters, rna_count_based); CPU. Floors on the same toy: > .9 (tests).
Usage (env scvi): python scripts/p3_posctl_f5.py -> runs/p3/posctl_f5.json"""

import json, sys, time
from pathlib import Path

sys.path.insert(0, "tests")
from test_f5_floor import _toy
from spatialfloors.floors import f5
from spatialfloors.models import transfer
from spatialfloors.metrics.f5 import macro_f1_gated
from spatialfloors.provenance import stamp

Xr, yr, Xq, cq, _ = _toy()
rep = {}
for name, fn in (
    ("scanvi", lambda: transfer.scanvi(Xr, yr, Xq, accelerator="cpu", max_epochs=100, scanvi_epochs=20)),
    ("tangram", lambda: transfer.tangram(Xr, yr, Xq, device="cpu")),
):
    t = time.time()
    p = f5.coarse(fn())
    f = macro_f1_gated(cq, p, ("epithelial", "immune_nonmyeloid", "myeloid"))[0]
    rep[name] = dict(macro_f1=f, s=round(time.time() - t, 1), passed=f >= 0.9)
    print(name, rep[name], flush=True)
Path("runs/p3/posctl_f5.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
