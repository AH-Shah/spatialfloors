"""Positive control for F2 learned models: sim.domain_toy (2 domains differing
only in neighbourhood composition). Pass: CellCharter (n_layers 3, n_latent 10; scVI on CPU, 100 epochs) ARI >= .80
at k = 2, clustering seeds 0-2 (mean). Reference: floor .98, expression-only k-means .25 (tests/test_f2_floor.py).
Novae cannot read synthetic genes (pretrained gene vocabulary); its control is the authors' tutorial run.
Usage (env cellcharter): python scripts/p3_posctl_f2.py -> runs/p3/posctl_f2.json"""

import json, time
from pathlib import Path
import numpy as np
from sklearn.metrics import adjusted_rand_score
from spatialfloors.sim import domain_toy
from spatialfloors.models import cellcharter
from spatialfloors.provenance import stamp

C, xy, fov, dom = domain_toy()
t = time.time()
aris = []
for s in range(3):
    lab, info = cellcharter.fit(C, xy, 2, 10, 3, accelerator="cpu", seed=s, max_epochs=100)
    aris.append(float(adjusted_rand_score(dom, lab)))
rep = dict(model="cellcharter", ari=aris, mean=float(np.mean(aris)), info=info, s=round(time.time() - t, 1))
rep["pass"] = rep["mean"] >= 0.8
print(json.dumps(rep), flush=True)
Path("runs/p3/posctl_f2.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
