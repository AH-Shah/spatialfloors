"""F2 v3 acceptance helper (PLAN §19): per-region recall (Hungarian match) averaged over all GMM seeds, from the labels
scripts/v3_f2.py saved (the run's per_region uses seed 0 only). Prints, per unit and partition, ari_bal and the mean
recall per region, and the largest drop vs v2. Usage (env floors): python scripts/v3_f2_regions.py <unit> ..."""
import json, os, sys
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn import metrics as M
from spatialfloors.stress.f2 import load, Panel

for u in sys.argv[1:]:
    P = Panel(load(u)); O = os.environ.get("SF_OUT", "runs/v3"); z = np.load(f"{O}/f2_{u}_labels.npz"); rep = json.loads(Path(f"{O}/f2_{u}.json").read_text())["partitions"]; regs = np.unique(P.y)
    def rec(L):
        out = []
        for l in L:
            C = M.cluster.contingency_matrix(P.y, np.asarray(l).astype(str)[P.m]); r, c = linear_sum_assignment(-C); d = {regs[i]: C[i, j] / C[i].sum() for i, j in zip(r, c)}
            out.append([d.get(x, 0.) for x in regs])
        return np.mean(out, 0)
    R = {k: rec(z[k]) for k in z.files}; base = R["v2"]
    print("==", u, list(regs))
    for k, v in R.items(): print(f"  {k:26} bal {rep[k]['ari_bal']:.3f} bnd {rep[k]['ari_boundary']:.3f} stab {rep[k]['seed_stability']:.3f} | rec {np.round(v, 2).tolist()} | worst drop vs v2 {float((v - base).min()):+.2f}")
