"""Evaluate the NCEM MIBI-TOF reproduction with the authors' own GridSearchContainer (as model_benchmarks_mibitof_cancer.ipynb):
per gs_id, best run by val r_squared_linreg, its test r_squared_linreg per CV, and best per radius (um = px * 400 / 1024).
Published targets: ED ~.69, ED-NCEM GCN ~.70-.715 at 13-46 um, linear interactions ~.80.
Env external/ncem_env. Usage: python scripts/ncem_repro_eval.py <gs_id ...> -> stdout"""
import sys
import numpy as np
import ncem

gs = ncem.train.GridSearchContainer("runs/p3/ncem_repro/gs/", gs_ids=sys.argv[1:], lateral_resolution=400 / 1024); gs.load_gs()
t = gs.summary_table; print(t.columns.tolist()[:40])
for g in sys.argv[1:]:
    s = t[t.gs_id == g]
    if not len(s): continue
    m = s.groupby("run_id")[["val_r_squared_linreg", "test_r_squared_linreg"]].mean(); best = m.val_r_squared_linreg.idxmax()
    print(g, "best (val)", best, "test per cv", np.round(s[s.run_id == best].test_r_squared_linreg.values, 4), "mean", round(m.loc[best, "test_r_squared_linreg"], 4))
    if "um_radius" in s:
        for r, sr in s.groupby("um_radius"):
            mr = sr.groupby("run_id")[["val_r_squared_linreg", "test_r_squared_linreg"]].mean(); b = mr.val_r_squared_linreg.idxmax()
            print(f"  radius {r:7.1f} um  test {mr.loc[b, 'test_r_squared_linreg']:.4f}  (val {mr.loc[b, 'val_r_squared_linreg']:.4f}, {len(mr)} runs)")
