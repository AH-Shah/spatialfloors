"""Print the run parameters of the val-best run per gs_id of the NCEM MIBI-TOF reproduction (the configuration ported into
F3). Env external/ncem_env. Usage: python scripts/ncem_repro_best_params.py <gs_id ...> -> stdout"""
import pickle, sys
import ncem

gs = ncem.train.GridSearchContainer("runs/p3/ncem_repro/gs/", gs_ids=sys.argv[1:], lateral_resolution=400 / 1024); gs.load_gs()
t = gs.summary_table
for g in sys.argv[1:]:
    s = t[t.gs_id == g]; m = s.groupby("run_id")[["val_r_squared_linreg", "test_r_squared_linreg"]].mean()
    for b in m.sort_values("val_r_squared_linreg", ascending=False).index[:3]:
        mid = s[s.run_id == b].model_id.iloc[0]; d = s[s.run_id == b].iloc[0]
        print(g, b, round(m.loc[b, "val_r_squared_linreg"], 4), round(m.loc[b, "test_r_squared_linreg"], 4))
        p = pickle.load(open(f"runs/p3/ncem_repro/gs/{g.upper()}/results/{mid}_runparams.pickle", "rb"))
        print("  ", {k: v for k, v in p.items() if k not in ("data_set", "gs_id", "model_id", "graph_covar_selection")})
