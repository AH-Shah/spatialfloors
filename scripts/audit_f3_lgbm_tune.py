"""Audit (lock v2, user 2026-11-02 'tune LightGBM properly'): is LightGBM's F3 deficit a hyperparameter-range artefact?
The locked config (15 leaves, lr .05) is the most regularised corner of the P3 grid {15, 63} x {.05, .1}, with LightGBM's
default min_child_samples 20, no row / feature subsampling and a fixed 100 rounds: the optimum may lie outside the grid.
Every variant: lr .05 (or .01), up to 3000 rounds, early stopping (50) on ~15% of whole training FOVs (f3.groups: spatially
separate inner validation), one model per target PC, on the v2 floor's features (human T+ctx ell 25 pc10, mouse ctx
ell 100 pc50), same block_cv as the floor. Variants span leaves x min_child_samples x L2 x subsampling, plus
linear_tree (ridge-like leaves: a LightGBM that contains the linear floor).
Usage (env scvi): python scripts/audit_f3_lgbm_tune.py <unit> [variant ...] -> runs/p3v2/audit_lgbm_tune_<unit>.json"""
import json, os, sys, time
from pathlib import Path
import anndata as ad, numpy as np
from spatialfloors.floors import f3
from spatialfloors.metrics import f3 as m3
from spatialfloors.provenance import stamp

B = dict(subsample=1., subsample_freq=0, colsample_bytree=1., reg_lambda=0.)
SS = dict(subsample=.8, subsample_freq=1, colsample_bytree=.8)
V = {"l4_m200": dict(num_leaves=4, min_child_samples=200), "l8_m200": dict(num_leaves=8, min_child_samples=200),
     "l16_m200": dict(num_leaves=16, min_child_samples=200), "l8_m50": dict(num_leaves=8, min_child_samples=50),
     "l8_m800": dict(num_leaves=8, min_child_samples=800), "l8_m200_l2": dict(num_leaves=8, min_child_samples=200, reg_lambda=10.),
     "l8_m200_ss": dict(num_leaves=8, min_child_samples=200, **SS), "l31_m800_l2_ss": dict(num_leaves=31, min_child_samples=800, reg_lambda=10., **SS),
     "lin_l4_m800": dict(num_leaves=4, min_child_samples=800, linear_tree=True, linear_lambda=10.),
     "l4_m200_lr.01": dict(num_leaves=4, min_child_samples=200, learning_rate=.01), "l8_m800_lr.01": dict(num_leaves=8, min_child_samples=800, learning_rate=.01),
     "l8_m200_ss_lr.01": dict(num_leaves=8, min_child_samples=200, learning_rate=.01, **SS),
     "lin_l4_m800_lr.01": dict(num_leaves=4, min_child_samples=800, linear_tree=True, linear_lambda=10., learning_rate=.01),
     # round 2: rounds as a grid choice (inner early stopping on ~15% FOVs stops at round 1 on mouse: 16% constant PCs)
     "fx_l15_m20_r100": dict(num_leaves=15, min_child_samples=20, rounds=100), "fx_l8_m200_ss_r100": dict(num_leaves=8, min_child_samples=200, rounds=100, **SS),
     "fx_l8_m200_ss_r300": dict(num_leaves=8, min_child_samples=200, rounds=300, **SS), "fx_l4_m200_ss_r300": dict(num_leaves=4, min_child_samples=200, rounds=300, **SS),
     "fx_l8_m50_ss_r300": dict(num_leaves=8, min_child_samples=50, rounds=300, **SS), "fx_l8_m200_ss_lr.02_r500": dict(num_leaves=8, min_child_samples=200, rounds=500, learning_rate=.02, **SS),
     "fx_l15_m200_ss_l2_r300": dict(num_leaves=15, min_child_samples=200, reg_lambda=10., rounds=300, **SS), "fx_lin_l4_m200_r100": dict(num_leaves=4, min_child_samples=200, linear_tree=True, linear_lambda=10., rounds=100)}
u = sys.argv[1]; sel = sys.argv[2:] or list(V); out = Path(f"runs/p3v2/audit_lgbm_tune_{u}.json")
sp = "mouse" if u.startswith("allen") else "human"; F_ = dict(human=("T+ctx", 25., 10), mouse=("ctx", 100., 50))[sp]
a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); p = json.loads(Path(f"manifests/panels/{u}.json").read_text()); g = a.var_names.get_indexer(p["genes"])
D = f3.prepare(a.X[:, g], a.obs[["x_um", "y_um"]].values, a.obs.fov.values); nb, w = f3.kernel(D, F_[1]); T = f3.id_mix(D, nb, w)
rep = json.loads(out.read_text()) if out.exists() else dict(unit=u, features=F_, grid={})


def run(prm):
    from lightgbm import LGBMRegressor, early_stopping
    res = {}
    for r in f3.receivers(D):
        rows, Y, dep = f3.target(D, r); ctx = f3.context_fn(D, r, rows, nb, w, F_[2]); grp = f3.groups(D["fov"][rows], D["xy"][rows])
        X = {"ctx": ctx, "T+ctx": lambda q, tr: np.c_[T[rows], ctx(q, tr)]}[F_[0]]; nr = []

        def fp(Xtr, Ytr, Xte, tr, te, q):
            fv = np.unique(grp[tr]); v = np.isin(grp[tr], np.random.default_rng(0).choice(fv, max(1, round(.15 * len(fv))), replace=False)); P = []
            q_ = dict(prm); R = q_.pop("rounds", None)  # fixed rounds (grid choice, no inner early stopping; all training rows)
            for j in range(Ytr.shape[1]):
                m = LGBMRegressor(n_estimators=R or 3000, n_jobs=4, random_state=0, verbose=-1, **{**B, "learning_rate": .05, **q_})
                if R: m.fit(Xtr, Ytr[:, j]); nr.append(R)
                else: m.fit(Xtr[~v], Ytr[~v, j], eval_set=[(Xtr[v], Ytr[v, j])], callbacks=[early_stopping(50, verbose=False)]); nr.append(int(m.best_iteration_ or 0))
                P.append(m.predict(Xte))
            return np.column_stack(P)
        fp.wants_rows = True
        s, fo = m3.block_cv(X, Y, D["xy"][rows], fp, f3.RF_UM, blocks=D["blocks"][rows], split=D["split"], covars=dep, detail=(det := {})); pp = np.array(det["per_pc"], float); res[r] = (s, fo, [np.nan], nr, float(np.nan_to_num(pp).mean()), float(np.isnan(pp).mean()))
    return dict(unit=f3.unit_score({k: v[:3] for k, v in res.items()}), per_receiver={int(k): dict(score=v[0], rounds_median=float(np.median(v[3])), score_nan0=v[4], nan_frac=v[5]) for k, v in res.items()}, unit_nan0=float(np.mean([v[4] for v in res.values()])))


for k in sel:
    t = time.time(); rep["grid"][k] = dict(run(V[k]), s=round(time.time() - t)); rr = [v["rounds_median"] for v in rep["grid"][k]["per_receiver"].values()]
    print(u, k, round(rep["grid"][k]["unit"], 4), "nan0", round(rep["grid"][k]["unit_nan0"], 4), "nanfrac", round(np.mean([v["nan_frac"] for v in rep["grid"][k]["per_receiver"].values()]), 3), "median rounds", np.median(rr), rep["grid"][k]["s"], "s", flush=True)
    out.write_text(json.dumps(stamp(rep), indent=1, default=float))
print("TUNE_DONE", u)
