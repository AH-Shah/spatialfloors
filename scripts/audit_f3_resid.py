"""Audit 2026-11-02 (pre-lock review, development units only): is LightGBM's F3 deficit vs the ridge floor an
under-tuned learner? Variants on the floor's own features [T, ctx] (ell 25 um, 50 PCs), same block_cv:
  lgbm_es    LightGBM 15 leaves, lr .05, up to 1000 rounds, early stopping (30) on a random 15% of training rows
  lgbm_resid ridge floor prediction + lgbm_es fitted on the ridge's training residuals (boosting from the floor)
Locked ridge / LightGBM scores come from runs/p3/f3_<unit>_{ridge,lightgbm}.json.
Usage (env scvi): python scripts/audit_f3_resid.py <unit ...> -> runs/p3/audit_f3_resid_<unit>.json"""
import json, sys, time
from pathlib import Path
import anndata as ad, numpy as np
from spatialfloors.floors import f3
from spatialfloors.floors.f3 import groups
from spatialfloors.provenance import stamp


def lgbm_es(Xtr, Ytr, Xte, base_tr=None, base_te=None, v=None):
    from lightgbm import LGBMRegressor, early_stopping
    rng = np.random.default_rng(0); v = rng.random(len(Xtr)) < .15 if v is None else v; out = []
    for j in range(Ytr.shape[1]):
        y = Ytr[:, j] - (base_tr[:, j] if base_tr is not None else 0)
        m = LGBMRegressor(num_leaves=15, learning_rate=.05, n_estimators=1000, n_jobs=4, random_state=0, verbose=-1)
        m.fit(Xtr[~v], y[~v], eval_set=[(Xtr[v], y[v])], callbacks=[early_stopping(30, verbose=False)])
        out.append(m.predict(Xte) + (base_te[:, j] if base_te is not None else 0))
    return np.column_stack(out)


def resid(Xtr, Ytr, Xte):
    from sklearn.linear_model import RidgeCV
    r = RidgeCV(alphas=f3.ALPHAS).fit(Xtr, Ytr); return lgbm_es(Xtr, Ytr, Xte, r.predict(Xtr), r.predict(Xte))


def score_fov(D, resid_on_ridge, ell=25., feats="T+ctx", locked=False):
    """As f3.score_floor (50 PCs) but early stopping on whole training FOVs (~15%; f3.groups): spatially separate
    inner validation (SF_AUDIT_FOV=1; the random-row split leaks spatial autocorrelation). locked: the P3 LightGBM
    (15 leaves, lr .05, 100 rounds, no early stopping) on these features. ell / feats: SF_AUDIT_ELL, SF_AUDIT_FEATS."""
    from sklearn.linear_model import RidgeCV
    from spatialfloors.metrics import f3 as m3
    nb, w = f3.kernel(D, ell); T = f3.id_mix(D, nb, w); out = {}
    for r in f3.receivers(D):
        rows, Y, dep = f3.target(D, r); ctx = f3.context_fn(D, r, rows, nb, w, 50); fov = groups(D["fov"][rows], D["xy"][rows])
        X = (lambda q, tr_: np.c_[T[rows], ctx(q, tr_)]) if feats == "T+ctx" else ctx

        def fp(Xtr, Ytr, Xte, tr, te, q):
            if locked:
                from lightgbm import LGBMRegressor
                return np.column_stack([LGBMRegressor(num_leaves=15, learning_rate=.05, n_estimators=100, n_jobs=4, random_state=0, verbose=-1).fit(Xtr, Ytr[:, j]).predict(Xte) for j in range(Ytr.shape[1])])
            fv = np.unique(fov[tr]); hold = np.random.default_rng(0).choice(fv, max(1, round(.15 * len(fv))), replace=False); v = np.isin(fov[tr], hold)
            if not resid_on_ridge: return lgbm_es(Xtr, Ytr, Xte, v=v)
            m = RidgeCV(alphas=f3.ALPHAS).fit(Xtr, Ytr); return lgbm_es(Xtr, Ytr, Xte, m.predict(Xtr), m.predict(Xte), v=v)
        fp.wants_rows = True
        s, fo = m3.block_cv(X, Y, D["xy"][rows], fp, f3.RF_UM, blocks=D["blocks"][rows], split=D["split"], covars=dep)
        out[r] = (s, fo, [])
    return out


import os
FOV = os.environ.get("SF_AUDIT_FOV") == "1"
for u in sys.argv[1:]:
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); p = json.loads(Path(f"manifests/panels/{u}.json").read_text()); g = a.var_names.get_indexer(p["genes"])
    D = f3.prepare(a.X[:, g], a.obs[["x_um", "y_um"]].values, a.obs.fov.values); rep = dict(unit=u)
    if FOV:
        ell, feats = float(os.environ.get("SF_AUDIT_ELL", 25)), os.environ.get("SF_AUDIT_FEATS", "T+ctx"); tag = "" if (ell, feats) == (25., "T+ctx") else f"_{feats}_ell{int(ell)}"
        vs = (("lgbm_es_fov", False, False), ("lgbm_resid_fov", True, False)) + ((("lgbm_locked", False, True),) if tag else ())
        for k, rr, lk in vs:
            t = time.time(); rep[k] = dict(unit=f3.unit_score(score_fov(D, rr, ell, feats, lk)), s=round(time.time() - t)); print(u, k + tag, rep[k], flush=True)
        Path(f"runs/p3/audit_f3_resid_fov_{u}{tag}.json").write_text(json.dumps(stamp(rep), indent=1, default=float)); continue
    for k, fp in (("lgbm_es", lgbm_es), ("lgbm_resid", resid)):
        t = time.time(); s = f3.score_floor(D, 25., 50, fit_predict=fp); rep[k] = dict(unit=f3.unit_score(s), s=round(time.time() - t))
        print(u, k, rep[k], flush=True)
    for m in ("ridge", "lightgbm"):
        q = Path(f"runs/p3/f3_{u}_{m}.json")
        if q.exists(): rep[m + "_grid"] = {c: v["unit"] for c, v in json.loads(q.read_text())["grid"].items()}
    Path(f"runs/p3/audit_f3_resid_{u}.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
