"""F3 niche -> receiver state on one crop (shared design floors.f3.prepare; quadrant CV); every config of the model's grid,
or only SF_CONFIG. SF_F3_FEATS=v3: the frozen receiver ids (runs/p2/f3_ids) and the v3 feature groups.
  ridge       v3 floor ('v3', alpha 1000; 'v3_a<alpha>'), fixed-alpha floors '<features>:ell<um>:pc<n>:a<alpha>', or the
              early RidgeCV grid ell {10, 25 um} x context PCs {30, 50}                                         (env floors)
  lightgbm    leaves {15, 63} x lr {.05, .1}, 100 rounds, one model per target PC, 4 threads, on the floor's features
  lgbm_resid  ridge + LightGBM on its residual (secondary)
  gat         GATConv x layers {1, 2}, hidden 128 (MPS; descriptive)                                        (env scvi)
Usage: python scripts/p3_f3.py <unit> <model> -> $SF_OUT/f3_<unit>_<model>.json (default runs/p3)"""
import json, sys, time
from pathlib import Path
import anndata as ad, numpy as np
from spatialfloors.floors import f3
from spatialfloors.provenance import stamp

import os
u, model = sys.argv[1], sys.argv[2]; OUT = Path(os.environ.get("SF_OUT", "runs/p3")); OUT.mkdir(parents=True, exist_ok=True); ONLY = os.environ.get("SF_CONFIG")  # P4: one locked config
a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad")
p = json.loads(Path(f"manifests/panels/{u}.json").read_text()); g = a.var_names.get_indexer(p["genes"]); assert (g >= 0).all()
fz = Path(f"runs/p2/f3_ids/{u}.npy"); fz = fz if os.environ.get("SF_F3_FEATS") == "v3" and fz.exists() else None  # v3: the frozen receiver ids the v3 floor uses
D = f3.prepare(a.X[:, g], a.obs[["x_um", "y_um"]].values, a.obs.fov.values, ids=None if fz is None else np.load(fz))
rep = dict(unit=u, model=model, K=D["K"], receivers={int(r): int((D["ids"] == r).sum()) for r in f3.receivers(D)}, grid={})
pack = lambda s: dict(unit=f3.unit_score(s), per_receiver={int(k): dict(score=v[0], folds=v[1], pc1=v[2][0]) for k, v in s.items()})


def lgbm(leaves, lr):
    def fp(Xtr, Ytr, Xte):
        from lightgbm import LGBMRegressor
        return np.column_stack([LGBMRegressor(num_leaves=leaves, learning_rate=lr, n_estimators=100, n_jobs=4, random_state=0, verbose=-1)
                                .fit(Xtr, Ytr[:, j]).predict(Xte) for j in range(Ytr.shape[1])])
    return fp


if model == "gat":
    from spatialfloors.models import gat
    cfgs = {f"{c}_L{L}": dict(conv=c, layers=L) for c in ("gat",) for L in (1, 2)}
    score = lambda kw: gat.score(D, 25., kw["layers"], 128, 50, device="mps", conv=kw["conv"], ckpt=OUT / f"f3_{u}_gat_{kw['conv']}_L{kw['layers']}.partial.json")[0]  # resume per receiver
else:
    # PLAN §18 (2026-11-02): floor configs '<features>:ell<um>:pc<n>:a<alpha>' (fixed-alpha ridge; human T+ctx:ell25:pc10:a1000,
    # mouse ctx:ell100:pc50:a1000; ':' because '|' separates job-file fields); LightGBM / lgbm_resid read the floor's features from SF_F3_FEATS ('<features>|ell<um>|pc<n>',
    # default the P3 'T+ctx|ell25|pc50'); lgbm_resid = ridge (alpha SF_F3_ALPHA) + LightGBM on its residual (secondary, not a test pair)
    import re
    V3 = os.environ.get("SF_F3_FEATS") == "v3"  # PLAN §19: floors.f3.score_v3 (ridge = the v3 floor: lock + same-lineage rings, alpha 1000; LightGBM = every v3 group)
    fs = re.match(r"^(T\+ctx|ctx|T)[|:]ell(\d+)[|:]pc(\d+)$", "T+ctx|ell25|pc50" if V3 else os.environ.get("SF_F3_FEATS", "T+ctx|ell25|pc50")); F_ = dict(features=fs[1], ell_um=float(fs[2]), n_pc=int(fs[3]))
    if V3: F_ = dict(features="v3", species="mouse" if u.startswith("allen") else "human")
    m_ = re.match(r"^(T\+ctx|ctx|T)[|:]ell(\d+)[|:]pc(\d+)[|:]a(.+)$", ONLY or "")
    ridge_a = lambda a: (lambda Xtr, Ytr, Xte: __import__("sklearn.linear_model", fromlist=["Ridge"]).Ridge(alpha=a).fit(Xtr, Ytr).predict(Xte))
    if model == "ridge" and m_: cfgs = {ONLY: dict(features=m_[1], ell_um=float(m_[2]), n_pc=int(m_[3]), fit_predict=ridge_a(float(m_[4])))}
    elif model == "ridge": cfgs = {f"ell{int(e)}_pc{n}": dict(ell_um=e, n_pc=n) for e in (10., 25.) for n in (30, 50)}
    elif model == "lightgbm": cfgs = {f"leaves{l}_lr{lr}": dict(F_, fit_predict=lgbm(l, lr)) for l in (15, 63) for lr in (.05, .1)}
    else:  # lgbm_resid: up to 1000 rounds, early stopping (30) on ~15% of whole training FOVs (f3.groups); fixed 100 rounds
        # on the residual overfit (278N .064 vs ridge .085), early stopping on spatial groups ties the ridge (audit_f3_resid.py)
        from lightgbm import LGBMRegressor, early_stopping
        from sklearn.linear_model import Ridge
        from spatialfloors.metrics import f3 as m3
        A = float(os.environ.get("SF_F3_ALPHA", 1000)); cfgs = {"resid_leaves15_lr0.05_es": {}}

        def score(kw):
            nb, w = f3.kernel(D, F_["ell_um"]); T = f3.id_mix(D, nb, w); out = {}
            for r in f3.receivers(D):
                rows, Y, dep = f3.target(D, r); ctx = f3.context_fn(D, r, rows, nb, w, F_["n_pc"]); grp = f3.groups(D["fov"][rows], D["xy"][rows])
                X = {"T": lambda q, tr: T[rows], "ctx": ctx, "T+ctx": lambda q, tr: np.c_[T[rows], ctx(q, tr)]}[F_["features"]]

                def fp(Xtr, Ytr, Xte, tr, te, q):
                    fv = np.unique(grp[tr]); v = np.isin(grp[tr], np.random.default_rng(0).choice(fv, max(1, round(.15 * len(fv))), replace=False))
                    rg = Ridge(alpha=A).fit(Xtr, Ytr); R, P = Ytr - rg.predict(Xtr), rg.predict(Xte)
                    for j in range(Ytr.shape[1]):
                        m = LGBMRegressor(num_leaves=15, learning_rate=.05, n_estimators=1000, n_jobs=4, random_state=0, verbose=-1)
                        m.fit(Xtr[~v], R[~v, j], eval_set=[(Xtr[v], R[v, j])], callbacks=[early_stopping(30, verbose=False)]); P[:, j] += m.predict(Xte)
                    return P
                fp.wants_rows = True
                s, fo = m3.block_cv(X, Y, D["xy"][rows], fp, f3.RF_UM, blocks=D["blocks"][rows], split=D["split"], covars=dep); out[r] = (s, fo, [np.nan])
            return out
    rep["features"] = F_ if model != "ridge" or V3 else None
    if V3 and model == "ridge": cfgs = {ONLY: dict(fit_predict=None, alpha=float(ONLY[4:]))} if (ONLY or "").startswith("v3_a") else {"v3": dict(fit_predict=None)}  # tier 2: v3_a<alpha>
    if V3 and model in ("ridge", "lightgbm"): score = lambda kw: f3.score_v3(D, F_["species"], kw["fit_predict"], **{k: kw[k] for k in ("alpha",) if k in kw})
    elif model != "lgbm_resid": score = lambda kw: f3.score_floor(D, **kw)
cfgs = {k: v for k, v in cfgs.items() if ONLY in (None, k)}; assert cfgs, f"no config {ONLY}"
for k, kw in cfgs.items():
    t = time.time(); rep["grid"][k] = dict(pack(score(kw)), s=round(time.time() - t, 1)); print(u, model, k, round(rep["grid"][k]["unit"], 4), rep["grid"][k]["s"], flush=True)
(OUT / f"f3_{u}_{model}.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
