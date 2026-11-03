"""Audit 2026-11-02 (pre-lock review, development units only): the F3 ridge floor picks alpha with RidgeCV's efficient
leave-one-out, which leaks spatial autocorrelation (neighbouring cells share niche and noise), so it may under-regularise;
the landscape's best Allen point uses a fixed alpha 1000 at ell 100 (.128 vs locked .091). Variants, same block_cv:
ridge alpha chosen by grouped CV over whole training FOVs (5 groups) instead of LOO, at ell {25, 100} um x features
{T+ctx, ctx}, 50 context PCs. Usage (env floors): python scripts/audit_f3_floor.py <unit ...> -> runs/p3/audit_f3_floor_<unit>.json"""
import json, sys, time
from pathlib import Path
import anndata as ad, numpy as np
from spatialfloors.floors import f3
from spatialfloors.floors.f3 import groups
from spatialfloors.metrics import f3 as m3
from spatialfloors.provenance import stamp


def ridge_fov(Xtr, Ytr, Xte, groups):
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import GroupKFold
    g = np.unique(groups, return_inverse=True)[1]; cv = GroupKFold(min(5, g.max() + 1)); best, ba = -np.inf, None
    for a in np.logspace(-2, 6, 17):  # grid extended to 1e6: grouped CV picked the old 1e4 edge on gse1_232
        sc = []
        for i, j in cv.split(Xtr, groups=g):
            p = Ridge(alpha=a).fit(Xtr[i], Ytr[i]).predict(Xtr[j]); sc.append(-((p - Ytr[j]) ** 2).mean())
        if np.mean(sc) > best: best, ba = np.mean(sc), a
    return Ridge(alpha=ba).fit(Xtr, Ytr).predict(Xte), ba


def score(D, ell, feats):
    nb, w = f3.kernel(D, ell); T = f3.id_mix(D, nb, w); out, alphas = {}, []
    for r in f3.receivers(D):
        rows, Y, dep = f3.target(D, r); ctx = f3.context_fn(D, r, rows, nb, w, 50); fov = groups(D["fov"][rows], D["xy"][rows])
        X = (lambda q, tr: np.c_[T[rows], ctx(q, tr)]) if feats == "T+ctx" else ctx

        def fp(Xtr, Ytr, Xte, tr, te, q):
            p, a = ridge_fov(Xtr, Ytr, Xte, fov[tr]); alphas.append(float(a)); return p
        fp.wants_rows = True
        s, fo = m3.block_cv(X, Y, D["xy"][rows], fp, f3.RF_UM, blocks=D["blocks"][rows], split=D["split"], covars=dep); out[r] = (s, fo, [])
    return f3.unit_score(out), float(np.median(alphas))


for u in sys.argv[1:]:
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); p = json.loads(Path(f"manifests/panels/{u}.json").read_text()); g = a.var_names.get_indexer(p["genes"])
    D = f3.prepare(a.X[:, g], a.obs[["x_um", "y_um"]].values, a.obs.fov.values); rep = dict(unit=u)
    for ell in (25., 100.):
        for feats in ("T+ctx", "ctx"):
            t = time.time(); s, am = score(D, ell, feats); k = f"{feats}|ell{int(ell)}|fovcv"; rep[k] = dict(unit=s, median_alpha=am, s=round(time.time() - t)); print(u, k, rep[k], flush=True)
    Path(f"runs/p3/audit_f3_floor_{u}.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
