"""F3 v3 design comparison (PLAN §19, development units; after the full penalty search lost to v2 on Allen, .113 vs .125):
designs with less selection noise, each with a null that permutes every v3-only group within FOV x quadrant.
  lock         v2 (ridge alpha 1000 on the locked features)
  fixed_self   lock + same-lineage ring 50-150 um (ridge alpha 1000, one penalty)
  fixed_sf     lock + rings 50-150 and 150-250 um (alpha 1000)
  onoff        GroupRidge at alpha 1000, each v3 group (T, self, far, typed) on (m 1) or off; inner leave-one-training-
               quadrant-out CV (16 combinations instead of 1024)
  *_null       the same with the v3-only groups permuted
Usage (env floors): python scripts/v3_f3_alt.py <unit> -> runs/v3/f3alt_<unit>.json"""
import json, sys, time
from pathlib import Path
import anndata as ad, numpy as np
from sklearn.linear_model import Ridge
from spatialfloors.floors import f3
from spatialfloors.metrics import f3 as m3
from spatialfloors.provenance import stamp

u = sys.argv[1]; t0 = time.time(); log = lambda *s: print(f"[{time.time() - t0:6.0f}s]", *s, flush=True); mouse = u.startswith("allen")
a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); g = a.var_names.get_indexer(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"]); assert (g >= 0).all()
D = f3.prepare(a.X[:, g], a.obs[["x_um", "y_um"]].values, a.obs.fov.values, ids=np.load(f"runs/p2/f3_ids/{u}.npy")); del a
L = f3.V3_LOCK["mouse" if mouse else "human"]; K = f3.kernel(D, L["ell"]); T = f3.id_mix(D, *K); R = {k: f3.ring(D, *v) for k, v in f3.RINGS_V3.items()}
cv = lambda Xf, fp, rows, Y, dep: m3.block_cv(Xf, Y, D["xy"][rows], fp, f3.RF_UM, blocks=D["blocks"][rows], split=D["split"], covars=dep)[0]
ridge = lambda Xtr, Ytr, Xte: Ridge(alpha=1000.).fit(Xtr, Ytr).predict(Xte)


def perm_fq(rows, seed=0):
    grp = np.char.add(D["fov"][rows].astype(str), D["blocks"][rows].astype(str)); p = np.arange(len(rows)); rng = np.random.default_rng(seed)
    for v in np.unique(grp):
        i = np.flatnonzero(grp == v); p[i] = rng.permutation(i)
    return p


def fixed(rows, Ec, Es, rings, perm=None):
    lk = f3.lock_part(rows, L["feat"], K, T, Ec); p = np.arange(len(rows)) if perm is None else perm
    return lambda q, tr: np.hstack([lk(q, tr)] + [f3.bs.kernel_mean(*R[k], Es(q), rows)[p] for k in rings])


rep = dict(unit=u, receivers={})
for r in f3.receivers(D):
    rows, Y, dep = f3.target(D, r); Ec = f3.svd_fn(D, np.flatnonzero(~D["mask"][r]), L["npc"]); Es = f3.svd_fn(D, np.flatnonzero(D["mask"][r]), f3.V3_NPC_OWN); pm = perm_fq(rows); o = {}
    o["lock"] = cv(f3.lock_part(rows, L["feat"], K, T, Ec), ridge, rows, Y, dep)
    for nm, rg in (("fixed_self", ("self",)), ("fixed_sf", ("self", "far"))):
        o[nm] = cv(fixed(rows, Ec, Es, rg), ridge, rows, Y, dep); o[nm + "_null"] = cv(fixed(rows, Ec, Es, rg, pm), ridge, rows, Y, dep)
    F = f3.v3_parts(D, r, rows, L["feat"], K, T, R, Ec, Es); G = f3.GroupRidge(F, D["blocks"][rows], alphas=(1000.,), mults=(1., None)); o["onoff"] = cv(F, G, rows, Y, dep)
    FN = f3.v3_parts(D, r, rows, L["feat"], K, T, R, Ec, Es, shuffle=pm, shuffle_typed=True); GN = f3.GroupRidge(FN, D["blocks"][rows], alphas=(1000.,), mults=(1., None)); o["onoff_null"] = cv(FN, GN, rows, Y, dep)
    rep["receivers"][int(r)] = dict(n=int(len(rows)), **o, choices=G.choices); log("receiver", r, {k: round(v, 4) for k, v in o.items()})
recs = list(rep["receivers"].values()); names = [k for k in recs[0] if k not in ("n", "choices")]
rep["unit_mean"] = {k: float(np.nanmean([x[k] for x in recs])) for k in names}; log("unit", {k: round(v, 4) for k, v in rep["unit_mean"].items()})
(Path("runs/v3") / f"f3alt_{u}.json").write_text(json.dumps(stamp(rep), indent=1, default=float)); print("V3F3ALT_DONE", u, flush=True)
