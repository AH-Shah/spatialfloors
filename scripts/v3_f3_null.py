"""F3 v3 null of every v3-only group (PLAN §19 acceptance): per receiver, lock (v2) and v3 with BOTH the same-lineage
rings and the sender-typed context permuted within FOV x quadrant (scripts/v3_f3.py's v3_null permutes the rings only,
so it keeps the typed context's real signal). Same design and GroupRidge as v3; a v3_null_all near lock means the
penalty search itself adds nothing. Usage (env floors): python scripts/v3_f3_null.py <unit> -> runs/v3/f3null_<unit>.json"""
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


def perm_fq(rows, seed=0):
    grp = np.char.add(D["fov"][rows].astype(str), D["blocks"][rows].astype(str)); p = np.arange(len(rows)); rng = np.random.default_rng(seed)
    for v in np.unique(grp):
        i = np.flatnonzero(grp == v); p[i] = rng.permutation(i)
    return p


rep = dict(unit=u, receivers={})
for r in f3.receivers(D):
    rows, Y, dep = f3.target(D, r); Ec = f3.svd_fn(D, np.flatnonzero(~D["mask"][r]), L["npc"]); Es = f3.svd_fn(D, np.flatnonzero(D["mask"][r]), f3.V3_NPC_OWN)
    lock = cv(f3.lock_part(rows, L["feat"], K, T, Ec), lambda Xtr, Ytr, Xte: Ridge(alpha=1000.).fit(Xtr, Ytr).predict(Xte), rows, Y, dep)
    FN = f3.v3_parts(D, r, rows, L["feat"], K, T, R, Ec, Es, shuffle=perm_fq(rows), shuffle_typed=True); GN = f3.GroupRidge(FN, D["blocks"][rows])
    rep["receivers"][int(r)] = dict(n=int(len(rows)), lock=lock, v3_null_all=cv(FN, GN, rows, Y, dep), choices=GN.choices)
    log("receiver", r, {k: round(v, 4) for k, v in rep["receivers"][int(r)].items() if k in ("lock", "v3_null_all")})
recs = list(rep["receivers"].values()); rep["unit_mean"] = {k: float(np.nanmean([x[k] for x in recs])) for k in ("lock", "v3_null_all")}
log("unit", rep["unit_mean"]); (Path("runs/v3") / f"f3null_{u}.json").write_text(json.dumps(stamp(rep), indent=1, default=float)); print("V3F3NULL_DONE", u, flush=True)
