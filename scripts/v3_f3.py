"""F3 floor v3 on one unit (PLAN §19; development units for the choice, validation units once with the chosen floor).
Shared F3 design (frozen receiver ids runs/p2/f3_ids, bleed-safe FOV-centred targets, quadrant CV + edge buffer, block_cv).
Per receiver: T (type mix only, human: decomposition), lock (v2 floor: ridge alpha 1000 on the locked features), v3
(floors.f3 v3_parts + GroupRidge: lock + same-lineage rings 50-150 / 150-250 um on the receiver's own genes + sender-
lineage context, one penalty per group by inner grouped CV), v3_null (both rings permuted within FOV x quadrant).
Usage (env floors, OMP_NUM_THREADS=3): [SF_OUT=runs/v3] python scripts/v3_f3.py <unit> -> <out>/f3_<unit>.json"""
import json, os, sys, time
from pathlib import Path
import anndata as ad, numpy as np
from sklearn.linear_model import Ridge
from spatialfloors.floors import f3
from spatialfloors.metrics import f3 as m3
from spatialfloors.provenance import stamp

u = sys.argv[1]; t0 = time.time(); log = lambda *s: print(f"[{time.time() - t0:6.0f}s]", *s, flush=True)
OUT = Path(os.environ.get("SF_OUT", "runs/v3")); OUT.mkdir(parents=True, exist_ok=True); mouse = u.startswith("allen")
a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); g = a.var_names.get_indexer(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"]); assert (g >= 0).all()
D = f3.prepare(a.X[:, g], a.obs[["x_um", "y_um"]].values, a.obs.fov.values, ids=np.load(f"runs/p2/f3_ids/{u}.npy")); del a
LOCK = dict(feat="ctx", ell=100., npc=50) if mouse else dict(feat="T+ctx", ell=25., npc=10)
K = f3.kernel(D, LOCK["ell"]); T = f3.id_mix(D, *K); R = {k: f3.ring(D, *v) for k, v in f3.RINGS_V3.items()}
log("design", D["K"], "ids; ring coverage", {k: round(float((v[1].sum(1) > 0).mean()), 3) for k, v in R.items()})
ridge = lambda Xtr, Ytr, Xte: Ridge(alpha=1000.).fit(Xtr, Ytr).predict(Xte)


def cv(Xf, fp, rows, Y, dep):
    return m3.block_cv(Xf, Y, D["xy"][rows], fp, f3.RF_UM, blocks=D["blocks"][rows], split=D["split"], covars=dep)[0]


def fov_quadrant_perm(rows, seed=0):
    grp = np.char.add(D["fov"][rows].astype(str), D["blocks"][rows].astype(str)); perm = np.arange(len(rows)); rng = np.random.default_rng(seed)
    for v in np.unique(grp):
        i = np.flatnonzero(grp == v); perm[i] = rng.permutation(i)
    return perm


rep = dict(unit=u, lock=LOCK, rings=f3.RINGS_V3, alphas=f3.V3_ALPHAS, mults=[m if m is not None else "drop" for m in f3.V3_MULTS], receivers={})
for r in f3.receivers(D):
    rows, Y, dep = f3.target(D, r); Ec = f3.svd_fn(D, np.flatnonzero(~D["mask"][r]), LOCK["npc"]); Es = f3.svd_fn(D, np.flatnonzero(D["mask"][r]), f3.V3_NPC_OWN)
    out = {}
    if not mouse: out["T"] = cv(lambda q, tr: T[rows], ridge, rows, Y, dep)
    out["lock"] = cv(f3.lock_part(rows, LOCK["feat"], K, T, Ec), ridge, rows, Y, dep)
    F3 = f3.v3_parts(D, r, rows, LOCK["feat"], K, T, R, Ec, Es); GR = f3.GroupRidge(F3, D["blocks"][rows]); out["v3"] = cv(F3, GR, rows, Y, dep)
    FN = f3.v3_parts(D, r, rows, LOCK["feat"], K, T, R, Ec, Es, shuffle=fov_quadrant_perm(rows)); GN = f3.GroupRidge(FN, D["blocks"][rows]); out["v3_null"] = cv(FN, GN, rows, Y, dep)
    rep["receivers"][int(r)] = dict(n=int(len(rows)), lineage=int(D["lin"][rows[0]]), **out, choices=GR.choices, choices_null=GN.choices)
    log("receiver", r, len(rows), {k: round(v, 4) for k, v in out.items()}, [(c["alpha"], tuple(c["mult"].values())[1:]) for c in GR.choices])
recs = list(rep["receivers"].values()); names = [k for k in ("T", "lock", "v3", "v3_null") if k in recs[0]]
rep["unit_mean"] = {k: float(np.nanmean([x[k] for x in recs])) for k in names}
rep["n_v3_gt_lock"] = int(sum(x["v3"] > x["lock"] for x in recs)); rep["n_receivers"] = len(recs)
log("unit", {k: round(v, 4) for k, v in rep["unit_mean"].items()}, "v3>lock", rep["n_v3_gt_lock"], "/", len(recs))
(OUT / f"f3_{u}.json").write_text(json.dumps(stamp(rep), indent=1, default=float)); print("V3F3_DONE", u, flush=True)
