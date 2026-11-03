"""F3 NCEM in the authors' setup (models.ncem.grid: the authors' ED-NCEM and linear interactions NCEM at radius 13 / 19 /
46 um, and our tuning on their training: type / type+ctx labels at higher capacity), with the engineering floor of each
(ridge on exactly its inputs, models.ncem.input_ridge) on the shared F3 design (receiver ids frozen by f3_ids_dump.py).
Development units; validation units only with SF_CONFIG = the selected configs. Env external/ncem_env, PYTHONPATH=src, CPU. Usage: external/ncem_env/bin/python
scripts/p3_f3_ncem_auth.py <unit> -> <SF_OUT or runs/p3v3>/f3_<unit>_ncemA.json (SF_CONFIG: comma list of configs)"""
import json, os, sys, time
from pathlib import Path
import anndata as ad, numpy as np
from spatialfloors.floors import f3
from spatialfloors.models import ncem
from spatialfloors.provenance import stamp

u = sys.argv[1]; OUT = Path(os.environ.get("SF_OUT", "runs/p3v3")); OUT.mkdir(parents=True, exist_ok=True); fn = OUT / f"f3_{u}_ncemA.json"
a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); g = a.var_names.get_indexer(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"]); assert (g >= 0).all()
D = f3.prepare(a.X[:, g], a.obs[["x_um", "y_um"]].values, a.obs.fov.values, ids=np.load(f"runs/p2/f3_ids/{u}.npy"))
NPC = 50 if u.startswith("allen") else 10  # context PCs = the relocked F3 floor's (PLAN §18)
G = ncem.grid(); ONLY = os.environ.get("SF_CONFIG"); G = {k: v for k, v in G.items() if ONLY is None or k in ONLY.split(",")}; assert G
rep = json.loads(fn.read_text()) if fn.exists() else dict(unit=u, model="ncemA", K=D["K"], receivers={int(r): int((D["ids"] == r).sum()) for r in f3.receivers(D)}, grid={}, ridge={})
pack = lambda s: dict(unit=f3.unit_score(s), per_receiver={int(k): dict(score=v[0], folds=v[1], pc1=v[2][0]) for k, v in s.items()})
for R in sorted({kw["radius_um"] for kw in G.values()} - {None}):
    for mode in ("mean", "presence"):
        k = f"{mode}_r{R:g}"
        if k not in rep["ridge"]: t = time.time(); rep["ridge"][k] = dict(pack(ncem.input_ridge(D, R, mode)), s=round(time.time() - t, 1)); print(u, "ridge", k, round(rep["ridge"][k]["unit"], 4), flush=True)
for k, kw in G.items():
    if k in rep["grid"]: continue
    kw = dict(kw, n_pc=NPC) if kw.get("h0") == "type+ctx" else kw
    t = time.time(); s, ep = ncem.score(D, 25., **kw); p = pack(s); tr = list(ep.values())
    p.update(s=round(time.time() - t, 1), unit_nan0=float(np.mean([v["nan0"] for v in tr])), nan_frac=float(np.mean([v["nan_frac"] for v in tr])),
             epochs_median=float(np.median([e for v in tr for e in v["epochs"]])), val_drop_frac=float(np.mean([l[1] < l[0] for v in tr for l in v["val"] if l])), train={int(r): v for r, v in ep.items()})
    rep["grid"][k] = p; fn.write_text(json.dumps(stamp(rep), indent=1, default=float))
    print(u, "ncemA", k, round(p["unit"], 4), "nan0", round(p["unit_nan0"], 4), "nan", round(p["nan_frac"], 3), "ep", p["epochs_median"], p["s"], flush=True)
fn.write_text(json.dumps(stamp(rep), indent=1, default=float)); print("NCEMA_DONE", u)
