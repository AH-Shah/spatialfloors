"""F3 development stress panel (PLAN §13): is a niche -> state score real coupling, or an artefact of the readout?
Development crops only; same shared design and quadrant CV as P3 (floors.f3), predictions kept (block_cv detail full).
Per model / variant and receiver: top-10 PC mean Spearman (primary), PC1, top-3, median per-gene Spearman (predictions
mapped back to z-scored genes; N2 row MP3: the 10-PC mean dilutes one-axis signals), PC-space R2 (top 10).
Variants: floor T+ctx at its P3-best config; T only (type mix: MP3 found the floor signal is mostly the mix); ctx only;
features shuffled within FOV x quadrant (null: must fall to ~0; N2 rows NX1 / NX2: structure-only learners and one-seed
gains misled); LightGBM at its P3-best config (+ shuffled); GAT best (models.gat, --gat, MPS).
Unit summary: mean over receivers, worst receiver, receivers above the shuffled null + .01.
Usage (env scvi): python scripts/f3_stress.py <unit> [--gat] -> runs/p3/f3_stress_<unit>.json"""
import json, sys, time
from pathlib import Path
import anndata as ad, numpy as np
from scipy.stats import spearmanr
from spatialfloors.floors import f3
from spatialfloors.provenance import stamp

u = sys.argv[1]; a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); p = json.loads(Path(f"manifests/panels/{u}.json").read_text())
g = a.var_names.get_indexer(p["genes"]); D = f3.prepare(a.X[:, g], a.obs[["x_um", "y_um"]].values, a.obs.fov.values); log = lambda *s: print(*s, flush=True)
best = lambda m: max(json.loads(Path(f"runs/p3/f3_{u}_{m}.json").read_text())["grid"].items(), key=lambda kv: kv[1]["unit"])[0]
fc = best("ridge"); ell, npc = float(fc.split("_")[0][3:]), int(fc.split("_pc")[1])


def lgbm(leaves, lr):
    def fp(Xtr, Ytr, Xte):
        from lightgbm import LGBMRegressor
        return np.column_stack([LGBMRegressor(num_leaves=leaves, learning_rate=lr, n_estimators=100, n_jobs=4, random_state=0, verbose=-1)
                                .fit(Xtr, Ytr[:, j]).predict(Xte) for j in range(Ytr.shape[1])])
    return fp


def panel(out):
    """out: {receiver: (score, folds, per_pc, detail)} -> per-receiver stress metrics."""
    res = {}
    for r, (s, _, _, det) in out.items():
        pc1, top3, gene, r2 = [], [], [], []
        for f in det.get("folds", []):
            P, T = f["pred"], f["true"]; rho = [spearmanr(P[:, j], T[:, j])[0] for j in range(T.shape[1])]
            pc1.append(rho[0]); top3.append(np.nanmean(rho[:3])); r2.append(1 - ((P - T) ** 2).sum() / ((T - T.mean(0)) ** 2).sum())
            G = P @ f["comp"]; gene.append(np.nanmedian([spearmanr(G[:, j], f["Yz"][:, j])[0] for j in range(G.shape[1]) if f["Yz"][:, j].std() > 0]))
        res[int(r)] = dict(top10=float(s), pc1=float(np.nanmean(pc1)), top3=float(np.nanmean(top3)), gene_median=float(np.nanmean(gene)), r2=float(np.nanmean(r2)))
    return res


runs = {"floor": dict(ell_um=ell, n_pc=npc), "floor_T": dict(ell_um=ell, n_pc=npc, features="T"), "floor_ctx": dict(ell_um=ell, n_pc=npc, features="ctx"),
        "floor_shuffled": dict(ell_um=ell, n_pc=npc, shuffle_seed=0)}
lc = best("lightgbm"); lv, lr = int(lc.split("_")[0][6:]), float(lc.split("_lr")[1])
runs["lightgbm"] = dict(ell_um=25., n_pc=50, fit_predict=lgbm(lv, lr)); runs["lightgbm_shuffled"] = dict(ell_um=25., n_pc=50, fit_predict=lgbm(lv, lr), shuffle_seed=0)
rep = dict(unit=u, floor_config=fc, lightgbm_config=lc, variants={})
for name, kw in runs.items():
    t = time.time(); rep["variants"][name] = dict(per_receiver=panel(f3.score_floor(D, full=True, **kw)), s=round(time.time() - t, 1))
    log(name, {k: round(float(np.nanmean([v[k] for v in rep["variants"][name]["per_receiver"].values()])), 4) for k in ("top10", "pc1", "top3", "gene_median", "r2")}, rep["variants"][name]["s"])
if "--gat" in sys.argv:
    from spatialfloors.models import gat
    gc = best("gat"); L = int(gc.split("_L")[1]); t = time.time()
    rep["variants"]["gat"] = dict(config=gc, per_receiver=panel(gat.score(D, 25., L, 128, 50, device="mps", conv="gat", full=True)[0]), s=round(time.time() - t, 1))
null = rep["variants"]["floor_shuffled"]["per_receiver"]
for name, v in rep["variants"].items():
    pr = v["per_receiver"]; v["summary"] = {k: float(np.nanmean([x[k] for x in pr.values()])) for k in ("top10", "pc1", "top3", "gene_median", "r2")}
    v["summary"].update(worst_top10=float(np.nanmin([x["top10"] for x in pr.values()])), above_null=float(np.mean([pr[r]["top10"] > null[r]["top10"] + .01 for r in pr if r in null])))
Path(f"runs/p3/f3_stress_{u}.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
for name, v in rep["variants"].items(): log("SUMMARY", name, {k: round(x, 4) for k, x in v["summary"].items()})
