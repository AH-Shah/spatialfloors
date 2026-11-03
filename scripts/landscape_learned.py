"""Learned comparators on the landscape perturbations, F1 / F3 / F5 (PLAN §14 (e)): each model at its P3-best config,
no re-tuning (F2: landscape_learned_f2.py). perturb "none" with seed s = a fresh model seed; otherwise perturbation seed s,
model seed 0. Inputs are built exactly as the matching landscape script builds its panel, so the outputs are scored there:
  f1 scvi             emb  -> runs/landscape/emb/f1_<unit>_scvi_<perturb>_s<seed>.npz       (landscape_f1.py learned)
  f5 scanvi/tangram   labels -> runs/landscape/labels/f5_<unit>_<model>_<perturb>_s<seed>.npz (landscape_f5.py learned)
  f3 gat              panel  -> runs/landscape/learned_f3/<unit>_gat_<perturb>_s<seed>.json    (landscape_f3.py learned)
Usage (pixi env scvi, MPS): python scripts/landscape_learned.py <fam> <unit> <model> <perturb> <seed>"""
import glob, json, sys, time
from pathlib import Path
import numpy as np
from spatialfloors.stress import perturb

fam, u, model, name, s = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]); ps, ms = (0, s) if name == "none" else (s, 0); t = time.time()
best = lambda f, key="unit": max(json.loads(Path(f).read_text())["grid"].items(), key=lambda kv: kv[1][key])[0]
if fam == "f1":
    import anndata as ad, scipy.sparse as sp, scvi, torch
    from spatialfloors.stress.f1 import load, Panel
    cfg = max(glob.glob(f"runs/p3/f1_{u}_scvi_nl*.json"), key=lambda f: next(iter(json.loads(Path(f).read_text())["scores"].values()))["score"])
    nl, lik = Path(cfg).stem.split("_scvi_")[1].split("_"); d0 = load(u)
    if name.startswith("genes"): P = Panel(d0, 0, genes=np.flatnonzero(np.isin(np.arange(len(d0["g"])), perturb.subsample(len(d0["g"]), .75, ps))), X=d0["X"])
    elif name == "none": P = Panel(d0, 0)
    else: X, c, _ = perturb.apply(name, d0["X"], d0["fov"], ps); P = Panel(d0, 0, cells=c if name.startswith(("sub", "fov")) else None, X=X)
    b = ad.AnnData(sp.csr_matrix(P.Xa[:, P.g]).astype(np.float32)); scvi.model.SCVI.setup_anndata(b); scvi.settings.seed = ms
    m = scvi.model.SCVI(b, n_latent=int(nl[2:]), gene_likelihood=lik); m.train(accelerator="mps", devices=1, enable_progress_bar=False)
    out = Path("runs/landscape/emb"); out.mkdir(parents=True, exist_ok=True); np.savez_compressed(out / f"f1_{u}_scvi_{name}_s{s}.npz", emb=m.get_latent_representation())
    torch.mps.empty_cache(); info = f"{nl}_{lik}"
elif fam == "f5":
    import torch
    from spatialfloors.floors import f5
    from spatialfloors.models import transfer
    from spatialfloors.stress.f5 import load, inputs
    d0 = load(u); Xq, genes, *_ = inputs(d0, name, ps); Xr, yr = f5.reference(genes, counts=True); info = best(f"runs/p3/f5_{u}_{model}.json", "heldout_auroc")
    if model == "scanvi": nl, lk = info.split("_"); lab = transfer.scanvi(Xr, yr, Xq, n_latent=int(nl[2:]), likelihood=lk, seed=ms)
    else: lab = transfer.tangram(Xr, yr, Xq, density=info, seed=ms)
    out = Path("runs/landscape/labels"); out.mkdir(parents=True, exist_ok=True); np.savez_compressed(out / f"f5_{u}_{model}_{name}_s{s}.npz", labels=np.asarray(lab).astype(str))
    torch.mps.empty_cache()
elif fam == "f3":
    import anndata as ad, scipy.sparse as sp, torch
    from spatialfloors.floors import f3
    from spatialfloors.models import gat
    from spatialfloors.provenance import stamp
    from spatialfloors.stress.f3 import panel, summary
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); g = a.var_names.get_indexer(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"])
    X0, XY, FOV = sp.csr_matrix(a.X[:, g]), a.obs[["x_um", "y_um"]].values.astype(float), a.obs.fov.astype(str).values
    X, c, _ = perturb.apply(name, X0, FOV, ps); ids = f3.prepare(X0, XY, FOV, seed=0)["ids"][c]; D = f3.prepare(X, XY[c], FOV[c], ids=ids, seed=0)  # receivers fixed (landscape_f3)
    info = best(f"runs/p3/f3_{u}_gat.json"); L = int(info.split("_L")[1])
    per = panel(gat.score(D, 25., L, 128, 50, device="mps", conv="gat", seed=ms, full=True)[0])
    out = Path("runs/landscape/learned_f3"); out.mkdir(parents=True, exist_ok=True)
    (out / f"{u}_gat_{name}_s{s}.json").write_text(json.dumps(stamp(dict(unit=u, model="gat", config=info, perturb=name, perturb_seed=ps, model_seed=ms, **summary(per), per_receiver=per)), default=float))
    torch.mps.empty_cache()
else: raise SystemExit(fam)
print(fam, u, model, info, name, s, round(time.time() - t, 1), flush=True)
