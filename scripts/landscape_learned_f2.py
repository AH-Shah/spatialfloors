"""F2 learned comparators on the landscape perturbations (PLAN §14 (e)): the model at its P3-best config, no re-tuning.
perturb "none" with seed s = a fresh model seed (embedding seed variance); otherwise perturbation seed s (stress.perturb),
model seed 0. Clustering seeds 0..4. -> runs/landscape/labels/f2_<unit>_<model>_<perturb>_s<seed>.npz (labels: seeds x cells)
Usage: python scripts/landscape_learned_f2.py <unit> cellcharter|novae <config> <perturb> <seed>
  (env cellcharter for cellcharter, MPS; env novae for novae zero_shot, CPU)"""
import sys, time
from pathlib import Path
import numpy as np
from spatialfloors.stress import perturb
from spatialfloors.stress.f2 import load
from spatialfloors.metrics.f2 import scored_regions

u, model, cfg, name, s = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]); d = load(u)
X, c, g = perturb.apply(name, d["X"], d["fov"], s if name != "none" else 0); xy = d["xy"][c]; mseed = s if name == "none" else 0
k = len(np.unique(scored_regions(d["region"][c], d["present"])[1])); t = time.time()
if model == "cellcharter":
    from spatialfloors.models import cellcharter
    nl, L = int(cfg.split("_")[0][2:]), int(cfg.split("_L")[1]); labs, _ = cellcharter.fit(X, xy, k, nl, L, accelerator="mps", seed=mseed, gmm_seeds=list(range(5))); lab = np.stack(labs)
elif model == "novae":
    import anndata as ad, novae
    assert cfg == "zero_shot"
    rev = {"human": ("prism-oncology/novae-human-0", "b8c0a5d7612bac6bc719ab57ed3cd16ad814728c"), "mouse": ("prism-oncology/novae-mouse-0", "3a23be90d5c6d87595ce3d5da522a6854f05564e")}
    import json
    genes = np.asarray(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"])[g]
    b = ad.AnnData(X.copy(), var=dict(index=genes)); b.var_names = genes; b.obsm["spatial"] = xy; novae.utils.spatial_neighbors(b, radius=80)
    mdl = novae.Novae.from_pretrained(*rev["mouse" if u.startswith("allen") else "human"][:1], revision=rev["mouse" if u.startswith("allen") else "human"][1])
    mdl.compute_representations(b, zero_shot=True, accelerator="cpu"); lab = b.obs[mdl.assign_domains(b, n_domains=k)].astype(str).values[None]
else: raise SystemExit(model)
out = Path("runs/landscape/labels"); out.mkdir(parents=True, exist_ok=True)
np.savez_compressed(out / f"f2_{u}_{model}_{name}_s{s}.npz", labels=lab); print(u, model, cfg, name, s, lab.shape, round(time.time() - t, 1), flush=True)
