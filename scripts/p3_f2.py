"""F2 domains on one crop with region labels (k = number of present regions, metrics.f2.scored_regions); every config of
the model's grid, or only SF_CONFIG. Partitions saved per config (seeds x cells) for the balanced-ARI primary.
  floor           full-covariance GMM (10 seeds) on [Expr PCA 50, neighbour mean], ell {100, 200, 400 um} x lambda {.8, 1};
                  'v3:<name>' = a multi-scale floor (floors.simple.v3_config)                                   (env floors)
  cellcharter     scVI n_latent {10, 30} x aggregation layers {1, 3}; scVI seeds SF_SCVI_SEEDS (default 0) x GMM seeds 0-9 (MPS)
  cellcharter_pca CellCharter on PCA input, 3 layers (published floor variant)
  novae           {zero_shot, fine_tune}; novae-human-0 / novae-mouse-0 (pinned), radius 80 um, assign_domains(k)    (CPU)
Usage: python scripts/p3_f2.py <unit> <model> -> $SF_OUT/f2_<unit>_<model>.json (default runs/p3)"""
import json, sys, time
from pathlib import Path
import anndata as ad, numpy as np
from sklearn.metrics import adjusted_rand_score
from spatialfloors.metrics.f2 import scored_regions, gmm_labels, SEEDS
from spatialfloors.data.units import index
from spatialfloors.provenance import stamp

import os
u, model = sys.argv[1], sys.argv[2]; OUT = Path(os.environ.get("SF_OUT", "runs/p3")); OUT.mkdir(parents=True, exist_ok=True); ONLY = os.environ.get("SF_CONFIG")  # P4: one locked config
want = lambda c: ONLY in (None, c); a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad")
p = json.loads(Path(f"manifests/panels/{u}.json").read_text()); g = a.var_names.get_indexer(p["genes"]); assert (g >= 0).all()
m, regl = scored_regions(a.obs.region.astype(str).values, index()[u]["present"]); k = len(np.unique(regl)); X = a.X[:, g]; xy = a.obs[["x_um", "y_um"]].values
rep = dict(unit=u, model=model, k=k, n_scored=int(m.sum()), grid={}); LAB = {}; ari = lambda lab: float(adjusted_rand_score(regl, np.asarray(lab)[m]))
if model == "floor" and ONLY and ONLY.startswith("v3:"):  # PLAN §19: F2 floor v3 configs (floors.simple.v3_config), GMM seeds as v2
    from spatialfloors.floors.simple import domain_features_v3, v3_config
    t = time.time(); F, _ = domain_features_v3(X, xy, **v3_config(ONLY[3:])); LAB[ONLY] = L_ = np.stack([gmm_labels(F, k, s) for s in SEEDS]); v = [ari(l) for l in L_]
    rep["grid"][ONLY] = dict(ari=float(np.mean(v)), ari_sd=float(np.std(v)), s=round(time.time() - t, 1))
elif model == "floor":
    from spatialfloors.floors.simple import domain_features
    for ell in (100., 200., 400.):  # 400: lock v2 (PLAN §18, G2 rule over dev + validation); floor v3 (PLAN §14, 2026-10-29): GMM, 50 PCs, grid around the landscape plateau centre ell 200 lam .8 (scripts/f2_choose.py);
        for lam in (.8, 1.):   # v2 (pc30, lam {.5, .8}) in runs/p3/f2_<u>_floor_v2.json, k-means v1 in runs/p3/f2_<u>_floor_kmeans.json
            if not want(f"ell{int(ell)}_lam{lam}"): continue
            t = time.time(); F, _ = domain_features(X, xy, a.obs.fov.values, 50, ell, lam=lam); LAB[f"ell{int(ell)}_lam{lam}"] = L_ = np.stack([gmm_labels(F, k, s) for s in SEEDS]); v = [ari(l) for l in L_]
            rep["grid"][f"ell{int(ell)}_lam{lam}"] = dict(ari=float(np.mean(v)), ari_sd=float(np.std(v)), s=round(time.time() - t, 1))
elif model == "cellcharter":
    from spatialfloors.models import cellcharter
    SSEEDS = [int(x) for x in os.environ.get("SF_SCVI_SEEDS", "0").split(",")]  # PLAN §18: scVI seeds averaged (P4 / pre-lock 0,1,2)
    for nl in (10, 30):
        for L in (1, 3):
            if not want(f"nl{nl}_L{L}"): continue
            t = time.time(); labs = []
            for ss in SSEEDS:
                l_, info = cellcharter.fit(X, xy, k, nl, L, accelerator="mps", seed=ss, gmm_seeds=list(SEEDS)); labs += list(l_)
            info = dict(info, scvi_seeds=SSEEDS); LAB[f"nl{nl}_L{L}"] = np.stack(labs); v = [ari(l) for l in labs]
            rep["grid"][f"nl{nl}_L{L}"] = dict(ari=float(np.mean(v)), ari_sd=float(np.std(v)), s=round(time.time() - t, 1), **info)
elif model == "cellcharter_pca":  # published floor variant (PLAN §3), reported, not a test pair
    from spatialfloors.models import cellcharter
    if want("pc50_L3"):
        t = time.time(); labs, info = cellcharter.fit(X, xy, k, n_layers=3, seed=0, gmm_seeds=list(SEEDS), rep="pca"); LAB["pc50_L3"] = np.stack(labs); v = [ari(l) for l in labs]
        rep["grid"]["pc50_L3"] = dict(ari=float(np.mean(v)), ari_sd=float(np.std(v)), s=round(time.time() - t, 1), **info)
elif model == "novae":
    import novae
    rev = {"human": ("prism-oncology/novae-human-0", "b8c0a5d7612bac6bc719ab57ed3cd16ad814728c"), "mouse": ("prism-oncology/novae-mouse-0", "3a23be90d5c6d87595ce3d5da522a6854f05564e")}
    sp = "mouse" if u.startswith("allen") else "human"; b = ad.AnnData(X.copy(), obs=a.obs[[]].copy(), var=a.var.iloc[g].copy()); b.obsm["spatial"] = xy
    novae.utils.spatial_neighbors(b, radius=80)
    for mode in ("zero_shot", "fine_tune"):
        if not want(mode): continue
        t = time.time(); mdl = novae.Novae.from_pretrained(rev[sp][0], revision=rev[sp][1])
        if mode == "fine_tune": mdl.fine_tune(b, accelerator="cpu")
        mdl.compute_representations(b, zero_shot=(mode == "zero_shot"), accelerator="cpu"); key = mdl.assign_domains(b, n_domains=k)
        Z = np.asarray(b.obsm["novae_latent"]); sv = np.linalg.svd(Z - Z.mean(0), compute_uv=False); pr = sv ** 2 / (sv ** 2).sum()
        lab = b.obs[key].astype(str).values; LAB[mode] = lab[None]; share = float(np.unique(lab, return_counts=True)[1].max() / len(lab))
        rep["grid"][mode] = dict(ari=ari(lab), effective_rank=float(np.exp(-(pr * np.log(pr + 1e-12)).sum())), largest_domain_share=share, s=round(time.time() - t, 1),
                                 collapse=bool(np.exp(-(pr * np.log(pr + 1e-12)).sum()) < 5 or share > .5), checkpoint=rev[sp])
for c, v in rep["grid"].items(): print(u, model, c, round(v["ari"], 3), v["s"], flush=True)
(OUT / "labels").mkdir(exist_ok=True); np.savez_compressed(OUT / "labels" / f"f2_{u}_{model}.npz", **LAB)  # all cells, seeds x cells (stress panel)
(OUT / f"f2_{u}_{model}.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
