"""F1 fit + score of one embedding on one crop. Every fit on X_a (count split seed 0 of the full crop); library sizes
over
the unit's filtered genes (PROTOCOL §2); Nicheformer inputs per platform (models.nicheformer.inputs).
  expr_embed <norm>_<scale>_pc<n> (floor; locked pearson_none_pc50) | expr_pca pc<n> (env floors)
  scvi nl<n>_{nb,zinb} (env scvi, MPS) | nicheformer layers (final and penultimate, env nicheformer, MPS)
Usage: python scripts/f1_fit.py <unit> <model> <config> [--out DIR] -> DIR/f1_<unit>_<model>_<config>.json + _emb.npz
(DIR defaults to $SF_OUT, then runs/p3). The primary is re-scored from the saved embeddings by p3_f1_rescore.py."""

import os
import json, sys, time
from pathlib import Path
import anndata as ad, numpy as np, scipy.sparse as sp
from spatialfloors.data.genes import filter_reason
from spatialfloors.metrics.f1 import count_split, knn_poisson_score, chance_level
from spatialfloors.provenance import stamp

u, model, cfg = sys.argv[1:4]
OUT = Path(sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else os.environ.get("SF_OUT", "runs/p3"))
OUT.mkdir(parents=True, exist_ok=True)
tag = f"f1_{u}_{model}_{cfg}"
a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad")
species = "mouse" if u.startswith("allen") else "human"
_sp = Path(
    f"runs/p2/cache/{u}_split.npz"
)  # whole-slide arm (scripts/xw_prepare.py): the full-slide count split on the lean cache's columns
Xa_all, Xb_all = (
    (lambda A: (A, sp.csr_matrix(a.X - A)))(sp.load_npz(_sp).tocsr()) if _sp.exists() else count_split(a.X, seed=0)
)
keep = np.flatnonzero(filter_reason(a.var_names, species) == "")
Xa, Xb, vn = Xa_all[:, keep], Xb_all[:, keep], a.var_names[keep]
g = vn.get_indexer(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"])
assert (g >= 0).all()
t = time.time()
info = {}
if model == "expr_embed":  # F1 floor v2 (locked: pearson_z_pc100)
    from spatialfloors.floors.simple import expr_embed

    nm, scl, pc = cfg.split("_")
    E = {cfg: expr_embed(Xa[:, g], int(pc[2:]), nm, scl)}
elif model == "expr_pca":
    from spatialfloors.floors.simple import expr_pca

    E = {cfg: expr_pca(Xa[:, g], int(cfg.replace("pc", "")))}
elif model == "scvi":
    import scvi, torch

    nl, lik = cfg.split("_")
    scvi.settings.seed = 0
    b = ad.AnnData(sp.csr_matrix(Xa[:, g]), var=a.var.iloc[keep[g]].copy())
    scvi.model.SCVI.setup_anndata(b)
    m = scvi.model.SCVI(b, n_latent=int(nl.replace("nl", "")), gene_likelihood=lik)
    m.train(accelerator="mps", devices=1, enable_progress_bar=False)
    E = {cfg: m.get_latent_representation()}
    h = m.history["elbo_train"].values.astype(float)
    info = dict(epochs=len(h), elbo_last=float(h[-1, 0]), elbo_finite=bool(np.isfinite(h).all()))
    torch.mps.empty_cache()
elif model == "nicheformer":
    from spatialfloors.models.nicheformer import embed, inputs

    q = ad.AnnData(
        sp.csr_matrix(Xa_all), obs=a.obs[[]].copy(), var=a.var.copy()
    )  # all genes of X_a (the model's own vocabulary decides)
    r, _ = embed(
        q,
        0,
        layers=(-1, -2),
        pool="masked",
        save_emb=OUT / f"{tag}_emb.npz",
        grad_check=False,
        **inputs(u, a.var_names),
    )  # zero-shot: no backward step
    z = dict(np.load(OUT / f"{tag}_emb.npz"))
    E = {"final": z["-1_masked"], "penultimate": z["-2_masked"]}
    info = {
        k: r[k]
        for k in (
            "vocab_share_genes",
            "vocab_share_counts",
            "tokens_nonpad_median",
            "masked_cells_per_s",
            "emb_finite",
            "vocab_gate_passed",
        )
    }
else:
    raise SystemExit(f"unknown model {model}")
rep = dict(
    unit=u,
    model=model,
    config=cfg,
    n_cells=int(a.n_obs),
    n_panel=len(g),
    fit_s=round(time.time() - t, 1),
    info=info,
    scores={},
    L="filtered genes",
)
for k, v in E.items():
    rep["scores"][k] = dict(
        score=float(knn_poisson_score(v, Xa, Xb, g)),
        chance=chance_level(v, Xa, Xb, g, n_perm=1),
        finite=bool(np.isfinite(v).all()),
        dim=int(v.shape[1]),
    )
if model != "nicheformer":
    np.savez(OUT / f"{tag}_emb.npz", **E)
print(tag, json.dumps(rep["scores"]), rep["fit_s"], flush=True)
(OUT / f"{tag}.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
