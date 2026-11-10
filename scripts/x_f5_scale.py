"""Diagnostic: why scANVI's F5 score drops on whole slides while the floor
does not (runs/xw/X_wholeslide.md). Locked config nl30_nb, seed 0, same inputs and metric as scripts/p3_f5.py on the
whole-slide cache gse1_<sid>w. The registered whole-slide run (scVI/scANVI epochs from the scvi-tools heuristic on
ref + all cells) is the reference. Two variants, each scored on every cell of the slide:
  A_epochs  all cells, but crop-equivalent epochs: scVI = heuristic(n_ref + 50,000), scANVI = min(10, max(2, scVI/3))
  B_ratio   train on ref + 50,000 random query cells (rng 0; crop-size ratio, default epochs), predict every cell
A recovers -> the epoch heuristic; B recovers but A not -> the query:reference ratio.
Usage (env scvi): python scripts/x_f5_scale.py 221 -> runs/xw/diag/X_f5_scale_<sid>.json"""

import json, os, sys, time
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp, scvi, torch
from scvi.model._utils import get_max_epochs_heuristic
from spatialfloors.floors import f5
from spatialfloors.metrics.f5 import unit_inputs, score_classes
from spatialfloors.provenance import stamp

sid = sys.argv[1]
u = f"gse1_{sid}w"
OUT = Path("runs/xw/diag")
OUT.mkdir(parents=True, exist_ok=True)
MK = f5.for_unit(u, json.loads(Path("manifests/panels/index.json").read_text())["f5_heldout_markers"])
held = {x for v in MK.values() for x in v}
a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad")
p = json.loads(Path(f"manifests/panels/{u}.json").read_text())
pg = set(f5.pelka_genes())
genes = [x for x in p["genes"] if x not in held and x in pg]
Xq = sp.csr_matrix(a[:, genes].X)
Xf, vf, gate, cls = unit_inputs(u, a)
Xr, yr = f5.reference(genes, counts=True)
reg = json.loads(Path(f"runs/xw/f5_{u}_scanvi.json").read_text())["grid"]["nl30_nb"]["heldout_auroc"]
flo = json.loads(Path(f"runs/xw/f5_{u}_floors.json").read_text())["grid"]["insitutype_spatial"]["heldout_auroc"]


def fit_predict(Xtr, Xpred, max_epochs=None, scanvi_epochs=None):
    A = ad.AnnData(sp.vstack([sp.csr_matrix(Xr), Xtr]).tocsr().astype(np.float32))
    A.obs["modality"] = ["ref"] * Xr.shape[0] + ["query"] * Xtr.shape[0]
    A.obs["label"] = pd.Categorical(list(map(str, yr)) + ["Unknown"] * Xtr.shape[0])
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(A, batch_key="modality", labels_key="label")
    m = scvi.model.SCVI(A, n_latent=30, gene_likelihood="nb")
    m.train(max_epochs=max_epochs, accelerator="mps", devices=1, enable_progress_bar=False)
    s = scvi.model.SCANVI.from_scvi_model(m, unlabeled_category="Unknown")
    s.train(max_epochs=scanvi_epochs, accelerator="mps", devices=1, enable_progress_bar=False)
    Q = ad.AnnData(Xpred.astype(np.float32))
    Q.obs["modality"] = pd.Categorical(["query"] * Q.n_obs, categories=A.obs.modality.astype("category").cat.categories)
    Q.obs["label"] = pd.Categorical(["Unknown"] * Q.n_obs, categories=A.obs.label.cat.categories)
    return np.asarray(s.predict(Q)), int(m.history["elbo_train"].shape[0]), int(s.history["elbo_train"].shape[0])


rep = dict(unit=u, n_cells=int(a.n_obs), n_ref=int(Xr.shape[0]), registered_whole=reg, floor_whole=flo, variants={})
E = get_max_epochs_heuristic(Xr.shape[0] + 50000)
E2 = int(min(10, max(2, round(E / 3))))
idx = np.sort(np.random.default_rng(0).choice(Xq.shape[0], 50000, replace=False))
VAR = os.environ.get("SF_VARIANTS", "B_ratio,A_epochs").split(",")  # variants to run
for name, kw in (("B_ratio", dict(Xtr=Xq[idx])), ("A_epochs", dict(Xtr=Xq, max_epochs=E, scanvi_epochs=E2))):
    if name not in VAR:
        continue
    t = time.time()
    lab, e1, e2 = fit_predict(Xpred=Xq, **kw)
    r = score_classes(f5.coarse(lab), Xf, vf, MK, gate, cls)
    rep["variants"][name] = dict(
        heldout_auroc=r["heldout_auroc"],
        per_class=r["per_class"],
        pred_share=r["pred_share"],
        scvi_epochs=e1,
        scanvi_epochs=e2,
        s=round(time.time() - t, 1),
    )
    print(
        u,
        name,
        round(r["heldout_auroc"], 4),
        "registered",
        round(reg, 4),
        "floor",
        round(flo, 4),
        "epochs",
        e1,
        e2,
        round(time.time() - t),
        flush=True,
    )
    torch.mps.empty_cache()
    (OUT / f"X_f5_scale_{sid}.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
print("X_F5_SCALE_DONE", flush=True)
