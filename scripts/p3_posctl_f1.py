"""P3 positive control for F1 scVI (decision 2026-10-28, option b): planted types (tests/test_metrics._groups: 5 types,
200 genes, Dirichlet .3 profiles, depth 300; n 5,000), count split seed 0, scVI fit on X_a (n_latent 10, nb, CPU, 200
epochs). Pass: F1 score >= the Expr PCA floor's - .05 and >= permutation chance + .20. Nicheformer cannot read synthetic
genes: its control is the vocabulary gate (>= 80%) plus F1 >= chance + .05 on every dev crop (P3 grid, runs/p3/f1_*).
Usage (env scvi): python scripts/p3_posctl_f1.py -> runs/p3/posctl_f1.json"""
import json, sys
from pathlib import Path
import anndata as ad, numpy as np, scipy.sparse as sp
sys.path.insert(0, "tests")
from test_metrics import _groups
from spatialfloors.metrics.f1 import count_split, knn_poisson_score, chance_level
from spatialfloors.floors.simple import expr_pca
from spatialfloors.provenance import stamp
import scvi

X, g, _ = _groups(n=5000, depth=300, seed=0); a_, b_ = count_split(X, 0); genes = np.arange(X.shape[1])
b = ad.AnnData(sp.csr_matrix(a_, dtype=np.float32)); scvi.model.SCVI.setup_anndata(b); scvi.settings.seed = 0
m = scvi.model.SCVI(b, n_latent=10); m.train(max_epochs=200, accelerator="cpu", enable_progress_bar=False); Z = m.get_latent_representation()
E = expr_pca(a_, 30); rep = dict(scvi=float(knn_poisson_score(Z, a_, b_, genes)), floor=float(knn_poisson_score(E, a_, b_, genes)), chance=chance_level(Z, a_, b_, genes, n_perm=1))
rep["pass"] = rep["scvi"] >= rep["floor"] - .05 and rep["scvi"] >= rep["chance"] + .2
print(json.dumps(rep), flush=True); Path("runs/p3/posctl_f1.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
