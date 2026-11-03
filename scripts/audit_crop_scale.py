"""Audit 2026-11-02 (pre-lock review): is the G7 crop-vs-whole difference on gse1_231 a training-data effect or a
scoring-population effect? Whole-slide outputs (trained / clustered on all 298k cells) are restricted to the crop's 50k
cells and scored there with the crop's own scoring rules, then compared with the crop-trained D (G7 'D crop') and the
whole-slide D (G7 'D whole').
F1: scVI nl30_nb and floor pearson_z_pc100 embeddings of gse1_231w; whole-slide count split (seed 0, the split the
embeddings were fitted on) restricted to crop cells; kNN among crop cells; panel genes of gse1_231w.
F2: whole-slide labels (floor ell200_lam0.8, CellCharter nl10_L3; 10 seeds) restricted to crop cells; balanced ARI vs
the crop's scored regions (as p5_analyze.py).
Usage (env floors): python scripts/audit_crop_scale.py -> runs/p3/audit_crop_scale.json"""
import json
from pathlib import Path
import anndata as ad, numpy as np
from spatialfloors.data.genes import filter_reason
from spatialfloors.metrics.f1 import count_split, knn_poisson_score
from spatialfloors.metrics.f2 import scored_regions
from spatialfloors.stress.f2 import Panel
from spatialfloors.provenance import stamp

out = {}
c = ad.read_h5ad("runs/p2/cache/gse1_231.h5ad", backed="r"); cid, creg = c.obs_names.copy(), c.obs.region.astype(str).values; c.file.close()
w = ad.read_h5ad("runs/p2/cache/gse1_231w.h5ad"); pos = w.obs_names.get_indexer(cid); assert (pos >= 0).all(), "crop cells missing from whole slide"
keep = np.flatnonzero(filter_reason(w.var_names, "human") == ""); Xa, Xb = count_split(w.X, seed=0); Xa, Xb = Xa[:, keep], Xb[:, keep]
g = w.var_names[keep].get_indexer(json.loads(Path("manifests/panels/gse1_231w.json").read_text())["genes"]); assert (g >= 0).all(); del w
E = {"scvi": np.load("runs/p3/f1_gse1_231w_scvi_nl30_nb_emb.npz")["nl30_nb"], "floor": np.load("runs/p3/f1_gse1_231w_expr_embed_pearson_z_pc100_emb.npz")["pearson_z_pc100"]}
Xac, Xbc = Xa[pos], Xb[pos]
f1 = {m: dict(on_crop_cells=float(knn_poisson_score(e[pos], Xac, Xbc, g))) for m, e in E.items()}
f1["D_whole_trained_on_crop_cells"] = f1["scvi"]["on_crop_cells"] - f1["floor"]["on_crop_cells"]; out["F1"] = f1; print("F1", f1, flush=True)
idx = json.loads(Path("manifests/crops/index.json").read_text()); mask, reg = scored_regions(creg, idx["gse1_231"]["present"])
P = Panel.__new__(Panel); P.y = reg; f2 = {}
for m, key in (("floor", "ell200_lam0.8"), ("cellcharter", "nl10_L3")):
    L = np.asarray(np.load(f"runs/p3/labels/f2_gse1_231w_{m}.npz", allow_pickle=True)[key]).astype(str)[:, pos]
    f2[m] = dict(on_crop_cells=float(np.mean([P.ari_bal(l[mask]) for l in L])))
f2["D_whole_clustered_on_crop_cells"] = f2["cellcharter"]["on_crop_cells"] - f2["floor"]["on_crop_cells"]; out["F2"] = f2; print("F2", f2, flush=True)
Path("runs/p3/audit_crop_scale.json").write_text(json.dumps(stamp(out), indent=1))
