"""Audit 2026-11-02 (pre-lock review): on the whole slide (gse1_231w, 298k cells) scVI beats the F1 floor (+.0095) while
it loses on the crop (-.010); audit_crop_scale.py shows the floor's PCA, not scVI, changes with population size. Does
the floor recover with more PCs at whole-slide scale? Floor = expr_embed (Pearson residuals, z) at n_pc {100, 200, 300},
whole-slide count split seed 0, panel of gse1_231w, kNN over all cells (as G7 'D whole').
Usage (env floors, OMP 4): python scripts/audit_f1_scale.py -> runs/p3/audit_f1_scale.json"""
import json, sys, time
from pathlib import Path
import anndata as ad, numpy as np
from spatialfloors.data.genes import filter_reason
from spatialfloors.floors.simple import expr_embed
from spatialfloors.metrics.f1 import count_split, knn_poisson_score
from spatialfloors.provenance import stamp

w = ad.read_h5ad("runs/p2/cache/gse1_231w.h5ad"); keep = np.flatnonzero(filter_reason(w.var_names, "human") == "")
Xa, Xb = count_split(w.X, seed=0); Xa, Xb = Xa[:, keep], Xb[:, keep]; del w
g = np.asarray(json.loads(Path("manifests/panels/gse1_231w.json").read_text())["genes"]); vn = ad.read_h5ad("runs/p2/cache/gse1_231w.h5ad", backed="r").var_names[keep]
gi = vn.get_indexer(g); assert (gi >= 0).all(); out = dict(scvi_nl30_nb=float(knn_poisson_score(np.load("runs/p3/f1_gse1_231w_scvi_nl30_nb_emb.npz")["nl30_nb"], Xa, Xb, gi)))
print(out, flush=True)
for nm, scl, n in [x.split("_") for x in (sys.argv[1:] or ["pearson_z_100", "pearson_z_200", "pearson_z_300"])]:
    t = time.time(); n = int(n); E = expr_embed(Xa[:, gi], n, nm, scl); out[f"floor_{nm}_{scl}_pc{n}"] = float(knn_poisson_score(E, Xa, Xb, gi)); print(nm, scl, n, out[f"floor_{nm}_{scl}_pc{n}"], round(time.time() - t), flush=True); del E
Path(f"runs/p3/audit_f1_scale{'_' + '-'.join(sys.argv[1:]) if sys.argv[1:] else ''}.json").write_text(json.dumps(stamp(out), indent=1))
