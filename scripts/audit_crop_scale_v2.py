"""Audit (lock v2 review): F2 G7 decomposition at the v2 configs, as scripts/audit_crop_scale.py (F2 part): whole-slide
labels of gse1_231w (floor ell400_lam0.8, 10 seeds; CellCharter nl10_L3, scVI seeds 0-2 x GMM seeds) restricted to the
crop's cells, balanced ARI vs the crop's scored regions; compared with the crop D and the whole-slide D (runs/p3v2/P5_gate.json).
Usage (env floors): python scripts/audit_crop_scale_v2.py -> runs/p3v2/audit_crop_scale_v2.json"""
import json
from pathlib import Path
import anndata as ad, numpy as np
from spatialfloors.metrics.f2 import scored_regions
from spatialfloors.stress.f2 import Panel
from spatialfloors.provenance import stamp

c = ad.read_h5ad("runs/p2/cache/gse1_231.h5ad", backed="r"); cid, creg = c.obs_names.copy(), c.obs.region.astype(str).values; c.file.close()
w = ad.read_h5ad("runs/p2/cache/gse1_231w.h5ad", backed="r"); pos = w.obs_names.get_indexer(cid); w.file.close(); assert (pos >= 0).all()
idx = json.loads(Path("manifests/crops/index.json").read_text()); mask, reg = scored_regions(creg, idx["gse1_231"]["present"])
P = Panel.__new__(Panel); P.y = reg; f2 = {}
for m, key in (("floor", "ell400_lam0.8"), ("cellcharter", "nl10_L3")):
    L = np.asarray(np.load(f"runs/p3v2/labels/f2_gse1_231w_{m}.npz", allow_pickle=True)[key]).astype(str)[:, pos]
    f2[m] = dict(on_crop_cells=float(np.mean([P.ari_bal(l[mask]) for l in L])))
f2["D_whole_clustered_on_crop_cells"] = f2["cellcharter"]["on_crop_cells"] - f2["floor"]["on_crop_cells"]
pu = json.loads(Path("runs/p3v2/P5_gate.json").read_text())["pairs"]["F2|cellcharter"]["per_unit"]; f2.update(D_crop=pu.get("gse1_231"), D_whole=pu.get("gse1_231w"))
print("F2", f2, flush=True); Path("runs/p3v2/audit_crop_scale_v2.json").write_text(json.dumps(stamp(dict(F2=f2)), indent=1))
