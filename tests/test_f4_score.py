"""F4 scoring pin (audit 2026-10-30): the learned models' counterfactual pseudobulk, taken as per-cell CP10k of the saved cf
(count scale) and scored with floors.l6.score against L6's control / target pseudobulks, must reproduce Cellina's own
official pearson (eval_loo JSON) on crc_232 (P0 run). Without the per-cell normalisation Epithelial drops from .335 to .048."""
import json
from pathlib import Path
import numpy as np, pytest

P0 = Path("runs/p0_cellina/datasets/crc"); LFC = Path("/Users/Hemato/project_3/antig/N2/trials/perturb/opt/crc_232_lfc.npz")


@pytest.mark.skipif(not (P0.exists() and LFC.exists()), reason="P0 Cellina outputs not on disk")
@pytest.mark.parametrize("model", ["cellina_mps", "cpa"])
def test_cf_rescore_matches_official(model):
    import anndata as ad
    from spatialfloors.floors import l6
    Z = np.load(LFC); n = 0
    for h in ("Endothelial", "Epithelial", "Fibroblast", "Myeloid", "T_cell"):
        f = P0 / f"crc_232/{h}/{model}_counterfactual_x_CRC.h5ad"; j = P0 / f"correlations/crc_232_{model}-cf_{h}_CRC.json"
        if not (f.exists() and j.exists()): continue
        X = np.asarray(ad.read_h5ad(f).X, float); mf = (X / (X.sum(1, keepdims=True) + 1e-8) * 1e4).mean(0); t = f"{h}_CRC"
        assert abs(l6.score(Z[f"{t}|mc"], Z[f"{t}|mt"], mf, Z[f"{t}|bleed"])["pearson"] - json.loads(j.read_text())["pearson"]) < 1e-6; n += 1
    assert n >= 3
