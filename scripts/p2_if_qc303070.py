"""P2 ground-truth QC for F5 on GSE303070 (model-free: IF metadata vs the authors' Manual_toplevel_pred labels; no model
or floor output is read). Rule fixed before running (2026-10-28): an IF channel is valid on a unit iff its AUROC for the
authors' class vs all other labelled cells >= .75: PanCK for Epi; CD45 for immune (B, Plasma, TCD4, TCD8, Tgd, TZBTB16,
NK, ILC, Macro, Mono, DC, Granulo, Mast); CD68 for Macro + Mono. The 3-class F5 gate is usable on a unit only if all
three channels are valid. Trigger: dev crops 148191A show IF ~ uncorrelated with RNA (runs/p2/gate_check_dev.json).
Usage (env floors): python scripts/p2_if_qc303070.py -> manifests/if_qc_gse303070.json"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from spatialfloors.data.units import index, RAW
from spatialfloors.provenance import stamp

IMM = {"B", "Plasma", "TCD4", "TCD8", "Tgd", "TZBTB16", "NK", "ILC", "Macro", "Mono", "DC", "Granulo", "Mast"}
M = pd.read_csv(RAW / "GSE303070/GSE303070_All_cells_metadata.csv", index_col=False, low_memory=False,
                usecols=["Barcode", "Manual_toplevel_pred", "Mean.PanCK", "Mean.CD45", "Mean.CD68"]).set_index("Barcode")
out = {}
for u, rec in index().items():
    if not u.startswith("gse303070_"): continue
    ids = pd.read_csv(f"manifests/crops/{u}.csv.gz").cell_id.astype(str).values; m = M.loc[ids]; lab = m.Manual_toplevel_pred.astype(str).values
    ok = ~np.isin(lab, ["QC_fail", "nan", ""])
    auc = lambda pos, col: float(roc_auc_score(pos[ok], m[col].values[ok])) if 0 < pos[ok].sum() < ok.sum() and m[col].notna().all() else None
    r = dict(role=rec.get("role"), PanCK_epi=auc(lab == "Epi", "Mean.PanCK"), CD45_immune=auc(np.isin(lab, list(IMM)), "Mean.CD45"),
             CD68_macro_mono=auc(np.isin(lab, ["Macro", "Mono"]), "Mean.CD68"))
    r["gate_usable"] = all(v is not None and v >= .75 for k, v in r.items() if k != "role"); out[u] = r
    print(u, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
Path("manifests/if_qc_gse303070.json").write_text(json.dumps(stamp(out), indent=1))
