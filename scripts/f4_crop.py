"""F4 crop rule (PLAN §4 amendment, as GSE1 110): slides too large for Cellina's whole-slide memory (> ~150k cells; Cellina
peaked 9.7 GB on crc_232's 131k) use their ~50k crop (manifests/crops/gse1_<sid>.csv.gz) for Cellina, CPA and L6 alike.
Writes the raw slide's rows (all obs / var untouched) to runs/f4/raw/cosmx_wt_<sid>.h5ad, which f4_run / f4_l6 / f4_score read.
Usage (env floors): python scripts/f4_crop.py <sid>"""
import os, sys
from pathlib import Path
import anndata as ad, pandas as pd

sid = sys.argv[1]; out = Path(os.environ.get("SF_F4_RUN", "runs/f4")) / "raw"; out.mkdir(parents=True, exist_ok=True)
ids = pd.read_csv(f"manifests/crops/gse1_{sid}.csv.gz").cell_id.astype(str).values
a = ad.read_h5ad(f"/Users/Hemato/project_3/data/raw/GSE1/cosmx_wt_{sid}.h5ad", backed="r"); b = a[ids].to_memory(); a.file.close()
b.write_h5ad(out / f"cosmx_wt_{sid}.h5ad"); print(sid, b.shape, b.obs.typ.value_counts().to_dict())
