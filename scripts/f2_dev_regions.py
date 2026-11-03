"""Region truth available for F2 development (no model, no metric): per-region cell counts on the whole GSE1 232 slide
and on every Allen Zhuang animal 4 section (development animal, PLAN §11). Usage: python scripts/f2_dev_regions.py"""
import anndata as ad, pandas as pd

R = "/Users/Hemato/project_3/data/raw"
a = ad.read_h5ad(f"{R}/GSE1/cosmx_wt_232.h5ad", backed="r"); cols = ["typ"]
print("232 whole slide", a.n_obs, "cells; region-like columns", cols)
for c in cols[:4]: print(c, a.obs[c].astype(str).str.extract(r"(REF|TVA|CRC)", expand=False).value_counts(dropna=False).head(8).to_dict())
b = ad.read_h5ad(f"{R}/allen_zhuang/WB_MERFISH_animal4_sagittal.h5ad", backed="r"); o = b.obs
sec = next(c for c in ("brain_section_label", "section", "slice") if c in o.columns); reg = "major_brain_region"
print("animal4", b.n_obs, "cells;", sec, reg)
t = pd.crosstab(o[sec].astype(str), o[reg].astype(str)); t["n"] = t.sum(1); print(t.loc[:, t.sum(0) > 2000].to_string())
