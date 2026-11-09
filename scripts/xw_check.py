"""Whole-slide arm: lean paths (SF_LEAN=1, scripts/xw_prepare.py caches) against the registered paths on validation slide
GSE1 231 as a whole slide (gse1_231wl = the lean rebuild of gse1_231w, checked identical by xw_prepare.py --check).
References: runs/p3v2 (F1 / F2 locks, unchanged since v2), runs/g7v3 (F3 / F5 v3, G7 check), runs/p2/f3_ids/gse1_231w.npy.
Pass rules (fixed before the run): F1 |d score| < .001; F2 |d balanced ARI| < .02 (GMM on PCA computed by another
solver); F3 ids identical, |d ridge| < .001, |d LightGBM| < .002; F5 InSituType log-likelihoods identical (< 1e-9) and
|d AUROC| < 1e-9. Usage (env floors): python scripts/xw_check.py runs/xw/check -> <dir>/XW_CHECK.json; exit 1 on a fail."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
from spatialfloors.data.units import index
from spatialfloors.metrics.f2 import scored_regions, ari_bal_saved

C = Path(sys.argv[1]); j = lambda p: json.loads(Path(p).read_text()); out = {}
r, n = j("runs/p3v2/f1_gse1_231w_expr_embed_pearson_none_pc50.json"), j(C / "f1_gse1_231wl_expr_embed_pearson_none_pc50.json")
out["F1 floor score"] = (r["scores"]["pearson_none_pc50"]["score"], n["scores"]["pearson_none_pc50"]["score"], .001)
reg = ad.read_h5ad("runs/p2/cache/gse1_231wl.h5ad", backed="r").obs.region.astype(str).values; m, rl = scored_regions(reg, index()["gse1_231w"]["present"])
out["F2 floor balanced ARI"] = (ari_bal_saved("runs/p3v2/labels/f2_gse1_231w_floor.npz", "ell400_lam0.8", m, rl), ari_bal_saved(C / "labels/f2_gse1_231wl_floor.npz", "ell400_lam0.8", m, rl), .02)
out["F3 ids identical"] = (1.0, float(np.array_equal(np.load("runs/p2/f3_ids/gse1_231w.npy"), np.load("runs/p2/f3_ids/gse1_231wl.npy"))), 1e-12)
out["F3 ridge v3"] = (j("runs/g7v3/f3_gse1_231w_ridge.json")["grid"]["v3"]["unit"], j(C / "f3_gse1_231wl_ridge.json")["grid"]["v3"]["unit"], .001)
out["F3 LightGBM v3"] = (j("runs/g7v3/f3_gse1_231w_lightgbm.json")["grid"]["leaves15_lr0.05"]["unit"], j(C / "f3_gse1_231wl_lightgbm.json")["grid"]["leaves15_lr0.05"]["unit"], .002)
L0, L1 = pd.read_csv("runs/g7v3/ist_gse1_231w/logliks.csv").values, pd.read_csv(C / "ist_gse1_231wl/logliks.csv").values; fin = np.isfinite(L0)
out["F5 logliks max abs diff"] = (0.0, float(np.max(np.abs(L0[fin] - L1[fin]))) if (np.isfinite(L1) == fin).all() else np.inf, 1e-9)
out["F5 floor AUROC"] = (j("runs/g7v3/f5_gse1_231w_floors.json")["grid"]["insitutype_spatial"]["heldout_auroc"], j(C / "f5_gse1_231wl_floors.json")["grid"]["insitutype_spatial"]["heldout_auroc"], 1e-9)
res = {k: dict(registered=a, lean=b, diff=abs(a - b), tol=t, ok=bool(abs(a - b) < t)) for k, (a, b, t) in out.items()}
for k, v in res.items(): print(f"{k:28s} registered {v['registered']:.6f} lean {v['lean']:.6f} |d| {v['diff']:.2e} tol {v['tol']} {'PASS' if v['ok'] else 'FAIL'}")
(C / "XW_CHECK.json").write_text(json.dumps(res, indent=1)); sys.exit(0 if all(v["ok"] for v in res.values()) else 1)
