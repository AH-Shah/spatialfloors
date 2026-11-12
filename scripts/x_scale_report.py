"""Exploratory scale curves: report. Per GSE1 confirmation slide and cell count n, D = learned - floor for F1 scVI,
F2 CellCharter, F3 LightGBM and F5 scANVI (locked configurations, registered metrics, scored as scripts/xw_analyze.py)
on the windows gse1_<sid>s<n> (runs/xs, scripts/x_scale_prepare.py) and on the whole slide gse1_<sid>w (runs/xw).
Per pair: mean D over slides at each n; per slide the OLS slope of D on log10(cells) over its available points; mean
slope over slides with a two-sided 95% t interval (df = slides - 1). Missing outputs are skipped and counted.
Usage (env floors, round-1 runs root as working directory): python scripts/x_scale_report.py runs/xs runs/xw
  -> <repo>/results/reports/X_scale.{json,md}"""

import json, sys
from pathlib import Path
import numpy as np, scipy.stats as st, anndata as ad
from spatialfloors.data.units import index
from spatialfloors.metrics.f2 import scored_regions, ari_bal_saved
from spatialfloors.provenance import commit

P = Path(__file__).resolve().parents[1]
XS, XW = Path(sys.argv[1]), Path(sys.argv[2])
SIDS = ("110", "120", "210", "221", "222", "242")
SIZES = (12500, 25000, 50000, 100000, 200000, 400000)
j = lambda p: json.loads(Path(p).read_text()) if Path(p).exists() else None


def scores(R, u):
    r = {}
    f1 = {m: j(R / f"f1_{u}_{m}.json") for m in ("expr_embed_pearson_none_pc50", "scvi_nl50_nb")}
    if all(f1.values()):
        r["F1|scvi"] = (f1["scvi_nl50_nb"]["scores"]["nl50_nb"]["score"], f1["expr_embed_pearson_none_pc50"]["scores"]["pearson_none_pc50"]["score"])
    if (R / f"labels/f2_{u}_floor.npz").exists() and (R / f"labels/f2_{u}_cellcharter.npz").exists():
        a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad", backed="r")
        m, rl = scored_regions(a.obs.region.astype(str).values, index()[u]["present"])
        a.file.close()
        r["F2|cellcharter"] = (ari_bal_saved(R / f"labels/f2_{u}_cellcharter.npz", "nl10_L3", m, rl), ari_bal_saved(R / f"labels/f2_{u}_floor.npz", "ell400_lam0.8", m, rl))
    f3 = {k: j(R / f"f3_{u}_{k}.json") for k in ("ridge", "lightgbm")}
    if all(f3.values()):
        r["F3|lightgbm"] = (f3["lightgbm"]["grid"]["leaves15_lr0.05"]["unit"], f3["ridge"]["grid"]["v3"]["unit"])
    f5 = {k: j(R / f"f5_{u}_{k}.json") for k in ("floors", "scanvi")}
    if all(f5.values()):
        r["F5|scanvi"] = (f5["scanvi"]["grid"]["nl30_nb"]["heldout_auroc"], f5["floors"]["grid"]["insitutype_spatial"]["heldout_auroc"])
    return r


rows = []  # pair, sid, n, learned, floor
for s in SIDS:
    pts = [(XS, f"gse1_{s}s{n}", n) for n in SIZES] + [(XW, f"gse1_{s}w", None)]
    for R, u, n in pts:
        if not Path(f"runs/p2/cache/{u}.h5ad").exists():
            continue
        cells = n or ad.read_h5ad(f"runs/p2/cache/{u}.h5ad", backed="r").n_obs
        for k, (le, fl) in scores(R, u).items():
            rows.append(dict(pair=k, sid=s, n=int(cells), point="whole" if n is None else str(n), learned=le, floor=fl, D=le - fl))
out = dict(commit=commit(), rows=rows, pairs={})
md = ["# X_scale: floor vs learned across cell counts (exploratory)", "", f"commit {out['commit']}", "",
      "Windows by the registered crop rule at each cell count, gene panel and count split of the whole slide; D = learned - floor.", ""]
for k in ("F1|scvi", "F2|cellcharter", "F3|lightgbm", "F5|scanvi"):
    R_ = [r for r in rows if r["pair"] == k]
    if not R_:
        continue
    by = {p: [r["D"] for r in R_ if r["point"] == p] for p in [str(n) for n in SIZES] + ["whole"]}
    sl = []
    for s in SIDS:
        q = [r for r in R_ if r["sid"] == s]
        if len(q) >= 3:
            sl.append(float(np.polyfit(np.log10([r["n"] for r in q]), [r["D"] for r in q], 1)[0]))
    t = dict(mean_D={p: (float(np.mean(v)), len(v)) for p, v in by.items() if v}, slopes=sl)
    if len(sl) >= 2:
        m, se = float(np.mean(sl)), float(np.std(sl, ddof=1) / np.sqrt(len(sl)))
        h = float(st.t.ppf(0.975, len(sl) - 1) * se)
        t.update(slope_mean=m, slope_ci=(m - h, m + h))
    out["pairs"][k] = t
    md += [f"## {k}", "", "| cells | slides | mean D |", "|---|---|---|"]
    md += [f"| {p} | {c} | {v:+.4f} |" for p, (v, c) in t["mean_D"].items()]
    if "slope_mean" in t:
        md += ["", f"Slope of D per log10(cells): mean {t['slope_mean']:+.4f}, 95% CI {t['slope_ci'][0]:+.4f}, {t['slope_ci'][1]:+.4f} ({len(sl)} slides)"]
    md.append("")
(P / "results/reports/X_scale.json").write_text(json.dumps(out, indent=1, default=float))
(P / "results/reports/X_scale.md").write_text("\n".join(md))
print("\n".join(md))
