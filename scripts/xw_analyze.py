"""Whole-slide arm analysis (results/WHOLESLIDE_PLAN.md; exploratory, no decision weight). Run once after every xw_*
job is done. Units: the six GSE1 confirmation slides as whole slides (gse1_<sid>w), one stratum (CosMx WTx), df 5.
Per pair (F1 scVI, F2 CellCharter, F3 LightGBM, F5 scANVI; locked configurations, registered metrics):
  D_whole = learned - floor on the whole slide; one-sided 95% bounds with the family margin (as tier 1);
  learned superior: lower bound > 0, Holm over the 4 pairs;
  scale effect Delta = D_whole - D_crop per slide (D_crop = the registered tier-1 unit value, P5_corrected.json),
  mean with a two-sided 90% CI; floor-lead change = 1 - mean D_whole / mean D_crop.
Secondary (no test weight): F3 LightGBM - v2 floor (no rings); F5 at scANVI's call rate (as scripts/x_f5_matched.py).
Usage (env floors): python scripts/xw_analyze.py runs/xw runs/p4/P5_corrected.json -> runs/xw/X_wholeslide.{json,md}"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, scipy.stats as st
from spatialfloors import stats
from spatialfloors.data.genes import filter_reason
from spatialfloors.data.units import index
from spatialfloors.floors import f5 as F5
from spatialfloors.metrics.f2 import scored_regions, ari_bal_saved
from spatialfloors.metrics.f5 import CLASSES, heldout_auroc
from spatialfloors.provenance import commit

R, T1 = Path(sys.argv[1]), json.loads(Path(sys.argv[2]).read_text())["pairs"]; SIDS = ("110", "120", "210", "221", "222", "242")
L = json.loads(Path("manifests/locked_v3.json").read_text())["locked"]; j = lambda p: json.loads(Path(p).read_text()) if Path(p).exists() else None
MK = json.loads(Path("manifests/panels/index.json").read_text())["f5_heldout_markers"]
rows = {}
for s in SIDS:
    u = f"gse1_{s}w"; r = {}
    f1 = {m: j(R / f"f1_{u}_{m}.json") for m in ("expr_embed_pearson_none_pc50", "scvi_nl50_nb")}
    if all(f1.values()): r["F1|scvi"] = (f1["scvi_nl50_nb"]["scores"]["nl50_nb"]["score"], f1["expr_embed_pearson_none_pc50"]["scores"]["pearson_none_pc50"]["score"])
    if (R / f"labels/f2_{u}_floor.npz").exists() and (R / f"labels/f2_{u}_cellcharter.npz").exists():
        reg = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad", backed="r").obs.region.astype(str).values; m, rl = scored_regions(reg, index()[u]["present"])
        r["F2|cellcharter"] = (ari_bal_saved(R / f"labels/f2_{u}_cellcharter.npz", "nl10_L3", m, rl), ari_bal_saved(R / f"labels/f2_{u}_floor.npz", "ell400_lam0.8", m, rl))
    f3 = {k: j(R / p) for k, p in (("r", f"f3_{u}_ridge.json"), ("l", f"f3_{u}_lightgbm.json"), ("v2", f"v2/f3_{u}_ridge.json"))}
    if f3["r"] and f3["l"]: r["F3|lightgbm"] = (f3["l"]["grid"]["leaves15_lr0.05"]["unit"], f3["r"]["grid"]["v3"]["unit"])
    if f3["v2"] and f3["l"]: r["F3 LightGBM - v2 floor"] = (f3["l"]["grid"]["leaves15_lr0.05"]["unit"], f3["v2"]["grid"]["T+ctx:ell25:pc10:a1000"]["unit"])
    f5 = {k: j(R / f"f5_{u}_{k}.json") for k in ("floors", "scanvi")}
    if all(f5.values()):
        r["F5|scanvi"] = (f5["scanvi"]["grid"]["nl30_nb"]["heldout_auroc"], f5["floors"]["grid"]["insitutype_spatial"]["heldout_auroc"])
        a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); keep = filter_reason(a.var_names, "human") == ""; X, vn, xy = sp.csr_matrix(a.X[:, keep]), list(map(str, a.var_names[keep])), a.obs[["x_um", "y_um"]].values.astype(float)
        LL = pd.read_csv(R / f"ist_{u}" / "logliks.csv"); from scipy.spatial import cKDTree
        ll = LL.values; s2 = np.sort(ll, 1); z0 = ll / max(float(np.median(s2[:, -1] - s2[:, -2])), 1e-6); nb = cKDTree(xy).query(xy, 16)[1][:, 1:]
        P = np.exp(z0 - z0.max(1, keepdims=True)); P /= P.sum(1, keepdims=True)
        for _ in range(5): z = z0 + P[nb].mean(1); z -= z.max(1, keepdims=True); P = np.exp(z); P /= P.sum(1, keepdims=True)
        names = np.asarray(LL.columns); cls = F5.coarse(names); lab = F5.coarse(np.load(R / f"labels/f5_{u}_scanvi.npz", allow_pickle=True)["nl30_nb"].astype(str))
        assert (names[P.argmax(1)] == np.load(R / f"labels/f5_{u}_floors.npz", allow_pickle=True)["insitutype_spatial"].astype(str)).all()
        fm = []
        for c in CLASSES:
            pr = np.full(len(lab), "x", dtype=object); pr[np.argsort(-P[:, cls == c].sum(1), kind="stable")[:int((lab == c).sum())]] = c
            fm.append(heldout_auroc(pr, X, vn, MK, classes=(c,))[1][c]["auroc"])
        r["F5 scANVI - floor at scANVI's call rate"] = (r["F5|scanvi"][0], float(np.mean(fm))); del a, X, P, nb
    rows[s] = r; print(s, {k: [round(x, 4) for x in v] for k, v in r.items()}, flush=True)
fam = {"F1|scvi": "F1", "F2|cellcharter": "F2", "F3|lightgbm": "F3", "F5|scanvi": "F5", "F3 LightGBM - v2 floor": "F3", "F5 scANVI - floor at scANVI's call rate": "F5"}
out = dict(commit=commit(), slides=list(SIDS), per_slide=rows, pairs={})
for k, f in fam.items():
    us = [s for s in SIDS if k in rows[s]]
    if len(us) < 2: continue
    D = np.array([rows[s][k][0] - rows[s][k][1] for s in us]); t = stats.family_test(D, None, margin=L[f]["margin"]); o = dict(n=len(us), est=t["est"], lower=t["lower"], upper=t["upper"], margin=L[f]["margin"],
        noninferior=t["noninferior"], p_superior=t["p_superior"], floor_mean=float(np.mean([rows[s][k][1] for s in us])), learned_mean=float(np.mean([rows[s][k][0] for s in us])), D_whole=dict(zip(us, D.tolist())))
    if k in T1:
        Dc = np.array([T1[k]["per_unit"][f"gse1_{s}"] for s in us]); dl = D - Dc; ci = st.t.ppf(.95, len(dl) - 1) * dl.std(ddof=1) / np.sqrt(len(dl))
        o.update(D_crop=dict(zip(us, Dc.tolist())), D_crop_mean=float(Dc.mean()), delta_mean=float(dl.mean()), delta_ci90=[float(dl.mean() - ci), float(dl.mean() + ci)],
                 floor_lead_change=float(1 - D.mean() / Dc.mean()) if Dc.mean() < 0 else None)
    out["pairs"][k] = o
tp = [k for k in out["pairs"] if "|" in k]; ph = stats.holm([out["pairs"][k]["p_superior"] for k in tp]) if tp else []
for k, h in zip(tp, ph): out["pairs"][k]["learned_superior_holm"] = bool(h)
md = ["# X_wholeslide: the six GSE1 confirmation slides as whole slides (exploratory; results/WHOLESLIDE_PLAN.md)", "", f"commit {commit()}", "",
      "| Pair | Slides | D whole | One-sided 95% bounds | δ | Non-inferior | Learned superior (Holm) | D crop (same slides) | Δ = whole − crop (90% CI) | Floor lead change |", "|---|---|---|---|---|---|---|---|---|---|"]
for k, o in out["pairs"].items():
    md.append(f"| {k} | {o['n']} | {o['est']:+.4f} | {o['lower']:+.4f}, {o['upper']:+.4f} | {o['margin']:.4f} | {o['noninferior'] if '|' in k else '–'} | {o.get('learned_superior_holm', '–')} | "
              + (f"{o['D_crop_mean']:+.4f} | {o['delta_mean']:+.4f} ({o['delta_ci90'][0]:+.4f}, {o['delta_ci90'][1]:+.4f}) | " + (f"{o['floor_lead_change']:+.0%}" if o["floor_lead_change"] is not None else "–") + " |" if "D_crop" in o else "– | – | – |"))
md += ["", "Per slide (learned, floor): " + "; ".join(f"{s}: " + ", ".join(f"{k.split('|')[0] if '|' in k else k} {v[0]:.3f}/{v[1]:.3f}" for k, v in rows[s].items()) for s in SIDS)]
(R / "X_wholeslide.json").write_text(json.dumps(out, indent=1, default=float)); (R / "X_wholeslide.md").write_text("\n".join(md) + "\n"); print("\n".join(md))
