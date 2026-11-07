"""Secondary outputs of the confirmatory study (PROTOCOL §2 "Secondary outputs", §4 Novae collapse rule), no decision weight.
Written after the tier-1 report; reads the same sealed outputs as p5_analyze.py and adds no fitting. Per comparison: D = first -
second named, per crop -> unit means (stats.aggregate_units) with the tier-1 one-sided 95% bounds (stats.family_test), descriptive.
  F2  multi-scale secondary floor (runs/p4/f2v3) and CellCharter on PCA input, balanced ARI as the primary; Novae collapse rule
      (effective rank < 5 or > 50% of cells in one domain) per run.
  F3  v2 floor (runs/p4/v2: other lineages only, no rings) against both learned models; ED-NCEM (the authors' second NCEM).
  F5  InSituType without the spatial prior (= v2 floor); logistic regression at its development-best configuration (lock rule,
      manifests/locked.json); IF-gated macro-F1 (GSE1 2-class gate); marker-positive recall and class shares per model.
  Per-crop primaries of every model entering F2 (balanced ARI is computed from saved labels, so it is not in the run JSONs).
Usage (env floors): python scripts/p5_secondary.py <runs dir> <confirmation crop ...> -> <runs dir>/P5_secondary.{json,md}"""
import json, sys
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
from spatialfloors import stats
from spatialfloors.data.genes import filter_reason
from spatialfloors.floors import f5 as F5
from spatialfloors.metrics.f2 import ari_bal_saved, scored_regions
from spatialfloors.metrics.f5 import CLASSES

R = Path(sys.argv[1]); U = sys.argv[2:]
L = json.loads(Path("manifests/locked_v3.json").read_text())["locked"]; L1 = json.loads(Path("manifests/locked.json").read_text())["locked"]
sp_ = lambda u: "mouse" if u.startswith("allen") else "human"; cfg = lambda f, m, u: L[f][m][sp_(u)]["config"]
donors = {f"gse303070_{d['slide'].replace('Tumor', 'Tumour')}_{d['pref']}": d["donor_cells"] for d in json.loads(Path("manifests/gse303070_donors.json").read_text())["check_vs_fov_csv"]}
unit_of = lambda c: donors.get(c, c)
idx = json.loads(Path("manifests/crops/index.json").read_text()); MK = json.loads(Path("manifests/panels/index.json").read_text())["f5_heldout_markers"]
F2U = [u for u in U if u.startswith(("gse1_", "allen_"))]; F3U = [u for u in U if not u.startswith("spatch")]; F5U = [u for u in U if not u.startswith(("allen", "spatch"))]
jl = lambda p: json.loads(Path(p).read_text()) if Path(p).exists() else None
out = dict(runs=str(R), comparisons={}, f2_per_crop={}, novae_collapse={}, f5_models={})


def compare(name, fam, vals):
    """vals: crop -> (first, second); D = first - second, unit means, one-sided 95% bounds (descriptive)."""
    vals = {c: v for c, v in vals.items() if v[0] is not None and v[1] is not None and np.isfinite(v[0]) and np.isfinite(v[1])}
    if not vals: return
    units, D, strata = stats.aggregate_units({c: a - b for c, (a, b) in vals.items()}, {c: unit_of(c) for c in vals})
    t = stats.family_test(D, strata, margin=L[fam]["margin"]) if len(D) > 1 else dict(est=float(D[0]), lower=np.nan, upper=np.nan)
    m = {k: float(stats.aggregate_units({c: v[i] for c, v in vals.items()}, {c: unit_of(c) for c in vals})[1].mean()) for i, k in enumerate(("first", "second"))}
    out["comparisons"][name] = dict(family=fam, n_units=len(D), est=t["est"], lower=t["lower"], upper=t["upper"], mean_first=m["first"], mean_second=m["second"],
                                    per_unit=dict(zip(units, map(float, D))))


# F2: balanced ARI of every saved partition on the scored cells
for u in F2U:
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad", backed="r"); mask, reg = scored_regions(a.obs.region.astype(str).values, idx[u]["present"]); a.file.close()
    s = {}
    for name, path, key in (("floor", R / f"labels/f2_{u}_floor.npz", cfg("F2", "floor", u)), ("floor_multiscale", R / f"f2v3/labels/f2_{u}_floor.npz", L["F2"]["floor"][sp_(u)]["secondary_v3"]),
                            ("cellcharter", R / f"labels/f2_{u}_cellcharter.npz", cfg("F2", "cellcharter", u)), ("cellcharter_pca", R / f"labels/f2_{u}_cellcharter_pca.npz", "pc50_L3"),
                            ("novae_zero_shot", R / f"labels/f2_{u}_novae.npz", "zero_shot")):
        if path.exists() and key in np.load(path, allow_pickle=True).files: s[name] = ari_bal_saved(path, key, mask, reg)
    out["f2_per_crop"][u] = s
    n = jl(R / f"f2_{u}_novae.json")
    if n: z = n["grid"]["zero_shot"]; out["novae_collapse"][u] = dict(effective_rank=z["effective_rank"], largest_domain_share=z["largest_domain_share"], collapse=bool(z["collapse"]),
                                                                       rule=bool(z["effective_rank"] < 5 or z["largest_domain_share"] > .5), test_pair=sp_(u) == "mouse")
f2 = out["f2_per_crop"]; g = lambda u, k: f2[u].get(k)
compare("F2 CellCharter - multi-scale floor", "F2", {u: (g(u, "cellcharter"), g(u, "floor_multiscale")) for u in F2U})
compare("F2 Novae - multi-scale floor (mouse)", "F2", {u: (g(u, "novae_zero_shot"), g(u, "floor_multiscale")) for u in F2U if sp_(u) == "mouse"})
compare("F2 multi-scale floor - locked floor", "F2", {u: (g(u, "floor_multiscale"), g(u, "floor")) for u in F2U})
compare("F2 CellCharter on PCA - locked floor", "F2", {u: (g(u, "cellcharter_pca"), g(u, "floor")) for u in F2U})
compare("F2 CellCharter (scVI) - CellCharter on PCA", "F2", {u: (g(u, "cellcharter"), g(u, "cellcharter_pca")) for u in F2U})

# F3: v2 floor, ED-NCEM
for u in F3U:
    v3, v2, lg, nc = (jl(p) for p in (R / f"f3_{u}_ridge.json", R / f"v2/f3_{u}_ridge.json", R / f"f3_{u}_lightgbm.json", R / f"f3_{u}_ncemA.json"))
    if not (v3 and v2 and lg and nc): continue
    ed, li = L["F3"]["ncemA"][sp_(u)]["config"].split(","); s = lambda d, c: d["grid"][c]["unit"]
    out.setdefault("f3_per_crop", {})[u] = dict(floor_v3=s(v3, cfg("F3", "ridge", u)), floor_v2=s(v2, L["F3"]["ridge"][sp_(u)]["v2"]), lightgbm=s(lg, cfg("F3", "lightgbm", u)), ncem_linear=s(nc, li), ncem_ed=s(nc, ed))
f3 = out.get("f3_per_crop", {})
compare("F3 LightGBM - v2 floor", "F3", {u: (v["lightgbm"], v["floor_v2"]) for u, v in f3.items()})
compare("F3 NCEM - v2 floor", "F3", {u: (v["ncem_linear"], v["floor_v2"]) for u, v in f3.items()})
compare("F3 ED-NCEM - floor", "F3", {u: (v["ncem_ed"], v["floor_v3"]) for u, v in f3.items()})
compare("F3 ED-NCEM - NCEM (linear interactions)", "F3", {u: (v["ncem_ed"], v["ncem_linear"]) for u, v in f3.items()})

# F5: floor variants, IF-gated macro-F1, marker-positive recall, class shares
LR = max((c for c in L1["F5"]["floors"]["human"]["grid_means"] if c.startswith("logreg")), key=lambda c: L1["F5"]["floors"]["human"]["grid_means"][c])
MODELS = {"floor": ("floors", cfg("F5", "floors", "gse1")), "insitutype_no_prior": ("floors", "insitutype"), "logreg": ("floors", LR),
          "scanvi": ("scanvi", cfg("F5", "scanvi", "gse1")), "tangram": ("tangram", cfg("F5", "tangram", "gse1"))}  # F5 is human only
COARSE0 = dict(F5.COARSE)
for u in F5U:
    js = {f: jl(R / f"f5_{u}_{f}.json") for f in ("floors", "scanvi", "tangram")}
    if not all(js.values()): continue
    F5.COARSE.clear(); F5.COARSE.update(COARSE0)  # for_unit extends the fine type -> class map with the unit's own reference
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); mk = F5.for_unit(u, MK); keep = filter_reason(a.var_names, "human") == ""
    Xf, vf = sp.csr_matrix(a.X[:, keep]), list(map(str, a.var_names[keep])); lib = np.maximum(np.asarray(Xf.sum(1)).ravel(), 1); del a
    hs = {c: np.log1p(np.asarray(Xf[:, [vf.index(x) for x in mk[c] if x in vf]].sum(1)).ravel() / lib * 1e4) for c in CLASSES}
    top = {c: hs[c] >= np.quantile(hs[c], .9) for c in CLASSES}; mpos = {c: top[c] & ~np.any([top[o] for o in CLASSES if o != c], 0) & (hs[c] > 0) for c in CLASSES}
    for m, (f, c) in MODELS.items():
        r = js[f]["grid"][c]; lab = F5.coarse(np.load(R / f"labels/f5_{u}_{f}.npz", allow_pickle=True)[c].astype(str))  # saved labels are reference types
        out["f5_models"].setdefault(m, {})[u] = dict(config=c, heldout_auroc=r["heldout_auroc"], per_class_auroc={k: v["auroc"] for k, v in r["per_class"].items()}, class_share=r["pred_share"],
                                                     gated_macro_f1=r.get("gated_macro_f1"), marker_pos_recall={k: float((lab[mpos[k]] == k).mean()) if mpos[k].any() else None for k in CLASSES})
        assert np.allclose([float((lab == k).mean()) for k in CLASSES], [r["pred_share"][k] for k in CLASSES], atol=1e-3), f"{u} {m}: class map differs from the run's"
fm = out["f5_models"]; h = lambda m, u, k="heldout_auroc": fm[m][u][k]
compare("F5 scANVI - InSituType without prior (v2 floor)", "F5", {u: (h("scanvi", u), h("insitutype_no_prior", u)) for u in fm.get("floor", {})})
compare("F5 Tangram - InSituType without prior (v2 floor)", "F5", {u: (h("tangram", u), h("insitutype_no_prior", u)) for u in fm.get("floor", {})})
compare("F5 floor - InSituType without prior (spatial prior gain)", "F5", {u: (h("floor", u), h("insitutype_no_prior", u)) for u in fm.get("floor", {})})
compare(f"F5 scANVI - logistic regression ({LR})", "F5", {u: (h("scanvi", u), h("logreg", u)) for u in fm.get("floor", {})})
compare(f"F5 Tangram - logistic regression ({LR})", "F5", {u: (h("tangram", u), h("logreg", u)) for u in fm.get("floor", {})})
compare(f"F5 logistic regression ({LR}) - floor", "F5", {u: (h("logreg", u), h("floor", u)) for u in fm.get("floor", {})})
for m in ("scanvi", "tangram", "insitutype_no_prior", "logreg"):
    compare(f"F5 IF-gated macro-F1 (GSE1): {m} - floor", "F5", {u: (h(m, u, "gated_macro_f1"), h("floor", u, "gated_macro_f1")) for u in fm.get("floor", {}) if u.startswith("gse1_")})
out["f5_summary"] = {m: dict(n_crops=len(v), heldout_auroc=float(np.mean([x["heldout_auroc"] for x in v.values()])),
                             gated_macro_f1_gse1=float(np.mean([x["gated_macro_f1"] for u, x in v.items() if u.startswith("gse1_")])),
                             marker_pos_recall={c: float(np.mean([x["marker_pos_recall"][c] for x in v.values() if x["marker_pos_recall"][c] is not None])) for c in CLASSES},
                             marker_pos_recall_n_crops={c: sum(x["marker_pos_recall"][c] is not None for x in v.values()) for c in CLASSES},
                             class_share={c: float(np.mean([x["class_share"][c] for x in v.values()])) for c in (*CLASSES, "other")}) for m, v in fm.items()}

(R / "P5_secondary.json").write_text(json.dumps(out, indent=1, default=float))
md = ["Secondary outputs (PROTOCOL §2; no decision weight). D = first - second named; unit means with one-sided 95% bounds.", "",
      "| comparison | units | mean first | mean second | D | lower | upper |", "|---|---|---|---|---|---|---|"]
md += [f"| {k} | {v['n_units']} | {v['mean_first']:.3f} | {v['mean_second']:.3f} | {v['est']:+.4f} | {v['lower']:+.4f} | {v['upper']:+.4f} |" for k, v in out["comparisons"].items()]
nr = out["f5_summary"]["floor"]["marker_pos_recall_n_crops"]
md += ["", "F5 per model, mean over crops. Marker-positive recall: share of the cells in the top 10% of one class's held-out-marker score",
       f"(and of no other class's) that the model assigns to that class; defined on {nr['epithelial']} / {nr['immune_nonmyeloid']} / {nr['myeloid']} crops"
       " (undefined where another class's 90th percentile is 0). Class share: predicted fraction.", "",
       "| model | config | AUROC | IF-gated macro-F1 (GSE1) | recall epi / imm / mye | share epi / imm / mye / other |", "|---|---|---|---|---|---|"]
j = lambda d, ks: " / ".join(f"{d[c]:.2f}" for c in ks)
md += [f"| {m} | {MODELS[m][1]} | {v['heldout_auroc']:.3f} | {v['gated_macro_f1_gse1']:.3f} | {j(v['marker_pos_recall'], CLASSES)} | {j(v['class_share'], (*CLASSES, 'other'))} |"
       for m, v in out["f5_summary"].items()]
md += ["", "Novae zero-shot collapse rule (effective rank < 5 or > 50% of cells in one domain; test pair on mouse units only)", "",
       "| crop | effective rank | largest domain share | collapsed |", "|---|---|---|---|"]
md += [f"| {u} | {v['effective_rank']:.2f} | {v['largest_domain_share']:.2f} | {v['rule']} |" for u, v in out["novae_collapse"].items()]
md += ["", "F2 balanced ARI per crop", "", "| crop | " + " | ".join(("floor", "floor_multiscale", "cellcharter", "cellcharter_pca", "novae_zero_shot")) + " |", "|---" * 6 + "|"]
md += [f"| {u} | " + " | ".join(f"{v[k]:.3f}" if k in v else "-" for k in ("floor", "floor_multiscale", "cellcharter", "cellcharter_pca", "novae_zero_shot")) + " |" for u, v in f2.items()]
(R / "P5_secondary.md").write_text("\n".join(md) + "\n"); print(len(out["comparisons"]), "comparisons ->", R / "P5_secondary.md")
