"""Tier-1 analysis (PROTOCOL §4), written and tested on development and validation outputs before the confirmation
outputs were unsealed, then run once on runs/p4. Per (family, learned model) pair: per-crop primary of the locked floor
and the locked learned
config -> D = learned - floor -> one D per unit (stats.aggregate_units: GSE303070 donors, Allen F4 animals) ->
stats.family_test
(non-inferiority = one-sided 95% upper bound < margin; TOST equivalence; floor superior) -> Holm over pairs for learned
superiority; intersection-union over all pairs. Sensitivity: leave-one-cancer-out (pan-cancer families) and Wilcoxon
(superiority, non-inferiority). Secondary: the F3 decomposition.
Primaries: F1 rescored kNN-Poisson score (p3_f1_rescore.py), F2 balanced ARI, F3 top-10 Spearman,
F4 official and spike-free pearson (co-primary; seed mean over tasks), F5 held-out-marker AUROC.
Usage (env floors): python scripts/p5_analyze.py <runs dir> <unit ...> [retro:<sid|section> ...] -> <runs dir>/P5.json +
P5.md
(retro: = F4 retrospective slide / section: GSE1 120 210 221 242, Allen animal 2 sections)"""

import glob, json, sys
from pathlib import Path
import numpy as np, pandas as pd
from spatialfloors import stats

R = Path(sys.argv[1])
U = [x for x in sys.argv[2:] if not x.startswith("retro:")]
RETRO = [x[6:] for x in sys.argv[2:] if x.startswith("retro:")]
L = json.loads(Path(__import__("os").environ.get("SF_LOCKED", "manifests/locked_v3.json")).read_text())[
    "locked"
]  # SF_LOCKED: an earlier lock (locked.json, locked_v2.json)
sp = lambda u: "mouse" if u.startswith("allen") else "human"
cfg = lambda f, m, u: L[f][m][sp(u)]["config"]
NC = (
    json.loads(Path("manifests/locked_v3_ncem.json").read_text())["locked"] if "ncemA" in L["F3"] else None
)  # PROTOCOL §2: the authors' NCEM test configuration
F3M = (
    {"lightgbm": ("lightgbm", lambda u: cfg("F3", "lightgbm", u)), "ncem": ("ncemA", lambda u: NC[sp(u)]["config"])}
    if NC
    else {m: (m, lambda u, m=m: cfg("F3", m, u)) for m in ("gat", "lightgbm", "ncem") if "human" in L["F3"].get(m, {})}
)
CANCER = {
    "gse1": "CRC",
    "gse0": "CRC",
    "gse303070": "CRC",
    "spatch": "SPATCH (tissue unverified)",
    "allen": "mouse brain",
    "bruker_breast": "breast",
    "bruker_mobreast": "breast",
    "bruker_skin": "skin SCC",
    "bruker_kidney": "RCC",
}
cancer = lambda u: next(v for k, v in sorted(CANCER.items(), key=lambda kv: -len(kv[0])) if u.startswith(k))
donors = {
    f"gse303070_{d['slide'].replace('Tumor', 'Tumour')}_{d['pref']}": d["donor_cells"]
    for d in json.loads(Path("manifests/gse303070_donors.json").read_text())["check_vs_fov_csv"]
}
unit_of = lambda c: donors.get(c, c)
EXP = (
    json.loads(Path(__import__("os").environ["SF_EXPECT"]).read_text())
    if __import__("os").environ.get("SF_EXPECT")
    else None
)  # fail if a registered pair or unit is missing
rows = []  # family, model, crop, floor, learned


def grid(f, u, m, key):
    p = R / f"{f.lower()}_{u}_{m}.json"
    return json.loads(p.read_text())["grid"] if p.exists() else None


# F1
f1 = sorted(glob.glob(str(R / "f1_rescored*.csv")))  # development + validation tables
if f1:
    d = pd.concat([pd.read_csv(x) for x in f1]).drop_duplicates(["unit", "model", "config"], keep="last")
    for u in U:
        s = lambda m: d[(d.unit == u) & (d.model == m) & (d.config == cfg("F1", m, u))].score
        fl = s("expr_embed")
        for m in ("scvi", "nicheformer"):
            le = s(m)
            if len(fl) and len(le):
                rows.append(dict(fam="F1", model=m, crop=u, floor=float(fl.iloc[0]), learned=float(le.iloc[0])))
# F2 balanced ARI from saved labels (seeds x cells)
import anndata as ad
from spatialfloors.metrics.f2 import scored_regions, ari_bal, ari_bal_saved

idx = json.loads(Path("manifests/crops/index.json").read_text())
for u in U:
    fl = R / f"labels/f2_{u}_floor.npz"
    if not fl.exists():
        continue
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad", backed="r")
    mask, reg = scored_regions(a.obs.region.astype(str).values, idx[u]["present"])
    a.file.close()
    bal = lambda f, key: ari_bal_saved(f, key, mask, reg)
    flv = bal(fl, cfg("F2", "floor", u))
    from sklearn.cluster import KMeans  # geometry reference: secondary, per crop

    a2 = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad", backed="r")
    xy = a2.obs[["x_um", "y_um"]].values
    a2.file.close()
    k_ = len(np.unique(reg))
    blob = float(
        np.mean(
            [ari_bal(reg, KMeans(k_, n_init=10, random_state=s_).fit_predict(xy)[mask].astype(str)) for s_ in range(3)]
        )
    )
    rows.append(dict(fam="F2ref", model="xy_blob", crop=u, floor=flv, learned=blob))
    for m in ("cellcharter", "novae"):
        lf = R / f"labels/f2_{u}_{m}.npz"
        try:
            c = cfg("F2", m, u)
        except KeyError:
            continue  # Novae has no human lock (both configs collapsed, reported separately)
        if lf.exists() and c in np.load(lf, allow_pickle=True).files:
            rows.append(dict(fam="F2", model=m, crop=u, floor=flv, learned=bal(lf, c)))
# F3, F5
for f, fm, key, ms in (
    ("F3", "ridge", "unit", F3M),
    ("F5", "floors", "heldout_auroc", {m: (m, lambda u, m=m: cfg("F5", m, u)) for m in ("scanvi", "tangram")}),
):
    for u in U:
        g0 = grid(f, u, fm, key)
        if not g0:
            continue
        for m, (fn, c) in ms.items():
            g = grid(f, u, fn, key)
            if g:
                rows.append(dict(fam=f, model=m, crop=u, floor=g0[cfg(f, fm, u)][key], learned=g[c(u)][key]))
# F3 decomposition (PROTOCOL §2, secondary, no decision weight): NCEM vs a ridge on its own inputs ('presence' = the
# linear
# interactions design, same graph and folds), NCEM vs the v2 floor (other lineages only, no rings), v3 floor vs v2 floor
sec = []
for u in U if NC else []:
    nf, v3f, v2f = R / f"f3_{u}_ncemA.json", R / f"f3_{u}_ridge.json", R / f"v2/f3_{u}_ridge.json"
    if not nf.exists():
        continue
    n = json.loads(nf.read_text())
    c = NC[sp(u)]["config"]
    s_nc = n["grid"][c]["unit"]
    ir = n.get("ridge", {}).get(f"presence_r{c.rsplit('_r', 1)[1]}", {}).get("unit")
    v2 = json.loads(v2f.read_text())["grid"][L["F3"]["ridge"][sp(u)]["v2"]]["unit"] if v2f.exists() else None
    v3 = json.loads(v3f.read_text())["grid"][cfg("F3", "ridge", u)]["unit"] if v3f.exists() else None
    for m, fl, le in (
        ("ncem_vs_input_ridge", ir, s_nc),
        ("ncem_vs_v2_floor", v2, s_nc),
        ("v3_floor_vs_v2_floor", v2, v3),
    ):
        if fl is not None and le is not None:
            sec.append(dict(fam="F3", model=m, crop=u, floor=float(fl), learned=float(le)))
# F4: prospective units scored by f4_score.py (<R>/f4/f4_<sid>.csv); unit = slide (GSE1) or animal (MERFISH sections)
AL = json.loads(Path("manifests/allen.json").read_text())
sec_animal = {s: an for an, v in AL.items() if isinstance(v, dict) for s in v.get("f4_sections", [])}
for f in sorted(glob.glob(str(R / "f4/f4_*.csv"))):
    d = pd.read_csv(f)
    sid = str(d.sid.iloc[0])
    crop = f"gse1_{sid}" if sid.isdigit() else f"allen_{sec_animal[sid]}_{sid}"
    for co in ("pearson", "pearson_clean"):
        fl = d[d.model == "L6"][co].mean()
        for m in ("cellina", "cpa"):
            le = d[d.model == m].groupby("holdout")[co].mean().mean()  # seed mean per task, then task mean
            if np.isfinite(le):
                rows.append(dict(fam=f"F4_{co}", model=m, crop=crop, floor=float(fl), learned=float(le)))
# F4 retrospective (PROTOCOL §3): Cellina / CPA = the authors' camera-ready official pearson; L6 = our run on the
# same slide / section (<R>/f4/l6). Spike-free needs counterfactual outputs, so retrospective units enter F4_pearson
# only.
CR = Path("external/cellina-reproducibility/results")
for sid in RETRO:
    mer = not sid.isdigit()
    l6f = R / f"f4/l6/{sid if mer else 'crc_' + sid}.csv"
    if not l6f.exists():
        continue
    cam = pd.read_csv(CR / f"loo_summary_{'merfish' if mer else 'crc'}_DEG_50_camready.csv")
    cam = cam[cam.sid == (sid if mer else f"crc_{sid}")]
    l6v = pd.read_csv(l6f)
    l6v = l6v[l6v.model_name == "L6_twoway_pseudobulk"].set_index("holdout_celltype").pearson
    crop = f"allen_{sec_animal[sid]}_{sid}" if mer else f"gse1_{sid}"
    for m in ("cellina", "cpa"):
        c = cam[cam.model_name == m].set_index("holdout_celltype").pearson
        t = c.index.intersection(l6v.index)
        if len(t):
            rows.append(
                dict(
                    fam="F4_pearson",
                    model=m,
                    crop=crop,
                    floor=float(l6v[t].mean()),
                    learned=float(c[t].mean()),
                    retrospective=True,
                )
            )
df = pd.DataFrame(rows)
out = dict(runs=str(R), units=U, pairs={})
if df.empty:
    print("no pairs")
    sys.exit(1)
df["D"] = df.learned - df.floor
out["f2_blob_reference"] = {r.crop: dict(floor=r.floor, blob=r.learned) for r in df[df.fam == "F2ref"].itertuples()}
df = df[df.fam != "F2ref"]
for (f, m), g in df.groupby(["fam", "model"]):
    mg = L[f.split("_")[0]]["margin"]
    uo = (
        (lambda c: c.rsplit("_", 1)[0] if c.startswith("allen_") else c) if f.startswith("F4") else unit_of
    )  # F4 MERFISH unit = animal
    units, D, strata = stats.aggregate_units(dict(zip(g.crop, g.D)), {c: uo(c) for c in g.crop})
    t = stats.family_test(D, strata, margin=mg) if len(D) > 1 else dict(est=float(D[0]), note="one unit: descriptive")
    if EXP and f"{f}|{m}" in EXP["n_units"]:
        assert len(units) == EXP["n_units"][f"{f}|{m}"], (
            f"{f}|{m}: {len(units)} units scored, {EXP['n_units'][f'{f}|{m}']} registered"
        )
    um = {
        k: float(stats.aggregate_units(dict(zip(g.crop, g[k])), {c: uo(c) for c in g.crop})[1].mean())
        for k in ("floor", "learned")
    }  # unit means, as D
    t.update(
        margin=mg,
        n_units=len(D),
        per_unit=dict(zip(units, map(float, D))),
        floor_mean=um["floor"],
        learned_mean=um["learned"],
    )
    canc = {c: cancer(c) for c in g.crop}
    if len(set(canc.values())) > 1 and len(D) > 2:  # leave-one-cancer-out
        t["leave_one_cancer_out"] = {}
        for k in sorted(set(canc.values())):
            gg = g[[canc[c] != k for c in g.crop]]
            if len(gg) < 2:
                continue
            u2, D2, s2 = stats.aggregate_units(dict(zip(gg.crop, gg.D)), {c: uo(c) for c in gg.crop})
            if len(D2) > 1:
                t["leave_one_cancer_out"][k] = {
                    kk: v
                    for kk, v in stats.family_test(D2, s2, margin=mg).items()
                    if kk in ("est", "upper", "noninferior")
                }
    t["by_stratum"] = {
        k: float(np.mean(v)) for k, v in pd.Series(D, index=strata).groupby(level=0)
    }  # moderators, descriptive
    if (f, m) == ("F1", "nicheformer"):  # pretraining overlap: Nicheformer's notebooks use the Yao MERFISH animal
        gg = g[[not c.startswith("allen_yao") for c in g.crop]]
        if len(gg) < len(g) and len(gg) > 2:
            u2, D2, s2 = stats.aggregate_units(dict(zip(gg.crop, gg.D)), {c: uo(c) for c in gg.crop})
            t["without_pretraining_overlap"] = {
                kk: v for kk, v in stats.family_test(D2, s2, margin=mg).items() if kk in ("est", "upper", "noninferior")
            }
    out["pairs"][f"{f}|{m}"] = t
if EXP:
    assert set(EXP["n_units"]) <= set(out["pairs"]), (
        f"registered pairs missing: {sorted(set(EXP['n_units']) - set(out['pairs']))}"
    )
keys = [k for k, v in out["pairs"].items() if "p_superior" in v]
for k, r in zip(keys, stats.holm([out["pairs"][k]["p_superior"] for k in keys])):
    out["pairs"][k]["learned_superior_holm"] = bool(r)
out["floor_noninferior_all"] = all(v.get("noninferior", False) for k, v in out["pairs"].items() if "noninferior" in v)
out["secondary"] = {}
for (f, m), g in pd.DataFrame(sec).groupby(["fam", "model"]) if sec else []:
    units, D, strata = stats.aggregate_units(dict(zip(g.crop, g.learned - g.floor)), {c: unit_of(c) for c in g.crop})
    t = stats.family_test(D, strata, margin=L[f]["margin"]) if len(D) > 1 else dict(est=float(D[0]))
    out["secondary"][f"{f}|{m}"] = {k: t[k] for k in ("est", "lower", "upper") if k in t} | dict(
        n_units=len(D), per_unit=dict(zip(units, map(float, D)))
    )
NAME = __import__("os").environ.get("SF_P5NAME", "P5")
(R / f"{NAME}.json").write_text(json.dumps(out, indent=1, default=float))
md = [
    "| pair | units | floor | learned | D (est) | upper | margin | non-inferior | equivalent | floor superior |",
    "|---|---|---|---|---|---|---|---|---|---|",
]
for k, v in out["pairs"].items():
    md.append(
        f"| {k} | {v['n_units']} | {v['floor_mean']:.3f} | {v['learned_mean']:.3f} | {v['est']:+.4f} | {v.get('upper', float('nan')):+.4f} | {v['margin']:.4f} | {v.get('noninferior', '-')} | {v.get('equivalent', '-')} | {v.get('floor_superior', '-')} |"
    )
md += [
    f"\nfloor non-inferior on every pair: {out['floor_noninferior_all']}",
    "",
    "Secondary (no decision weight): F3 decomposition, D = first - second named",
    "",
    "| comparison | units | D (est) | lower | upper |",
    "|---|---|---|---|---|",
]
md += [
    f"| {k} | {v['n_units']} | {v['est']:+.4f} | {v.get('lower', float('nan')):+.4f} | {v.get('upper', float('nan')):+.4f} |"
    for k, v in out["secondary"].items()
]
(R / f"{NAME}.md").write_text("\n".join(md) + "\n")
print("\n".join(md[:3]), "...", len(out["pairs"]), "pairs,", len(out["secondary"]), "secondary ->", R / f"{NAME}.md")
