"""Replication-round analysis (registration v2): the registered round-1 analysis (scripts/p5_analyze.py, PROTOCOL §4)
with the crops of one role, their units (patient / donor) and the leave-one-cancer-out groups read from a unit table.
Per (family, learned model) pair: per-crop primary of the locked floor and the locked learned configuration ->
D = learned - floor -> one D per unit (stats.aggregate_units) -> stats.family_test (non-inferiority = one-sided 95% upper
bound < margin; TOST equivalence; floor superior) -> Holm over pairs for learned superiority; intersection-union over
all pairs. Sensitivity: leave-one-cancer-out, Wilcoxon. Secondary: the F3 decomposition. F4 (round 1 only): units and
retrospective slides by the round-1 rule.
Primaries: F1 rescored kNN-Poisson score, F2 balanced ARI, F3 top-10 Spearman, F4 official and spike-free pearson,
F5 held-out-marker AUROC.
Usage (env floors): python scripts/v2_analyze.py --units <units.tsv> --role <role> --runs <runs dir> [--root <dir>]
    [--out <stem>] [--loco <column>] [--expect <json>] [--retro <sid|section> ...]
  --root: directory holding runs/p2/cache, manifests/crops/index.json and external/ (default: working directory).
  Locked configurations and margins: this repository's manifests/locked_v3.json and locked_v3_ncem.json.
  --loco: unit-table column for leave-one-cancer-out (default: 'tissue' if present, else the round-1 labels).
  --out: output stem (default <runs dir>/V2) -> <stem>.json, <stem>.md.
Round-1 check (from the round-1 runs root): --units registration/units/units.tsv --role confirmation --runs runs/p4
    --expect manifests/expected_units_tier1.json --retro 120 210 221 242 C57BL6J-2.036 C57BL6J-2.039 C57BL6J-2.041
  gives results/reports/P5_corrected.md byte for byte and P5_corrected.json with every value identical; per-unit keys
  of the four Allen animals carry the table's unit name instead of the crop name (one crop per animal)."""

import argparse, glob, json
from pathlib import Path
import numpy as np, pandas as pd
from spatialfloors import stats

REPO = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--units", required=True)
ap.add_argument("--role", required=True)
ap.add_argument("--runs", required=True)
ap.add_argument("--root", default=".")
ap.add_argument("--out")
ap.add_argument("--loco")
ap.add_argument("--expect")
ap.add_argument("--retro", nargs="*", default=[])
A = ap.parse_args()
R, ROOT = Path(A.runs), Path(A.root)

T = pd.read_csv(A.units, sep="\t", dtype=str, keep_default_na=False)
T = T[T.role == A.role]
assert len(T), f"no crops with role {A.role}"
U = T.crop.tolist()
for c, s in zip(T.crop, T.stratum):
    assert stats.stratum(c) == s, f"{c}: table stratum {s}, stats.stratum {stats.stratum(c)}"
UNIT = dict(zip(T.crop, T.unit))
unit_of = lambda c: UNIT.get(c, c)
LOCO = A.loco or ("tissue" if "tissue" in T.columns else None)
CANCER_V1 = {  # round-1 registered labels, used when the table has no leave-one-cancer-out column
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
CMAP = dict(zip(T.crop, T[LOCO])) if LOCO else {}
cancer = lambda u: CMAP.get(u) or next(
    v for k, v in sorted(CANCER_V1.items(), key=lambda kv: -len(kv[0])) if u.startswith(k)
)

L = json.loads((REPO / "manifests/locked_v3.json").read_text())["locked"]
sp = lambda u: "mouse" if u.startswith("allen") else "human"
cfg = lambda f, m, u: L[f][m][sp(u)]["config"]
NC = json.loads((REPO / "manifests/locked_v3_ncem.json").read_text())["locked"]  # authors' NCEM test configuration
F3M = {"lightgbm": ("lightgbm", lambda u: cfg("F3", "lightgbm", u)), "ncem": ("ncemA", lambda u: NC[sp(u)]["config"])}
EXP = json.loads(Path(A.expect).read_text()) if A.expect else None  # fail if a registered pair or unit is missing
rows = []  # family, model, crop, floor, learned


def grid(f, u, m):
    p = R / f"{f.lower()}_{u}_{m}.json"
    return json.loads(p.read_text())["grid"] if p.exists() else None


# F1
f1 = sorted(glob.glob(str(R / "f1_rescored*.csv")))
if f1:
    d = pd.concat([pd.read_csv(x) for x in f1]).drop_duplicates(["unit", "model", "config"], keep="last")
    for u in U:
        s = lambda m: d[(d.unit == u) & (d.model == m) & (d.config == cfg("F1", m, u))].score
        fl = s("expr_embed")
        for m in ("scvi", "nicheformer"):
            le = s(m)
            if len(fl) and len(le):
                rows.append(dict(fam="F1", model=m, crop=u, floor=float(fl.iloc[0]), learned=float(le.iloc[0])))
# F2 balanced ARI from saved labels (seeds x cells); xy k-means = geometry reference, per crop
F2U = [u for u in U if (R / f"labels/f2_{u}_floor.npz").exists()]
if F2U:
    import anndata as ad
    from sklearn.cluster import KMeans
    from spatialfloors.metrics.f2 import scored_regions, ari_bal, ari_bal_saved

    idx = json.loads((ROOT / "manifests/crops/index.json").read_text())
for u in F2U:
    fl = R / f"labels/f2_{u}_floor.npz"
    a = ad.read_h5ad(ROOT / f"runs/p2/cache/{u}.h5ad", backed="r")
    mask, reg = scored_regions(a.obs.region.astype(str).values, idx[u]["present"])
    xy = a.obs[["x_um", "y_um"]].values
    a.file.close()
    bal = lambda f, key: ari_bal_saved(f, key, mask, reg)
    flv = bal(fl, cfg("F2", "floor", u))
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
            continue  # Novae has no human lock
        if lf.exists() and c in np.load(lf, allow_pickle=True).files:
            rows.append(dict(fam="F2", model=m, crop=u, floor=flv, learned=bal(lf, c)))
# F3, F5
for f, fm, key, ms in (
    ("F3", "ridge", "unit", F3M),
    ("F5", "floors", "heldout_auroc", {m: (m, lambda u, m=m: cfg("F5", m, u)) for m in ("scanvi", "tangram")}),
):
    for u in U:
        g0 = grid(f, u, fm)
        if not g0:
            continue
        for m, (fn, c) in ms.items():
            g = grid(f, u, fn)
            if g:
                rows.append(dict(fam=f, model=m, crop=u, floor=g0[cfg(f, fm, u)][key], learned=g[c(u)][key]))
# F3 decomposition (secondary, no decision weight): NCEM vs a ridge on its own inputs, NCEM vs the v2 floor (other
# lineages only, no rings), v3 floor vs v2 floor
sec = []
for u in U:
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
# F4 (round 1): prospective units scored by f4_score.py (<runs>/f4/f4_<sid>.csv); unit = slide (GSE1) or animal (MERFISH)
AL = json.loads((REPO / "manifests/allen.json").read_text())
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
# F4 retrospective: Cellina / CPA = the authors' camera-ready official pearson; L6 = our run on the same slide / section
CR = ROOT / "external/cellina-reproducibility/results"
for sid in A.retro:
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
    raise SystemExit("no pairs")
df["D"] = df.learned - df.floor
out["f2_blob_reference"] = {r.crop: dict(floor=r.floor, blob=r.learned) for r in df[df.fam == "F2ref"].itertuples()}
df = df[df.fam != "F2ref"]
for (f, m), g in df.groupby(["fam", "model"]):
    mg = L[f.split("_")[0]]["margin"]
    uo = (lambda c: c.rsplit("_", 1)[0] if c.startswith("allen_") else c) if f.startswith("F4") else unit_of
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
    t["by_stratum"] = {k: float(np.mean(v)) for k, v in pd.Series(D, index=strata).groupby(level=0)}  # descriptive
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
stem = Path(A.out) if A.out else R / "V2"
stem.with_suffix(".json").write_text(json.dumps(out, indent=1, default=float))
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
stem.with_suffix(".md").write_text("\n".join(md) + "\n")
print(len(out["pairs"]), "pairs,", len(out["secondary"]), "secondary ->", stem.with_suffix(".json"))
