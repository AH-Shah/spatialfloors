"""Exploratory: baseline ladder and assay sensitivity on the round-1 confirmation units. Each family gets one or more
naive baselines below the floor; the comparison can only be read as "learned ≈ floor" if it separates both from them.
  F1  log-normalised PCA 50 (f1_fit.py expr_pca pc50; new run)        < floor (Pearson-residual PCA), scVI, Nicheformer
  F2  k-means on coordinates (tier-1 report); GMM on expression PCA 50 without any neighbourhood (new run, mode f2)
  F3  ridge on the neighbourhood's type mix only (p3_f3.py T:ell25:pc10:a1000; new run); the floor without rings (v2)
  F4  the simpler rungs of the floor's own ladder (control mean, L1-L5; saved with the L6 runs)
  F5  logistic regression on the reference (secondary report); InSituType without the spatial prior
Tests: first - second per unit (stats.aggregate_units), registered test (stats.family_test, platform strata);
one-sided 95% bounds; "ahead" = lower bound > 0. F4 learned scores = L6 + the registered D of the unit.
Usage (cwd = the runs root; env floors):
  python scripts/x_ladder.py f2 <crop> <ladder dir>              -> <ladder dir>/f2_<crop>_expr_gmm.{json,npz}
  python scripts/x_ladder.py report <ladder dir> <runs>/p4 <out> -> <out>/X_ladder.{json,md}"""

import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from spatialfloors import stats
from spatialfloors.provenance import commit, stamp

PUB = Path(__file__).resolve().parents[1]
T_ONLY = "T:ell25:pc10:a1000"

if sys.argv[1] == "f2":
    import anndata as ad
    from spatialfloors.data.units import index
    from spatialfloors.floors.simple import expr_pca
    from spatialfloors.metrics.f2 import SEEDS, ari_bal, gmm_labels, scored_regions

    u, out = sys.argv[2], Path(sys.argv[3])
    out.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad")
    g = a.var_names.get_indexer(json.loads(Path(f"manifests/panels/{u}.json").read_text())["genes"])
    assert (g >= 0).all()
    m, regl = scored_regions(a.obs.region.astype(str).values, index()[u]["present"])
    k = len(np.unique(regl))
    E = expr_pca(a.X[:, g], 50)
    L = np.stack([gmm_labels(E, k, s) for s in SEEDS])
    v = [ari_bal(regl, l[m].astype(str)) for l in L]
    np.savez_compressed(out / f"f2_{u}_expr_gmm.npz", expr_gmm=L)
    rep = dict(unit=u, model="expr_gmm", k=k, n_scored=int(m.sum()), ari_bal=float(np.mean(v)), per_seed=v)
    (out / f"f2_{u}_expr_gmm.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
    print(u, "expr_gmm", round(rep["ari_bal"], 4), flush=True)
    sys.exit()

LAD, R4, OUT = (Path(x) for x in sys.argv[2:5])
J = lambda p: json.loads(Path(p).read_text())
T1 = J(PUB / "results/reports/P5_corrected.json")
S2 = J(PUB / "results/reports/P5_secondary.json")
LOCK = J(PUB / "manifests/locked_v3.json")["locked"]
donors = {
    f"gse303070_{d['slide'].replace('Tumor', 'Tumour')}_{d['pref']}": d["donor_cells"]
    for d in J(PUB / "manifests/gse303070_donors.json")["check_vs_fov_csv"]
}
unit_of = lambda c: donors.get(c, c)
species = lambda c: "mouse" if c.startswith("allen") else "human"


def test(a, b):
    """first - second on the crops both have, one value per unit, registered test."""
    crops = sorted(set(a) & set(b))
    units, D, S = stats.aggregate_units({c: a[c] - b[c] for c in crops}, {c: unit_of(c) for c in crops})
    t = stats.family_test(D, S if len(set(S)) > 1 else None)
    ma = stats.aggregate_units({c: a[c] for c in crops}, {c: unit_of(c) for c in crops})[1].mean()
    mb = stats.aggregate_units({c: b[c] for c in crops}, {c: unit_of(c) for c in crops})[1].mean()
    return dict(n_units=len(D), mean_first=float(ma), mean_second=float(mb), est=t["est"], lower=t["lower"], upper=t["upper"], ahead=bool(t["lower"] > 0))  # fmt: skip


fam = {}  # family -> list of (label, role, scores per crop)

# F1
f1 = pd.read_csv(R4 / "f1_rescored.csv").drop_duplicates(["unit", "model", "config"], keep="last")
pick = lambda m: {r.unit: r.score for r in f1[f1.model == m].itertuples() if r.config == LOCK["F1"][m][species(r.unit)]["config"]}  # fmt: skip
naive1 = {c: J(LAD / f"f1_{c}_expr_pca_pc50.json")["scores"]["pc50"]["score"] for c in pick("expr_embed") if (LAD / f"f1_{c}_expr_pca_pc50.json").exists()}  # fmt: skip
fam["F1"] = [
    ("log-normalised PCA 50", "naive", naive1),
    ("floor (Pearson-residual PCA 50)", "floor", pick("expr_embed")),
    ("scVI", "learned", pick("scvi")),
    ("Nicheformer", "learned", pick("nicheformer")),
]
assert abs(test(pick("scvi"), pick("expr_embed"))["est"] - T1["pairs"]["F1|scvi"]["est"]) < 1e-9

# F2
f2 = S2["f2_per_crop"]
naive2 = {c: J(LAD / f"f2_{c}_expr_gmm.json")["ari_bal"] for c in f2 if (LAD / f"f2_{c}_expr_gmm.json").exists()}
fam["F2"] = [
    ("k-means on coordinates", "naive", {c: v["blob"] for c, v in T1["f2_blob_reference"].items()}),
    ("GMM on expression PCA 50 (no neighbourhood)", "naive", naive2),
    ("floor (expression + neighbourhood mean)", "floor", {c: v["floor"] for c, v in f2.items()}),
    ("CellCharter", "learned", {c: v["cellcharter"] for c, v in f2.items()}),
    ("Novae zero-shot (mouse units)", "learned", {c: v["novae_zero_shot"] for c, v in f2.items() if species(c) == "mouse"}),
]  # fmt: skip
assert abs(test(fam["F2"][3][2], fam["F2"][2][2])["est"] - T1["pairs"]["F2|cellcharter"]["est"]) < 1e-9

# F3
f3 = S2["f3_per_crop"]
naive3 = {c: J(LAD / f"f3_{c}_ridge.json")["grid"][T_ONLY]["unit"] for c in f3 if (LAD / f"f3_{c}_ridge.json").exists()}
fam["F3"] = [
    ("ridge on the type mix only", "naive", naive3),
    ("floor without rings (v2)", "naive", {c: v["floor_v2"] for c, v in f3.items()}),
    ("floor (v3)", "floor", {c: v["floor_v3"] for c, v in f3.items()}),
    ("LightGBM", "learned", {c: v["lightgbm"] for c, v in f3.items()}),
    ("NCEM", "learned", {c: v["ncem_linear"] for c, v in f3.items()}),
]
assert abs(test(fam["F3"][3][2], fam["F3"][2][2])["est"] - T1["pairs"]["F3|lightgbm"]["est"]) < 1e-9

# F4: the ladder saved with each L6 run; unit = slide (GSE1) or animal (mean over its sections)
AL = J(PUB / "manifests/allen.json")
sec_unit = {s: f"allen_{an}" for an, v in AL.items() if isinstance(v, dict) for s in v.get("f4_sections", [])}
rung = {}
for f in sorted((R4 / "f4/l6").glob("*.csv")):
    sid = f.stem
    u = f"gse1_{sid[4:]}" if sid.startswith("crc_") else sec_unit[sid]
    for mname, v in pd.read_csv(f).groupby("model_name").pearson.mean().items():
        rung.setdefault(mname, {}).setdefault(u, []).append(float(v))
rung = {k: {u: float(np.mean(x)) for u, x in v.items()} for k, v in rung.items()}
l6 = rung["L6_twoway_pseudobulk"]
P4 = T1["pairs"]["F4_pearson|cellina"]
assert sorted(l6) == sorted(P4["per_unit"]) and abs(np.mean(list(l6.values())) - P4["floor_mean"]) < 1e-6
F4N = (
    ("control mean", "baseline_repro"),
    ("L1 shared shift", "L1_shared_shift"),
    ("L2 contamination", "L2_contamination"),
    ("L3 niche ridge", "L3_h_niche_ridge"),
    ("L4 niche match", "L4_niche_match"),
    ("L5 pooled ridge", "L5_pooled_ridge"),
)
fam["F4"] = (
    [(lab, "naive", rung[k]) for lab, k in F4N]
    + [("floor (L6 two-way pseudobulk)", "floor", l6)]
    + [
        (lab, "learned", {u: l6[u] + d for u, d in T1["pairs"][f"F4_pearson|{m}"]["per_unit"].items()})
        for lab, m in (("Cellina", "cellina"), ("CPA", "cpa"))
    ]
)

# F5
f5 = {m: {c: v["heldout_auroc"] for c, v in d.items()} for m, d in S2["f5_models"].items()}
fam["F5"] = [
    ("logistic regression on the reference", "naive", f5["logreg"]),
    ("InSituType without the spatial prior", "naive", f5["insitutype_no_prior"]),
    ("floor (InSituType with the spatial prior)", "floor", f5["floor"]),
    ("scANVI", "learned", f5["scanvi"]),
    ("Tangram", "learned", f5["tangram"]),
]
assert abs(test(f5["scanvi"], f5["floor"])["est"] - T1["pairs"]["F5|scanvi"]["est"]) < 1e-9

out = dict(commit=commit(), families={})
md = [
    "# X_ladder: baseline ladder and assay sensitivity (exploratory)",
    "",
    "Naive baselines below each floor. D = first − second per unit, one-sided 95% bounds; ahead = lower bound > 0.",
]
for f, rungs in fam.items():
    floor = next(s for _, role, s in rungs if role == "floor")
    o = dict(rungs=[], tests=[])
    md += ["", f"## {f}", "", "| Rung | Role | Units | Mean score |", "|---|---|---|---|"]
    for lab, role, s in rungs:
        n = len({unit_of(c) for c in s})
        mean = float(stats.aggregate_units(s, {c: unit_of(c) for c in s})[1].mean()) if s else None
        o["rungs"].append(dict(label=lab, role=role, n_units=n, mean=mean))
        md.append(f"| {lab} | {role} | {n} | {'–' if mean is None else format(mean, '.4f')} |")
    md += ["", "| Comparison | Units | D | Bounds | Ahead |", "|---|---|---|---|---|"]
    for nlab, _, ns in [r for r in rungs if r[1] == "naive" and r[2]]:
        for lab, role, s in [r for r in rungs if r[1] != "naive"]:
            t = test(s, ns) | dict(first=lab, second=nlab, role=role)
            o["tests"].append(t)
            md.append(
                f"| {lab} − {nlab} | {t['n_units']} | {t['est']:+.4f} | {t['lower']:+.4f}, {t['upper']:+.4f} | "
                f"{'yes' if t['ahead'] else 'no'} |"
            )
    o["floor_ahead_of_every_naive"] = all(t["ahead"] for t in o["tests"] if t["role"] == "floor")
    o["learned_ahead_of_every_naive"] = {
        lab: all(t["ahead"] for t in o["tests"] if t["first"] == lab) for lab, role, _ in rungs if role == "learned"
    }
    out["families"][f] = o
md += ["", "## Assay sensitivity", "", "| Family | Floor ahead of every naive rung | Learned models ahead of every naive rung |", "|---|---|---|"]  # fmt: skip
for f, o in out["families"].items():
    le = "; ".join(f"{k}: {'yes' if v else 'no'}" for k, v in o["learned_ahead_of_every_naive"].items())
    md.append(f"| {f} | {'yes' if o['floor_ahead_of_every_naive'] else 'no'} | {le} |")
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "X_ladder.json").write_text(json.dumps(out, indent=1, default=float))
(OUT / "X_ladder.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
print("X_LADDER_DONE", flush=True)
