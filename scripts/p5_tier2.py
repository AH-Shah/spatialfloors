"""Tier-2 diagnostic (PROTOCOL §5), run once after the tier-1 report (<runs>/P5.json) exists; it never changes a tier-1
verdict. Grids, locked configs and the unit subset: manifests/tier2.json. Per model and crop: best = highest primary
over
the model's grid (ties -> grid order), transfer gap T = best - locked (>= 0). Per (family, learned model) pair and crop:
capacity gap C = best learned - best floor = D + T_learned - T_floor (D = the tier-1 difference on the same crop). Crops
-> units as in tier 1 (stats.aggregate_units); unit means with one-sided 95% bounds from stats.family_test, descriptive
(selection on the scored unit, no multiplicity control). Models marked units "subset" in the manifest enter on the
subset
only. Reading of T: upper bound < margin -> transfers; lower bound > margin -> under-represented; otherwise
inconclusive.
Outputs are read from every run directory under <runs> (tier 1 at its root, tier 2 under <runs>/t2).
Usage (env floors): python scripts/p5_tier2.py <runs dir> <crop ...> -> <runs dir>/T2.json + T2.md"""

import glob, hashlib, json, os, sys
from pathlib import Path
import numpy as np, pandas as pd
from spatialfloors import stats

R = Path(sys.argv[1])
U = sys.argv[2:]
T1 = R / os.environ.get("SF_T1NAME", "P5.json")
assert T1.exists(), f"tier 2 runs only after the tier-1 report ({T1}) is written"
D2 = json.loads(Path("manifests/tier2.json").read_text())
G = D2["grids"]
SUB = (
    set(U) if os.environ.get("SF_T2_SUBSET") == "all" else set(D2["subset_crops"])
)  # "all": dry runs on development units
LK = json.loads(Path(D2["locked_file"]).read_text())["locked"]
MARGIN = {f: LK[f]["margin"] for f in G}
sp = lambda u: "mouse" if u.startswith("allen") else "human"
donors = {
    f"gse303070_{d['slide'].replace('Tumor', 'Tumour')}_{d['pref']}": d["donor_cells"]
    for d in json.loads(Path("manifests/gse303070_donors.json").read_text())["check_vs_fov_csv"]
}
FILES = {
    "F2": {"floor": ["floor"], "cellcharter": ["cellcharter"], "novae": ["novae"]},
    "F3": {"floor": ["ridge"], "lightgbm": ["lightgbm"], "ncem": ["ncemA"]},
    "F5": {"floor": ["floors", "prior"], "scanvi": ["scanvi"], "tangram": ["tangram"]},
}
F1M = {"floor": "expr_embed", "scvi": "scvi", "nicheformer": "nicheformer"}
depth = lambda p: (len(Path(p).parts), p)  # tier 1 (shallowest) first


def f1_table():
    fs = sorted(glob.glob(str(R / "**" / "f1_rescored.csv"), recursive=True), key=depth)
    return (
        pd.concat([pd.read_csv(f) for f in fs]).drop_duplicates(["unit", "model", "config"], keep="first")
        if fs
        else None
    )


F1 = f1_table()
_f2 = {}


def f2_mask(u):
    if u not in _f2:
        import anndata as ad
        from spatialfloors.metrics.f2 import scored_regions

        a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad", backed="r")
        _f2[u] = scored_regions(
            a.obs.region.astype(str).values, json.loads(Path("manifests/crops/index.json").read_text())[u]["present"]
        )
        a.file.close()
    return _f2[u]


def scores(f, m, u):
    """{config: primary} of model m on crop u over every output under R; the first file (shallowest) wins a
    duplicate."""
    out, g, alias = {}, G[f][m], {v[sp(u)]: k for k, v in G[f][m].get("alias", {}).items()}
    if f == "F1":
        if F1 is not None:
            for r in F1[(F1.unit == u) & (F1.model == F1M[m])].itertuples():
                out.setdefault(r.config, float(r.score))
    elif f == "F2":
        from spatialfloors.metrics.f2 import ari_bal_saved

        for fn in FILES[f][m]:
            for p in sorted(glob.glob(str(R / "**" / "labels" / f"f2_{u}_{fn}.npz"), recursive=True), key=depth):
                for k in np.load(p, allow_pickle=True).files:
                    if k in g["configs"] and k not in out:
                        mask, reg = f2_mask(u)
                        out[k] = ari_bal_saved(p, k, mask, reg)
    else:
        key = "unit" if f == "F3" else "heldout_auroc"
        for fn in FILES[f][m]:
            for p in sorted(glob.glob(str(R / "**" / f"{f.lower()}_{u}_{fn}.json"), recursive=True), key=depth):
                for k, v in json.loads(Path(p).read_text()).get("grid", {}).items():
                    out.setdefault(alias.get(k, k), float(v[key]))
    return {c: v for c, v in out.items() if c in g["configs"]}


rec = {}  # (family, model) -> crop -> record
for f, ms in G.items():
    for m, g in ms.items():
        for u in U:
            lk = g["locked"].get(sp(u))
            if lk is None or (g["units"] == "subset" and u not in SUB):
                continue
            s = scores(f, m, u)
            if lk not in s:
                continue
            best = max(s, key=lambda c: (s[c], -g["configs"].index(c)))
            rec.setdefault((f, m), {})[u] = dict(
                locked=s[lk],
                best=s[best],
                best_config=best,
                T=s[best] - s[lk],
                configs=s,
                complete=len(s) == len(g["configs"]),
            )


def summary(vals, margin):
    units, D, strata = stats.aggregate_units(vals, {c: donors.get(c, c) for c in vals})
    t = stats.family_test(D, strata, margin=margin) if len(D) > 1 else dict(est=float(D[0]), upper=np.nan, lower=np.nan)
    return dict(
        n_units=len(D), est=t["est"], lower=t["lower"], upper=t["upper"], per_unit=dict(zip(units, map(float, D)))
    )


out = dict(
    runs=str(R),
    tier1=str(T1),
    tier1_sha256=hashlib.sha256(T1.read_bytes()).hexdigest(),
    tier2_manifest_sha256=hashlib.sha256(Path("manifests/tier2.json").read_bytes()).hexdigest(),
    transfer={},
    capacity={},
)
for (f, m), r in rec.items():
    s = summary({u: v["T"] for u, v in r.items()}, MARGIN[f])
    s.update(
        reading="not estimable"
        if not np.isfinite(s["upper"])
        else "transfers"
        if s["upper"] < MARGIN[f]
        else "under-represented by development tuning"
        if s["lower"] > MARGIN[f]
        else "inconclusive",
        margin=MARGIN[f],
        locked_is_best=float(np.mean([v["best"] == v["locked"] for v in r.values()])),
        max_T=max(v["T"] for v in r.values()),
        incomplete=[u for u, v in r.items() if not v["complete"]],
        best_configs={u: v["best_config"] for u, v in r.items()},
    )
    out["transfer"][f"{f}|{m}"] = s
for f, ms in G.items():
    fl = rec.get((f, "floor"), {})
    for m in ms:
        if m == "floor" or (f, m) not in rec:
            continue
        both = sorted(set(fl) & set(rec[(f, m)]))
        if not both:
            continue
        C = summary({u: rec[(f, m)][u]["best"] - fl[u]["best"] for u in both}, MARGIN[f])
        D = summary({u: rec[(f, m)][u]["locked"] - fl[u]["locked"] for u in both}, MARGIN[f])
        C.update(
            reading="not estimable"
            if not np.isfinite(C["upper"])
            else "floor holds at best tuning"
            if C["upper"] < MARGIN[f]
            else "learned ahead at best tuning"
            if C["lower"] > 0
            else "inconclusive",
            margin=MARGIN[f],
            D_same_units=D["est"],
            T_learned=float(np.mean([rec[(f, m)][u]["T"] for u in both])),
            T_floor=float(np.mean([fl[u]["T"] for u in both])),
        )
        out["capacity"][f"{f}|{m}"] = C
NAME = os.environ.get("SF_T2NAME", "T2")
(R / f"{NAME}.json").write_text(json.dumps(out, indent=1, default=float))
md = ["| model | units | T mean | upper 95% | margin | locked = best | reading |", "|---|---|---|---|---|---|---|"]
md += [
    f"| {k} | {v['n_units']} | {v['est']:+.4f} | {v['upper']:+.4f} | {v['margin']:.4f} | {v['locked_is_best']:.2f} | {v['reading']} |"
    for k, v in out["transfer"].items()
]
md += [
    "",
    "| pair | units | C mean | lower 95% | upper 95% | D same units | T learned | T floor | reading |",
    "|---|---|---|---|---|---|---|---|---|",
]
md += [
    f"| {k} | {v['n_units']} | {v['est']:+.4f} | {v['lower']:+.4f} | {v['upper']:+.4f} | {v['D_same_units']:+.4f} | {v['T_learned']:+.4f} | {v['T_floor']:+.4f} | {v['reading']} |"
    for k, v in out["capacity"].items()
]
(R / f"{NAME}.md").write_text("\n".join(md) + f"\n\ntier-1 report {T1} sha256 {out['tier1_sha256']}\n")
print(len(out["transfer"]), "models,", len(out["capacity"]), "pairs ->", R / f"{NAME}.md")
