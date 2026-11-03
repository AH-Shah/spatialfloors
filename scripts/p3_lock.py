"""P3 lock (PLAN §6, §11): per family x model x species, the config with the highest mean development score (ties: first
in grid order; Novae: collapsed configs excluded); margins = stats.margin_from_dev on the locked floor configs.
Human configs from human dev units; mouse configs from Allen animal 4 when the model ran there, else human.
Refuses to lock if any expected P3 result is missing.
Floors (PLAN §14 choice rule, amended 2026-10-30): the floor config is the landscape plateau centre (FLOOR_LOCK), not the
argmax; it must be in the P3 grid. F2 margin on balanced ARI (row 359 decision) from runs/landscape/f2_<unit>.jsonl. Usage (env floors): python scripts/p3_lock.py -> manifests/locked.json"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from spatialfloors import stats
from spatialfloors.provenance import stamp

R = Path("runs/p3"); A4 = "allen_animal4_C57BL6J-4.003"
HUM = ["gse1_232", "gse0_S0", "gse303070_148191A_118N", "gse303070_148191A_278N", "gse303070_148191A_278T", "gse303070_148191A_46T"]
F5U = [u for u in HUM if u != "gse0_S0"]; F2H = ["gse1_232"]
FAM = dict(F1=dict(floor="expr_embed", learned=["scvi", "nicheformer"]), F2=dict(floor="floor", learned=["cellcharter", "novae"]),
           F3=dict(floor="ridge", learned=["gat", "lightgbm"]), F5=dict(floor="floors", learned=["scanvi", "tangram"]))
FLOOR_LOCK = dict(F1="pearson_z_pc100", F2="ell200_lam0.8", F3="ell25_pc50", F5="insitutype")  # SPEC; landscape_choose.py, f2_choose.py
F2_LAND = "grid|gmm_full|ell200|lam0.8|pc50"
missing, rows = [], []


def grid(fam, u, m):
    f = R / f"{fam.lower()}_{u}_{m}.json"
    if not f.exists(): missing.append(str(f)); return {}
    g = json.loads(f.read_text())["grid"]; key = dict(F2="ari", F3="unit", F5="heldout_auroc")[fam]
    return {c: v[key] for c, v in g.items() if not v.get("collapse", False) and not c.startswith("gatv2")}  # GATv2 dropped 2026-10-28


f1 = pd.read_csv(R / "f1_rescored.csv") if (R / "f1_rescored.csv").exists() else None
if f1 is None: missing.append("runs/p3/f1_rescored.csv")
for fam, d in FAM.items():
    for m in [d["floor"]] + d["learned"]:
        for sp, units in (("human", F5U if fam == "F5" else F2H if fam == "F2" else HUM), ("mouse", [A4] if fam != "F5" else [])):
            for u in units:
                if fam == "F1":
                    if f1 is None: continue
                    s = f1[(f1.unit == u) & (f1.model == m)]
                    if s.empty: missing.append(f"F1 {u} {m}")
                    for _, r in s.iterrows(): rows.append(dict(fam=fam, model=m, species=sp, unit=u, config=r.config, score=r.score, chance=r.chance))
                else:
                    for c, v in grid(fam, u, m).items(): rows.append(dict(fam=fam, model=m, species=sp, unit=u, config=c, score=v, chance=.5 if fam == "F5" else 0.))
if missing: print("MISSING", *missing, sep="\n"); sys.exit(2)
df = pd.DataFrame(rows); out = dict(rule="highest mean dev score; ties first in grid order; mouse from Allen animal 4 when run", locked={}, dev_table=df.to_dict("records"))
for (fam, m, sp), g in df.groupby(["fam", "model", "species"], sort=False):
    mean = g.groupby("config", sort=False).score.mean(); best = mean.idxmax()
    if m == FAM[fam]["floor"]: assert FLOOR_LOCK[fam] in mean.index, (fam, sp, FLOOR_LOCK[fam]); best = FLOOR_LOCK[fam]
    out["locked"].setdefault(fam, {}).setdefault(m, {})[sp] = dict(config=best, dev_mean=float(mean[best]), grid_means=mean.round(4).to_dict())
for fam, d in FAM.items():
    fl = out["locked"][fam][d["floor"]]; cfg = fl["human"]["config"]; g = df[(df.fam == fam) & (df.model == d["floor"]) & (df.config == cfg) & (df.species == "human")]
    if fam == "F2":  # balanced ARI at the locked config on every development unit
        sc = [next(r["ari_bal"] for r in map(json.loads, open(f"runs/landscape/f2_{u}.jsonl")) if r["id"] == F2_LAND) for u in F2H + [A4]]
        out["locked"][fam]["margin"] = stats.margin_from_dev(fam, sc, np.zeros(len(sc))); out["locked"][fam]["primary"] = "ari_bal"; continue
    out["locked"][fam]["margin"] = stats.margin_from_dev(fam, g.score.values, g.chance.values)
out["locked"]["F4"] = dict(margin=stats.margin_from_dev("F4", [0], [0]), note="authors' configs; L6 locked (SPEC F4)")
out["equivalence"] = ("per (family, learned model) pair, paired per-unit D = learned - floor, margin = locked margin: floor non-inferior "
    "if the one-sided 95% upper bound of mean D < margin (pre-registered primary, PLAN §5); equivalent if the 90% CI of mean D lies in "
    "(-margin, +margin) (TOST, reported); floor superior if the upper bound < 0; learned superior if the lower bound > 0 (Holm)")
Path("manifests/locked.json").write_text(json.dumps(stamp(out), indent=1, default=float))
print(json.dumps({f: {m: (v if not isinstance(v, dict) else {s: x["config"] for s, x in v.items()}) for m, v in d.items()} for f, d in out["locked"].items()}, indent=1))
