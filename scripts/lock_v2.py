"""Lock v2 (PLAN §18, user decisions 2026-11-02) -> manifests/locked_v2.json (locked.json = v1, tag p3-lock, kept).
Development units only for every choice (§11 rule: highest development mean; human configs from human dev units,
mouse from Allen animal 4). Changes vs v1: F1 floor pearson_none_pc50 (scale-robust, §18); scVI grid {nl10, nl30, nl50,
nl100} x nb; F2 floor ell400_lam0.8 (G2 lock v2), CellCharter scVI seeds 0-2; F3 floor per species
(T+ctx:ell25:pc10:a1000 / ctx:ell100:pc50:a1000), LightGBM re-gridded on the floor's features, NCEM added, GAT dropped
(descriptive). Margins unchanged. Reads runs/p3 (v1 outputs) and runs/p3v2 (pre-lock reruns).
Usage: python scripts/lock_v2.py [--write] -> prints the table; --write saves manifests/locked_v2.json"""
import copy, json, sys
from pathlib import Path
import numpy as np

V1 = json.loads(Path("manifests/locked.json").read_text()); L = copy.deepcopy(V1["locked"]); R2 = Path("runs/p3v2")
DEVH = ["gse1_232", "gse0_S0"] + [f"gse303070_148191A_{s}" for s in ("118N", "278N", "278T", "46T")]; DEVM = ["allen_animal4_C57BL6J-4.003"]
SP = {"human": DEVH, "mouse": DEVM}; miss = []


RS = R2 / "f1_rescored.csv"  # SPEC-exact library sizes (p3_f1_rescore.py), as the v1 grid means; raw json scores otherwise
RS = __import__("pandas").read_csv(RS) if RS.exists() else None


def f1(u, model, cfg):
    if RS is not None:
        s = RS[(RS.unit == u) & (RS.model == model) & (RS.config == cfg)].score
        if len(s): return float(s.iloc[0])
    p = R2 / f"f1_{u}_{model}_{cfg}.json"
    if not p.exists(): miss.append(str(p)); return np.nan
    return json.loads(p.read_text())["scores"][cfg]["score"]


def grid(fam, u, model, d=R2):
    p = d / f"{fam}_{u}_{model}.json"
    if not p.exists(): miss.append(str(p)); return {}
    return {k: v.get("unit", v.get("ari")) for k, v in json.loads(p.read_text())["grid"].items()}


def choose(means):
    ok = {k: v for k, v in means.items() if np.isfinite(v)}; return max(ok, key=ok.get) if ok else None


for s, U in SP.items():
    L["F1"]["expr_embed"][s] = dict(config="pearson_none_pc50", dev_mean=float(np.nanmean([f1(u, "expr_embed", "pearson_none_pc50") for u in U])), v1=V1["locked"]["F1"]["expr_embed"][s]["config"])
    gm = {k: V1["locked"]["F1"]["scvi"][s]["grid_means"][k] for k in ("nl10_nb", "nl30_nb")}
    gm.update({f"nl{n}_nb": float(np.nanmean([f1(u, "scvi", f"nl{n}_nb") for u in U])) for n in (50, 100)})
    c = choose(gm); L["F1"]["scvi"][s] = dict(config=c, dev_mean=gm[c], grid_means=gm)
    L["F2"]["floor"][s] = dict(config="ell400_lam0.8", note="G2 lock v2 (§17 / §18)", v1=V1["locked"]["F2"]["floor"][s]["config"])
    L["F2"]["cellcharter"][s] = dict(V1["locked"]["F2"]["cellcharter"][s], scvi_seeds=[0, 1, 2])
    fe = "T+ctx:ell25:pc10" if s == "human" else "ctx:ell100:pc50"
    rc = f"{fe}:a1000"; L["F3"]["ridge"][s] = dict(config=rc, dev_mean=float(np.nanmean([grid("f3", u, "ridge").get(rc, np.nan) for u in U])), v1=V1["locked"]["F3"]["ridge"][s]["config"])
    for m in ("lightgbm", "ncem"):
        gs = [grid("f3", u, m) for u in U]; keys = set.intersection(*(set(g) for g in gs)) if all(gs) else set()
        gm = {k: float(np.mean([g[k] for g in gs])) for k in keys}; c = choose(gm)
        L["F3"][m] = L["F3"].get(m, {}) if isinstance(L["F3"].get(m), dict) else {}
        L["F3"][m][s] = dict(config=c, dev_mean=gm.get(c, np.nan), grid_means=gm, features=fe)
    xg = [grid("f3", u, "lgbm_resid") for u in U]; L["F3"].setdefault("lgbm_resid", {})[s] = dict(config="resid_leaves15_lr0.05_es", secondary=True, features=fe,
                                                                                                   dev_mean=float(np.nanmean([g.get("resid_leaves15_lr0.05_es", np.nan) for g in xg])))
L["F3"]["gat"] = dict(dropped="PLAN §18 (user 2026-11-02): not a test pair; P3 / validation results descriptive", v1=V1["locked"]["F3"]["gat"])
out = dict(rule="PLAN §18 lock v2 (development units only; §11 highest dev mean for learned grids; floors per §14 / §18)", locked=L, v1_commit=V1["spatialfloors_commit"])
for f in ("F1", "F3"):
    for m, v in L[f].items():
        if isinstance(v, dict) and "human" in v: print(f, m, {s: (v[s].get("config"), round(v[s].get("dev_mean", np.nan), 4)) for s in ("human", "mouse") if s in v})
print("missing:", len(miss), miss[:6])
if "--write" in sys.argv:
    from spatialfloors.provenance import stamp
    Path("manifests/locked_v2.json").write_text(json.dumps(stamp(out), indent=1, default=float)); print("-> manifests/locked_v2.json")
