"""PLAN §19 lock: manifests/locked_v3.json = locked_v2.json + the v3 outcomes chosen on development (runs/v3):
  F1  v2 kept (no candidate passed); F2 v2 kept, best candidate run as a reported secondary floor; F3 floor v3 (lock +
  same-lineage rings, alpha 1000) with v2 reported, LightGBM on the v3 feature groups (its 4-config grid re-chosen on
  development), NCEM = the authors' models (rows 378-380); F4 L6 kept; F5 floor v3 = InSituType + spatial prior.
GAT stays dropped (§18). Margins unchanged. Usage: python3 scripts/v3_lock.py"""
import glob, json, subprocess
from pathlib import Path
import numpy as np

V2 = json.loads(Path("manifests/locked_v2.json").read_text()); L = json.loads(json.dumps(V2["locked"])); J = lambda p: json.loads(Path(p).read_text())
sp = lambda u: "mouse" if u.startswith("allen") else "human"; DEV = ["gse1_232", "gse0_S0", "gse303070_148191A_118N", "gse303070_148191A_278N", "gse303070_148191A_278T", "gse303070_148191A_46T", "allen_animal4_C57BL6J-4.003"]
alt = {u: J(f"runs/v3/f3alt_{u}.json")["unit_mean"] for u in DEV}
for s in ("human", "mouse"):
    us = [u for u in DEV if sp(u) == s]
    L["F3"]["ridge"][s] = dict(config="v3", v2=V2["locked"]["F3"]["ridge"][s]["config"], design="v2 lock features + same-lineage rings 50-150 / 150-250 um (receiver's masked genes, 10 SVD), ridge alpha 1000",
                               dev_mean=float(np.mean([alt[u]["fixed_sf"] for u in us])), dev_mean_v2=float(np.mean([alt[u]["lock"] for u in us])),
                               dev_null=float(np.mean([alt[u]["fixed_sf_null"] for u in us])), note="PLAN §19; scripts/v3_f3_alt.py 'fixed_sf'")
    G = [J(f"runs/v3/lgbm/f3_{u}_lightgbm.json")["grid"] for u in us]; means = {k: float(np.mean([g[k]["unit"] for g in G])) for k in G[0]}; best = max(means, key=means.get)
    L["F3"]["lightgbm"][s] = dict(config=best, features="v3", dev_mean=means[best], grid_means=means, v2=V2["locked"]["F3"]["lightgbm"][s]["config"], note="v3 feature groups (floors.f3.v3_parts, a superset of the floor's)")
L["F3"]["ncemA"] = {"human": dict(config="edA_e1_r19,linA_l20_r19", note="authors' ED-NCEM + linear NCEM, radius 19 um (row 379; validation row 380)"),
                    "mouse": dict(config="edA_e1_r46,linA_l20_r46", note="radius 46 um (row 379)")}
f5 = {u: J(f"runs/v3/f5_{u}.json")["summaries"] for u in DEV if sp(u) == "human" and u != "gse0_S0"}  # S0: no F5 (no IF, not an F5 unit)
L["F5"]["floors"]["human"] = dict(config="insitutype_spatial", v2="insitutype", dev_mean=float(np.mean([f5[u]["v2+spatial"]["heldout_auroc"] for u in f5])),
                                   dev_mean_v2=float(np.mean([f5[u]["v2"]["heldout_auroc"] for u in f5])), note="PLAN §19: InSituType ML + mean-field neighbour prior (beta 1, 15-NN) on its log-likelihoods")
for s in ("human", "mouse"): L["F2"]["floor"][s]["secondary_v3"] = "v3:ms50_200_800+c200"; L["F2"]["floor"][s]["v3_note"] = "PLAN §19: no candidate passed (region-recall rule); best run as a reported secondary floor"
L["F1"]["expr_embed"]["v3_note"] = "PLAN §19: v2 kept (no candidate passed); depth-matched reference reported"
L["F4"]["v3_note"] = "PLAN §19: L6 kept (L7 15/30/100 um, similarity-weighted, EB, rank 1 did not beat L6 on crc_232 whole slide or crop)"
commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
out = dict(rule="PLAN §19 (2026-11-03): v3 floors chosen on development units only; locked before p4-lock", locked=L, v2_file="manifests/locked_v2.json", spatialfloors_commit=commit)
Path("manifests/locked_v3.json").write_text(json.dumps(out, indent=1)); print("locked_v3:", {s: (L["F3"]["ridge"][s]["dev_mean"], L["F3"]["lightgbm"][s]["config"]) for s in ("human", "mouse")}, L["F5"]["floors"]["human"]["config"])
