"""Tier-2 design (PROTOCOL §5) -> manifests/tier2.json: each model's grid, its locked config, the configs tier 1 already
produces, and the seeded dataset-stratified unit subset for the GPU-trained models and NCEM. Written before any confirmation
result exists; p4_build.py builds the tier-2 jobs from it and p5_tier2.py analyses them.
Subset: per dataset ceil(n / 2) units drawn without replacement by numpy default_rng(seed) from the sorted unit list
(datasets in sorted order); seed = first 8 hex digits of SHA-256(manifests/locked_v3.json).
Usage: python scripts/t2_manifest.py"""
import hashlib, json, math
from pathlib import Path
import numpy as np

LOCK = Path("manifests/locked_v3.json"); L = json.loads(LOCK.read_text())["locked"]
idx = json.loads(Path("manifests/crops/index.json").read_text())
donor = {f"gse303070_{d['slide'].replace('Tumor', 'Tumour')}_{d['pref']}": d["donor_cells"] for d in json.loads(Path("manifests/gse303070_donors.json").read_text())["check_vs_fov_csv"]}
ALLEN = ["allen_animal1_C57BL6J-1.086", "allen_animal2_C57BL6J-2.039", "allen_animal3_C57BL6J-3.003", "allen_yao638850_C57BL6J-638850.36"]
crops = sorted(u for u, v in idx.items() if isinstance(v, dict) and v.get("role") in ("confirmation", "confirmation_F1")) + ALLEN
unit = lambda c: donor.get(c, c.split("_C57BL6J")[0] if c.startswith("allen") else c)
dataset = lambda c: {"gse1": "GSE1", "gse303070": "GSE303070", "bruker": "Bruker WTx", "spatch": "SPATCH", "allen": "Allen MERFISH"}[c.split("_")[0]]
units = {}
for c in crops: units.setdefault(dataset(c), {}).setdefault(unit(c), []).append(c)

sha = hashlib.sha256(LOCK.read_bytes()).hexdigest(); seed = int(sha[:8], 16); rng = np.random.default_rng(seed); subset = {}
for d in sorted(units):
    u = sorted(units[d]); subset[d] = sorted(rng.choice(u, math.ceil(len(u) / 2), replace=False).tolist())
sub_crops = sorted(c for d in subset for u in subset[d] for c in units[d][u])

sp = ("human", "mouse"); lk = lambda f, m: {s: L[f][m][s]["config"] for s in sp if s in L[f][m]}
NCEM = {s: v["config"] for s, v in json.loads(Path("manifests/locked_v3_ncem.json").read_text())["locked"].items()}  # PROTOCOL §2 (scripts/lock_ncem.py)
G = {
 "F1": {"floor": dict(configs=["pearson_none_pc50", "pearson_none_pc100", "pearson_z_pc50", "pearson_z_pc100"], locked=lk("F1", "expr_embed"), tier1={}, units="all"),
        "scvi": dict(configs=["nl10_nb", "nl30_nb", "nl50_nb", "nl100_nb"], locked=lk("F1", "scvi"), tier1={}, units="subset"),
        "nicheformer": dict(configs=["final", "penultimate"], locked=lk("F1", "nicheformer"), tier1={s: ["penultimate"] for s in sp}, units="all")},
 "F2": {"floor": dict(configs=["ell400_lam0.8", "ell400_lam1.0", "ell200_lam0.8", "v3:ms50_200_800+c200"], locked=lk("F2", "floor"),
                      tier1={s: ["v3:ms50_200_800+c200"] for s in sp}, units="all"),
        "cellcharter": dict(configs=["nl10_L1", "nl10_L3", "nl30_L1", "nl30_L3"], locked=lk("F2", "cellcharter"), tier1={}, units="subset"),
        "novae": dict(configs=["zero_shot", "fine_tune"], locked=lk("F2", "novae"), tier1={}, units="all")},  # 4 mouse units: a subset of 2 leaves df 1
 "F3": {"floor": dict(configs=["v3_a300", "v3", "v3_a3000", "v3_a10000"], locked=lk("F3", "ridge"), tier1={}, units="all"),  # v3 = alpha 1000; the v3 development alpha grid
        "lightgbm": dict(configs=["leaves15_lr0.05", "leaves15_lr0.1", "leaves63_lr0.05", "leaves63_lr0.1"], locked=lk("F3", "lightgbm"), tier1={}, units="all"),
        "ncem": dict(configs=["edA_e1_r19", "linA_l20_r19", "edA_e1_r46", "linA_l20_r46"], locked=NCEM, tier1={"human": ["edA_e1_r19"], "mouse": ["edA_e1_r46"]}, units="subset")},
 "F5": {"floor": dict(configs=["insitutype", "insitutype_spatial", "prior_b0.5", "prior_b2"], locked=lk("F5", "floors"), tier1={"human": ["insitutype"]}, units="all"),
        "scanvi": dict(configs=["nl10_nb", "nl10_zinb", "nl30_nb", "nl30_zinb"], locked=lk("F5", "scanvi"), tier1={}, units="subset"),
        "tangram": dict(configs=["uniform", "rna_count_based"], locked=lk("F5", "tangram"), tier1={}, units="subset")}}
for f, ms in G.items():
    for m, g in ms.items():
        for s, c in g["locked"].items(): assert c in g["configs"], (f, m, s, c)
out = dict(rule="PROTOCOL §5 (tier-2 diagnostic); computed only after the tier-1 report (runs/p4/P5.md) is written", locked_file=str(LOCK), locked_sha256=sha,
           seed=seed, seed_source="int(SHA-256(manifests/locked_v3.json)[:8], 16)", subset_rule="per dataset ceil(n/2) units, numpy default_rng(seed).choice on the sorted unit list, datasets in sorted order",
           units={d: {u: units[d][u] for u in sorted(units[d])} for d in sorted(units)}, subset=subset, subset_crops=sub_crops,
           families={"F1": "all F1 crops", "F2": "gse1 + allen", "F3": "gse1 + gse303070 + allen + bruker", "F5": "gse1 + gse303070 + bruker"},
           excluded={"F4": "no learned-model grid (authors' configurations only) and the floor is locked on both sides"}, grids=G)
Path("manifests/tier2.json").write_text(json.dumps(out, indent=1) + "\n")
print("seed", seed, "| subset", {d: len(v) for d, v in subset.items()}, "| crops", len(sub_crops), "-> manifests/tier2.json")
