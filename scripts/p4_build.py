"""Confirmation jobs (PROTOCOL §4-§6) -> runs/sched/p4_jobs.txt; scripts/swap_p4_block.py puts them in the pool's job file.
Locked configs from $SF_LOCK (manifests/locked_v3.json); every job runs through scripts/blind.sh (sealed logs), outputs
runs/p4. Gates (pseudo-jobs that fail until runs/sched/<gate>.ok exists and the .fail is removed):
  G7      learned models (configs frozen)
  P4LOCK  floors and their F4 L6 runs
  T2      tier 2 (manifests/tier2.json): only once the tier-1 report runs/p4/P5.md exists
Crop caches and F4 input files run at once (no metrics). F4 MERFISH: prospective Allen animals 1, 3, Yao (3 sections each,
Cellina / CPA x 5 types x 3 seeds, L6); retrospective animal 2: L6 only (authors' Cellina / CPA numbers).
Usage: SF_LOCK=manifests/locked_v3.json python scripts/p4_build.py --bruker"""
import json, sys
from pathlib import Path

import os
LOCK = os.environ.get("SF_LOCK", "manifests/locked_v3.json")  # earlier locks: manifests/locked.json (v1), locked_v2.json
L = json.loads(Path(LOCK).read_text())["locked"]; idx = json.loads(Path("manifests/crops/index.json").read_text()); V3 = "v3" in LOCK
GSE1 = [f"gse1_{s}" for s in (110, 120, 210, 221, 222, 242)]
G303 = [u for u, v in idx.items() if isinstance(v, dict) and v.get("role") == "confirmation" and u.startswith("gse303070_")]
ALLEN = ["allen_animal1_C57BL6J-1.086", "allen_animal2_C57BL6J-2.039", "allen_animal3_C57BL6J-3.003", "allen_yao638850_C57BL6J-638850.36"]
XEN = ["spatch_xenium5k"]; BRUKER = ["bruker_mobreast", "bruker_skin", "bruker_kidney"] if "--bruker" in sys.argv else []
F1U, F2U, F3U, F5U = GSE1 + G303 + XEN + ALLEN + BRUKER, GSE1 + ALLEN, GSE1 + G303 + ALLEN + BRUKER, GSE1 + G303 + BRUKER
sp = lambda u: "mouse" if u.startswith("allen") else "human"; cfg = lambda f, m, u: L[f][m][sp(u)]["config"]
short = lambda u: u.replace("gse303070_", "g3_").replace("allen_", "a_").split("_C57BL6J")[0]
rows = ["# P4 confirmation (scripts/p4_build.py; PLAN §17). Blind: sealed logs in runs/p4/sealed.",
        "G7|cpu|.1|1|runs/sched/G7.ok|-|test -e runs/sched/G7.ok", "P4LOCK|cpu|.1|1|runs/sched/P4LOCK.ok|-|test -e runs/sched/P4LOCK.ok"]


def job(id_, cls, est, cap, done, after, cmd, env=""):
    rows.append(f"{id_}|{cls}|{est}|{cap}|{done}|{after}|zsh scripts/blind.sh {id_} env SF_OUT=runs/p4 {env} {cmd}".replace("  ", " "))


for u in sorted(set(F1U + F2U + F3U + F5U)):
    if not Path(f"runs/p2/cache/{u}.h5ad").exists() and not u.startswith("bruker"):
        rows.append(f"c_{short(u)}|cpu|4|14|runs/p2/cache/{u}.h5ad|-|pixi run -e floors python scripts/p2_cache_crop.py {u}")
dep = lambda u, g: ",".join([g] + ([f"c_{short(u)}"] if any(r.startswith(f"c_{short(u)}|") for r in rows) else []))
P = "runs/p4"
for u in F1U:  # learned, GPU, after G7
    job(f"p4_1s_{short(u)}", "gpu", 2, 12, f"{P}/f1_{u}_scvi_{cfg('F1', 'scvi', u)}.json", dep(u, "G7"), f"pixi run -e scvi python scripts/f1_fit.py {u} scvi {cfg('F1', 'scvi', u)}")
    job(f"p4_1n_{short(u)}", "gpu", 4, 14, f"{P}/f1_{u}_nicheformer_layers.json", dep(u, "G7"), f"pixi run -e nicheformer python scripts/f1_fit.py {u} nicheformer layers")
for u in F5U:
    for m, e in (("scanvi", 4), ("tangram", 5)):
        job(f"p4_5{m[0]}_{short(u)}", "gpu", e, 12, f"{P}/f5_{u}_{m}.json", dep(u, "G7"), f"pixi run -e scvi python scripts/p3_f5.py {u} {m}", f"SF_CONFIG={cfg('F5', m, u)}")
for u in F2U:
    job(f"p4_2c_{short(u)}", "gpu", 4, 12, f"{P}/f2_{u}_cellcharter.json", dep(u, "G7"), f"pixi run -e cellcharter python scripts/p3_f2.py {u} cellcharter",
        f"SF_CONFIG={cfg('F2', 'cellcharter', u)}" + (" SF_SCVI_SEEDS=0,1,2" if V3 else ""))  # §18: scVI seeds 0-2 averaged
    job(f"p4_2n_{short(u)}", "cpu", 7, 14, f"{P}/f2_{u}_novae.json", dep(u, "G7"), f"pixi run -e novae python scripts/p3_f2.py {u} novae", "SF_CONFIG=zero_shot")
if V3:  # F3 receiver ids frozen once (env floors, label-free) so NCEM (ncem_env) and LightGBM (scvi) share the floor's receivers
    F3NEW = [u for u in F3U if not Path(f"runs/p2/f3_ids/{u}.npy").exists()]
    cdeps = ",".join(["G7"] + [r.split("|")[0] for r in rows if r.startswith("c_") and any(r.startswith(f"c_{short(u)}|") for u in F3NEW)])
    rows.append(f"p4_3ids|cpu|4|10|runs/p2/f3_ids/{F3NEW[-1]}.npy|{cdeps}|env OMP_NUM_THREADS=4 pixi run -e floors python scripts/f3_ids_dump.py {' '.join(F3NEW)}")
NT = "PYTHONPATH=src OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 TF_NUM_INTRAOP_THREADS=4 TF_NUM_INTEROP_THREADS=2 TF_CPP_MIN_LOG_LEVEL=2"
for u in F3U:
    if "human" in L["F3"].get("gat", {}):  # GAT dropped as a test pair (PLAN §18): no P4 job under the v2 / v3 locks
        job(f"p4_3g_{short(u)}", "gpu", 4, 12, f"{P}/f3_{u}_gat.json", dep(u, "G7"), f"pixi run -e scvi python scripts/p3_f3.py {u} gat", f"SF_CONFIG={cfg('F3', 'gat', u)}")
    job(f"p4_3l_{short(u)}", "cpu", 2, 10, f"{P}/f3_{u}_lightgbm.json", dep(u, "G7") + (",p4_3ids" if V3 else ""), f"pixi run -e scvi python scripts/p3_f3.py {u} lightgbm",
        f"SF_CONFIG={cfg('F3', 'lightgbm', u)}" + (" SF_F3_FEATS=v3" if L["F3"]["lightgbm"][sp(u)].get("features") == "v3" else ""))
    if "ncemA" in L["F3"]:  # the authors' NCEMs (PLAN §18 / rows 378-380), both configs in one job
        job(f"p4_3n_{short(u)}", "cpu", 3, 10, f"runs/sched/p4_3n_{short(u)}.done", dep(u, "G7") + ",p4_3ids", f"external/ncem_env/bin/python scripts/p3_f3_ncem_auth.py {u}",  # the json is written per config: done = the pool's own marker
            f"{NT} SF_CONFIG={cfg('F3', 'ncemA', u)}")
# F4 GSE1 prospective (110 crop rule, 222 whole slide), 3 seeds x 5 types; L6 on every GSE1 confirmation slide
R4 = "SF_F4_RUN=runs/p4/f4"
rows.append(f"p4_4crop_110|cpu|6|16|runs/p4/f4/raw/cosmx_wt_110.h5ad|G7|env {R4} pixi run -e floors python scripts/f4_crop.py 110")
for sid in (110, 222):
    for s in (0, 1, 2):
        for h in ("Endothelial", "Epithelial", "Fibroblast", "Myeloid", "T_cell"):
            aft = "G7" + (",p4_4crop_110" if sid == 110 else ""); hd = "TVA" if sid == 222 else "CRC"  # SPEC F4: 222 holds out TVA
            job(f"p4_4c_{sid}_{h}_s{s}", "gpu", 8, 16, f"{P}/f4/datasets/crc/correlations/crc_{sid}_cellina_s{s}-cf_{h}_{hd}.json", aft,
                f"pixi run -e cellina python scripts/f4_run.py cellina {sid} {s} {h}", f"{R4} SF_HOLDOUT={hd} SF_ACCEL=mps")
            job(f"p4_4p_{sid}_{h}_s{s}", "cpu", 6, 16, f"{P}/f4/datasets/crc/correlations/crc_{sid}_cpa_s{s}-cf_{h}_{hd}.json", aft,
                f"pixi run -e cpa python scripts/f4_run.py cpa {sid} {s} {h}", f"{R4} SF_HOLDOUT={hd} SF_ACCEL=cpu SF_THREADS=4")
# F4 MERFISH (per section, as the authors' camera-ready)
AL = json.loads(Path("manifests/allen.json").read_text()); MT = ["GABAergic_neuron", "astrocyte", "endothelial_cell", "glutamatergic_neuron", "oligodendrocyte"]
for an in ("animal1", "animal3", "yao638850", "animal2"):
    secs = AL[an]["f4_sections"]; rows.append(f"p4_4mprep_{an}|cpu|8|16|runs/p4/f4/raw/{secs[-1]}.h5ad|-|env {R4} pixi run -e floors python scripts/f4_merfish_prep.py {an}")
    for sec in secs:
        if an != "animal2":
            for s in (0, 1, 2):
                for h in MT:
                    hs = h.replace("_", " ")
                    job(f"p4_4c_{sec}_{h}_s{s}", "gpu", 6, 14, f"{P}/f4/datasets/merfish/correlations/{sec}_cellina_s{s}-cf_{hs}_Fiber_tracts.json", f"G7,p4_4mprep_{an}",
                        f"pixi run -e cellina python scripts/f4_run.py cellina {sec} {s} {h}", f"{R4} SF_ACCEL=mps")
                    job(f"p4_4p_{sec}_{h}_s{s}", "cpu", 5, 14, f"{P}/f4/datasets/merfish/correlations/{sec}_cpa_s{s}-cf_{hs}_Fiber_tracts.json", f"G7,p4_4mprep_{an}",
                        f"pixi run -e cpa python scripts/f4_run.py cpa {sec} {s} {h}", f"{R4} SF_ACCEL=cpu SF_THREADS=4")
        job(f"p4_4l_{sec}", "cpu", 6, 14, f"{P}/f4/l6/{sec}.csv", f"P4LOCK,p4_4mprep_{an}", f"pixi run -e floors python scripts/f4_l6.py {sec}", R4)
# floors: after P4LOCK (the validation gate may re-centre them; SF_CONFIG read from locked.json at build time -> rebuild after a lock v2)
for u in F1U:
    job(f"p4_1f_{short(u)}", "cpu", 4, 10, f"{P}/f1_{u}_expr_embed_{cfg('F1', 'expr_embed', u)}.json", dep(u, "P4LOCK"), f"pixi run -e floors python scripts/f1_fit.py {u} expr_embed {cfg('F1', 'expr_embed', u)}")
for u in F2U:
    job(f"p4_2f_{short(u)}", "cpu", 1.5, 10, f"{P}/f2_{u}_floor.json", dep(u, "P4LOCK"), f"pixi run -e floors python scripts/p3_f2.py {u} floor", f"OMP_NUM_THREADS=4 SF_CONFIG={cfg('F2', 'floor', u)}")
for u in F2U:  # CellCharter on PCA input: published floor variant (PLAN §3), reported, not a test pair
    job(f"p4_2cp_{short(u)}", "cpu", 3, 10, f"{P}/f2_{u}_cellcharter_pca.json", dep(u, "P4LOCK"), f"pixi run -e cellcharter python scripts/p3_f2.py {u} cellcharter_pca", "OMP_NUM_THREADS=4 SF_CONFIG=pc50_L3")
for u in F3U:
    c3 = cfg('F3', 'ridge', u)
    job(f"p4_3r_{short(u)}", "cpu", 2, 10, f"{P}/f3_{u}_ridge.json", dep(u, "P4LOCK") + (",p4_3ids" if c3 == "v3" else ""), f"pixi run -e floors python scripts/p3_f3.py {u} ridge",
        f"OMP_NUM_THREADS=4 SF_CONFIG={c3}" + (" SF_F3_FEATS=v3" if c3 == "v3" else ""))  # v3 reads the frozen receiver ids
    if c3 == "v3":  # the v2 floor, reported next to v3 (PLAN §19)
        job(f"p4_3r2_{short(u)}", "cpu", 2, 10, f"{P}/v2/f3_{u}_ridge.json", dep(u, "P4LOCK"), f"pixi run -e floors python scripts/p3_f3.py {u} ridge", f"OMP_NUM_THREADS=4 SF_OUT={P}/v2 SF_CONFIG={L['F3']['ridge'][sp(u)]['v2']}")
if V3:  # F2 v3 candidate (failed one §19 rule on development): run and reported as a secondary floor, not a test pair
    for u in F2U:
        c2 = L["F2"]["floor"][sp(u)]["secondary_v3"]
        job(f"p4_2f3_{short(u)}", "cpu", 2, 10, f"{P}/f2v3/f2_{u}_floor.json", dep(u, "P4LOCK"), f"pixi run -e floors python scripts/p3_f2.py {u} floor", f"OMP_NUM_THREADS=4 SF_OUT={P}/f2v3 SF_CONFIG={c2}")
for u in F5U:  # InSituType (v2) + InSituType with the spatial prior (v3, PLAN §19) + logistic regression (secondary): all floor configs
    job(f"p4_5f_{short(u)}", "cpu", 4, 10, f"{P}/f5_{u}_floors.json", dep(u, "P4LOCK"), f"pixi run -e floors python scripts/p3_f5.py {u} floors")
for sid in (110, 120, 210, 221, 222, 242):  # L6 (retrospective slides whole, authors' protocol; 110 crop)
    big = sid in (110, 120, 210, 242)
    job(f"p4_4l_{sid}", "cpu", 12 if big else 6, 20 if big else 14, f"{P}/f4/l6/crc_{sid}.csv", "P4LOCK" + (",p4_4crop_110" if sid == 110 else ""),
        f"pixi run -e floors python scripts/f4_l6.py {sid}", R4 + (" SF_HOLDOUT=TVA" if sid == 222 else ""))
# Tier 2 (PROTOCOL §5, manifests/tier2.json): each model's grid configs beyond those tier 1 produces, behind gate T2 (opened
# only once the tier-1 report runs/p4/P5.md exists); outputs runs/p4/t2/<family>_<model>_<config> (F1: runs/p4/t2/f1)
if V3 and Path("manifests/tier2.json").exists():
    T2 = json.loads(Path("manifests/tier2.json").read_text()); SUB = set(T2["subset_crops"]); O2 = f"{P}/t2"; t2rows = len(rows)
    rows += ["# Tier 2 (manifests/tier2.json; PROTOCOL §5): after the tier-1 report", f"T2|cpu|.1|1|runs/sched/T2.ok|-|test -e runs/sched/T2.ok -a -s {P}/P5.md"]
    FU = {"F1": F1U, "F2": F2U, "F3": F3U, "F5": F5U}
    ENV = {("F1", "floor"): ("floors", "cpu", 4, 10), ("F1", "scvi"): ("scvi", "gpu", 2, 12), ("F2", "floor"): ("floors", "cpu", 1.5, 10), ("F2", "cellcharter"): ("cellcharter", "gpu", 4, 12),
           ("F2", "novae"): ("novae", "cpu", 7, 14), ("F3", "floor"): ("floors", "cpu", 2, 10), ("F3", "lightgbm"): ("scvi", "cpu", 2, 10), ("F3", "ncem"): (None, "cpu", 3, 10),
           ("F5", "floor"): ("floors", "cpu", 3, 10), ("F5", "scanvi"): ("scvi", "gpu", 4, 12), ("F5", "tangram"): ("scvi", "gpu", 5, 12)}
    RUN = {"F2": ("p3_f2.py", {"floor": "floor"}), "F3": ("p3_f3.py", {"floor": "ridge", "lightgbm": "lightgbm"}), "F5": ("p3_f5.py", {})}
    f1ids = []
    for f, ms in T2["grids"].items():
        for m, g in ms.items():
            for u in FU[f]:
                s = sp(u); lk = g["locked"].get(s)
                if lk is None or (g["units"] == "subset" and u not in SUB): continue
                extra = [c for c in g["configs"] if c != lk and c not in g["tier1"].get(s, [])]
                if not extra: continue
                env, cls, est, cap = ENV[(f, m)]; aft = "T2," + ("P4LOCK" if m == "floor" else "G7")
                if f == "F1":
                    for c in extra:
                        f1ids.append(f"t2_1{m[0]}_{short(u)}_{c}"); mm = "expr_embed" if m == "floor" else m
                        job(f1ids[-1], cls, est, cap, f"{O2}/f1/f1_{u}_{mm}_{c}.json", aft, f"pixi run -e {env} python scripts/f1_fit.py {u} {mm} {c}", f"SF_OUT={O2}/f1 OMP_NUM_THREADS=4")
                elif (f, m) == ("F3", "ncem"):
                    job(f"t2_3n_{short(u)}", cls, est, cap, f"runs/sched/t2_3n_{short(u)}.done", aft + ",p4_3ids", f"external/ncem_env/bin/python scripts/p3_f3_ncem_auth.py {u}",
                        f"{NT} SF_OUT={O2}/f3_ncem SF_CONFIG={','.join(extra)}")
                elif (f, m) == ("F5", "floor"):
                    job(f"t2_5f_{short(u)}", cls, est, cap, f"{O2}/f5_floor_prior/f5_{u}_prior.json", aft + f",p4_5f_{short(u)}", f"pixi run -e {env} python scripts/t2_f5_prior.py {u}",
                        f"OMP_NUM_THREADS=4 SF_OUT={O2}/f5_floor_prior")
                else:
                    script, fn = RUN[f]; fn = fn.get(m, m)
                    for c in extra:
                        d = f"{O2}/{f.lower()}_{m}_{c}"; xe = (" SF_F3_FEATS=v3" if f == "F3" else "") + (" SF_SCVI_SEEDS=0,1,2" if m == "cellcharter" else "")
                        job(f"t2_{f[1]}{m[0]}_{short(u)}_{c}", cls, est, cap, f"{d}/{f.lower()}_{u}_{fn}.json", aft + (",p4_3ids" if f == "F3" else ""),
                            f"pixi run -e {env} python scripts/{script} {u} {fn}", f"OMP_NUM_THREADS=4 SF_OUT={d} SF_CONFIG={c}{xe}")
    job("t2_1rescore", "cpu", 4, 14, f"{O2}/f1/f1_rescored.csv", ",".join(f1ids), "pixi run -e floors python scripts/p3_f1_rescore.py", f"SF_OUT={O2}/f1")
    print(len(rows) - t2rows - 2, "tier-2 jobs", end="; ")
Path("runs/sched/p4_jobs.txt").write_text("\n".join(rows) + "\n")
print(len(rows), "lines;", sum(r.split("|")[1] == "gpu" for r in rows if not r.startswith("#")), "gpu jobs ->", "runs/sched/p4_jobs.txt")
