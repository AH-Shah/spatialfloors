"""Whole-slide arm job list (results/WHOLESLIDE_PLAN.md) for scripts/sched.sh -> runs/sched/xw_jobs.txt.
Largest slide first, so a memory problem shows on 120 (741k cells) before the others start. Locked configurations
(manifests/locked_v3.json), registered scripts, SF_LEAN=1, outputs in runs/xw.
Memory (GB): est at 741k cells from the gse1_231wl peaks (runs/xw/check: F1 floor 2.6, F2 floor 1.9, F3 ridge 4.6,
LightGBM 4.9, F5 floor 2.6 at 298k; F3 profiled per step) scaled to slide 120 (2.5x cells, 2.2x non-zeros per cell),
then x n_cells / 741k, at least 2.5. Every est <= 9.5, so each job fits beside the 2.5 GB reserve on an idle machine
(~12.4 GB free). cap (mem_guard kill limit) = est + 6, at most 16. Usage: python3 scripts/xw_jobs.py > runs/sched/xw_jobs.txt"""
N = {"120": 741155, "110": 694553, "210": 617498, "242": 474471, "221": 370074, "222": 140542}
EST_741K = {"prep": 6.5, "ids": 4.5, "1f": 7, "1s": 8.5, "2f": 4.5, "2c": 7, "3r": 9.5, "3l": 9.5, "3r2": 8, "5f": 8, "5s": 8}
E = "env SF_LEAN=1 SF_OUT=runs/xw OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4"
print("# whole-slide arm (exploratory): id|gpu or cpu|est_gb|cap_gb|done_file|after|command")
ids = []
for s, n in N.items():
    u = f"gse1_{s}w"
    J = [  # id, kind, done file, after, command
        ("prep", "cpu", f"manifests/panels/{u}.json", "-", f"{E} pixi run -e floors python scripts/xw_prepare.py {s}"),
        ("ids", "cpu", f"runs/p2/f3_ids/{u}.npy", "prep", f"{E} pixi run -e floors python scripts/f3_ids_dump.py {u}"),
        ("1f", "cpu", f"runs/xw/f1_{u}_expr_embed_pearson_none_pc50.json", "prep", f"{E} pixi run -e floors python scripts/f1_fit.py {u} expr_embed pearson_none_pc50"),
        ("1s", "gpu", f"runs/xw/f1_{u}_scvi_nl50_nb.json", "prep", f"{E} pixi run -e scvi python scripts/f1_fit.py {u} scvi nl50_nb"),
        ("2f", "cpu", f"runs/xw/f2_{u}_floor.json", "prep", f"{E} SF_CONFIG=ell400_lam0.8 pixi run -e floors python scripts/p3_f2.py {u} floor"),
        ("2c", "gpu", f"runs/xw/f2_{u}_cellcharter.json", "prep", f"{E} SF_CONFIG=nl10_L3 SF_SCVI_SEEDS=0,1,2 pixi run -e cellcharter python scripts/p3_f2.py {u} cellcharter"),
        ("3r", "cpu", f"runs/xw/f3_{u}_ridge.json", "ids", f"{E} SF_CONFIG=v3 SF_F3_FEATS=v3 pixi run -e floors python scripts/p3_f3.py {u} ridge"),
        ("3l", "cpu", f"runs/xw/f3_{u}_lightgbm.json", "ids", f"{E} SF_CONFIG=leaves15_lr0.05 SF_F3_FEATS=v3 pixi run -e scvi python scripts/p3_f3.py {u} lightgbm"),
        ("3r2", "cpu", f"runs/xw/v2/f3_{u}_ridge.json", "prep", f"{E} SF_OUT=runs/xw/v2 SF_CONFIG=T+ctx:ell25:pc10:a1000 pixi run -e floors python scripts/p3_f3.py {u} ridge"),
        ("5f", "cpu", f"runs/xw/f5_{u}_floors.json", "prep", f"{E} SF_CONFIG=insitutype_spatial pixi run -e floors python scripts/p3_f5.py {u} floors"),
        ("5s", "gpu", f"runs/xw/f5_{u}_scanvi.json", "prep", f"{E} SF_CONFIG=nl30_nb pixi run -e scvi python scripts/p3_f5.py {u} scanvi")]
    for i, k, d, after, cmd in J:
        est = max(2.5, round(EST_741K[i] * n / 741155, 1)); cap = min(16, round(est + 6)); jid = f"xw_{i}_{s}"; ids.append(jid)
        dep = f"xw_{after}_{s}" if after != "-" else ("-" if s == "120" else "xw_gate120")  # plan §4: other slides only after all of slide 120
        print(f"{jid}|{k}|{est}|{cap}|{d}|{dep}|{cmd}")
    if s == "120":  # gate: every slide-120 job done (a failure there stops the arm: scripts/xw_watch.sh)
        print(f"xw_gate120|cpu|0.5|2|runs/xw/GATE120.ok|{','.join(ids)}|touch runs/xw/GATE120.ok")
print(f"xw_analyze|cpu|4|10|runs/xw/X_wholeslide.json|{','.join(ids)}|pixi run -e floors python scripts/xw_analyze.py runs/xw runs/p4/P5_corrected.json")
