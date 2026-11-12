"""Exploratory scale curves: job list for scripts/sched.sh (format: id|gpu or cpu|est_gb|cap_gb|done_file|after|command).
Per GSE1 confirmation slide: one window-preparation job (scripts/x_scale_prepare.py; the slides run one after another
because they append to the manifest indexes), then per window the whole-slide arm's jobs (scripts/xw_jobs.py): F1 floor
and scVI, F2 floor and CellCharter, F3 receiver ids, floor and LightGBM, F5 floor and scANVI; locked configurations,
registered scripts, SF_LEAN=1, outputs in runs/xs. Smallest windows first. Memory estimates as xw_jobs.py, scaled by
cells. The code is this repository's (scripts/, PYTHONPATH src/); the working directory and the pixi environments are
the round-1 runs root's.
Usage: python3 scripts/x_scale_jobs.py <pixi executable> > <round-1 runs root>/runs/sched/xs_jobs.txt"""

import sys
from pathlib import Path

PUB = Path(__file__).resolve().parents[1]
PIXI = sys.argv[1]
N = {"120": 741155, "110": 694553, "210": 617498, "242": 474471, "221": 370074, "222": 140542}
SIZES = (12500, 25000, 50000, 100000, 200000, 400000)
EST_741K = {"ids": 4.5, "1f": 7, "1s": 8.5, "2f": 4.5, "2c": 7, "3r": 9.5, "3l": 9.5, "5f": 8, "5s": 8}
E = (
    f"env PATH={Path(PIXI).parent}:/usr/bin:/bin:/usr/sbin:/sbin PIXI_FROZEN=true SF_LEAN=1 SF_OUT=runs/xs SF_CROPS=manifests/crops PYTHONPATH={PUB}/src "
    "OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 "
    "OPENBLAS_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4"
)
S = f"{PUB}/scripts"
print("# exploratory scale curves: id|gpu or cpu|est_gb|cap_gb|done_file|after|command")
prev, jobs = "-", []
for s in N:
    last = max(n for n in SIZES if N[s] > 1.1 * n)  # the slide's largest window: its panel file is written last
    print(f"xs_prep_{s}|cpu|8|23|manifests/panels/gse1_{s}s{last}.json|{prev}|{E} {PIXI} run -e floors python {S}/x_scale_prepare.py {s}")
    prev = f"xs_prep_{s}"
for n in SIZES:
    for s, ns in N.items():
        if ns <= 1.1 * n:
            continue
        u = f"gse1_{s}s{n}"
        J = [
            ("ids", "cpu", f"runs/p2/f3_ids/{u}.npy", f"xs_prep_{s}", f"{PIXI} run -e floors python {S}/f3_ids_dump.py {u}"),
            ("1f", "cpu", f"runs/xs/f1_{u}_expr_embed_pearson_none_pc50.json", f"xs_prep_{s}",
             f"{PIXI} run -e floors python {S}/f1_fit.py {u} expr_embed pearson_none_pc50"),
            ("1s", "gpu", f"runs/xs/f1_{u}_scvi_nl50_nb.json", f"xs_prep_{s}", f"{PIXI} run -e scvi python {S}/f1_fit.py {u} scvi nl50_nb"),
            ("2f", "cpu", f"runs/xs/f2_{u}_floor.json", f"xs_prep_{s}",
             f"env SF_CONFIG=ell400_lam0.8 {PIXI} run -e floors python {S}/p3_f2.py {u} floor"),
            ("2c", "gpu", f"runs/xs/f2_{u}_cellcharter.json", f"xs_prep_{s}",
             f"env SF_CONFIG=nl10_L3 SF_SCVI_SEEDS=0,1,2 {PIXI} run -e cellcharter python {S}/p3_f2.py {u} cellcharter"),
            ("3r", "cpu", f"runs/xs/f3_{u}_ridge.json", f"xs_ids_{u}",
             f"env SF_CONFIG=v3 SF_F3_FEATS=v3 {PIXI} run -e floors python {S}/p3_f3.py {u} ridge"),
            ("3l", "cpu", f"runs/xs/f3_{u}_lightgbm.json", f"xs_ids_{u}",
             f"env SF_CONFIG=leaves15_lr0.05 SF_F3_FEATS=v3 {PIXI} run -e scvi python {S}/p3_f3.py {u} lightgbm"),
            ("5f", "cpu", f"runs/xs/f5_{u}_floors.json", f"xs_prep_{s}",
             f"env SF_CONFIG=insitutype_spatial {PIXI} run -e floors python {S}/p3_f5.py {u} floors"),
            ("5s", "gpu", f"runs/xs/f5_{u}_scanvi.json", f"xs_prep_{s}",
             f"env SF_CONFIG=nl30_nb {PIXI} run -e scvi python {S}/p3_f5.py {u} scanvi"),
        ]
        for i, k, d, after, cmd in J:
            print(f"xs_{i}_{u}|{k}|{max(2.5, round(EST_741K[i] * n / 741155, 1))}|23|{d}|{after}|{E} {cmd}")
