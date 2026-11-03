#!/bin/bash
# NCEM reproduction on the authors' MIBI-TOF data (Hartmann et al.; Zenodo 3951613), user 2026-11-02: run the authors'
# own grid-search scripts (theislab/ncem_benchmarks scripts/grid_searches/mibitof_cancer/{ed_ncem_gcn,ed,interactions,
# interaction_baseline}.sh) with their exact argument keys, sequentially on CPU instead of SLURM. Published targets
# (model_benchmarks_mibitof_cancer.ipynb, test r_squared_linreg, best on val): ED ~.69, ED-NCEM GCN ~.70-.715 (13-46 um),
# linear interactions ~.80. Usage: bash scripts/ncem_repro_mibi.sh {profile|ed_ncem_gcn|ed|interactions|interaction_baseline}
cd /Users/Hemato/project_3/antig/spatialfloors
PY=external/ncem_env/bin/python; S=external/ncem_benchmarks_src/ncem_benchmarks/scripts; DATA=runs/p3/ncem_repro/data; OUT=runs/p3/ncem_repro/gs
export OMP_NUM_THREADS=4 TF_NUM_INTRAOP_THREADS=4 TF_NUM_INTEROP_THREADS=2 TF_CPP_MIN_LOG_LEVEL=2
run_ed_ncem() {  # args: GS_KEY radius enc_depth  (ed_ncem_gcn.sh: LR 1, LAT 1, DR 1, L1 1, L2 2, EID 2, DID 1, DDE 1, CDE 1, CDI 0, CDR 1, CL2 2, BS S, NEVAL 10, TYPE_COND 1)
  mkdir -p $OUT/$1/results $OUT/$1/logs
  $PY $S/train_script_ed_ncem.py HARTMANN ADAM GCN IMAGE 1 1 1 1 2 2 $3 1 1 1 0 1 2 S $2 10 1 ED_NCEM $1 $DATA $OUT/$1
}
case $1 in
  profile) date +%T; run_ed_ncem profile_ED_NCEM_GCN 1 0; date +%T ;;
  ed_ncem_gcn) for ede in 0 1; do for rd in 0 1 2 3 4 5; do run_ed_ncem ED_NCEM_GCN_IMAGE_HARTMANN $rd $ede || echo FAIL ed_ncem $rd $ede; done; done ;;
  ed) mkdir -p $OUT/ED_IMAGE_HARTMANN/results  # ed.sh: LR 1, LAT 2+3, DR 1, L1 {1,4,6}, L2 1+4+5, EID 4, EDE 0, DID 4, DDE 0, BS S, RD 0, NEVAL 10+100, TYPE_COND 1
      for l1 in 1 4 6; do $PY $S/train_script_ed.py HARTMANN ADAM IMAGE 1 2+3 1 $l1 1+4+5 4 0 4 0 S 0 10+100 1 ED ED_IMAGE_HARTMANN $DATA $OUT/ED_IMAGE_HARTMANN || echo FAIL ed $l1; done ;;
  interactions) mkdir -p $OUT/INTERACTIONS_IMAGE_HARTMANN/results
      for rd in 0 1 2 3 4 5; do $PY $S/train_script_linear.py HARTMANN ADAM IMAGE 1 1 1+4+6 S $rd 10 INTERACTIONS INTERACTIONS_IMAGE_HARTMANN $DATA $OUT/INTERACTIONS_IMAGE_HARTMANN || echo FAIL inter $rd; done ;;
  interaction_baseline) mkdir -p $OUT/INTERACTIONS_BASELINE_IMAGE_HARTMANN/results
      $PY $S/train_script_linear.py HARTMANN ADAM IMAGE 1 1 1 S 0 10 INTERACTIONS_BASELINE INTERACTIONS_BASELINE_IMAGE_HARTMANN $DATA $OUT/INTERACTIONS_BASELINE_IMAGE_HARTMANN || echo FAIL ib ;;
esac
echo REPRO_DONE $1 $(date +%H:%M)
