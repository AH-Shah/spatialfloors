"""F4 runner (P4 / validation): Cellina / CPA LOO on one GSE1 slide at the authors' configs (external/cellina-reproducibility@6bb060b),
seed s (train_loo / eval_loo DEFAULT_SEED; split and model seed alike). Generalises scripts/p0_cellina_loo.py (same device and
MPS workarounds, unchanged); outputs under runs/f4/datasets/crc/crc_<sid>/<holdout>/<model>_s<seed>_*.
Only change: the trainer device. scvi-tools 1.2 maps accelerator='auto' on Apple silicon to CPU, so Cellina gets
accelerator='mps' (CPA's scvi 0.20.3 maps use_gpu=True to MPS itself). Logged in the output JSON.
Usage (in env cellina or cpa): python scripts/f4_run.py {cellina,cpa} <sid> <seed> [holdout ...]"""
import os, sys, json, time, resource
from pathlib import Path
SF = Path(__file__).resolve().parents[1]; REPO = SF / "external/cellina-reproducibility"; RUN = SF / os.environ.get("SF_F4_RUN", "runs/f4")  # P4: runs/p4/f4
os.environ["DATA_ROOT"] = str(RUN); os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1"); os.chdir(REPO); sys.path[:0] = [str(REPO / "scripts")]
mc, SID, SEED = sys.argv[1], sys.argv[2], int(sys.argv[3]); holdouts = sys.argv[4:] or ["Endothelial", "Epithelial", "Fibroblast", "Myeloid", "T_cell"]
import torch
device_note = "authors' config"
NAME = f"{mc}_s{SEED}"
ACC = os.environ.get("SF_ACCEL", "mps"); MAXEP = os.environ.get("SF_MAX_EPOCHS"); TRAIN_ONLY = bool(os.environ.get("SF_TRAIN_ONLY"))
if mc == "cellina":
    import configs.cellina_config as C; C.TRAIN_ARGS["accelerator"] = ACC; C.TRAIN_ARGS["devices"] = 1; device_note = f"TRAIN_ARGS accelerator={ACC}, devices=1 (was devices=[0], accelerator auto -> CPU on Mac)"
    if MAXEP: C.TRAIN_ARGS["max_epochs"] = int(MAXEP); device_note += f"; DIAGNOSTIC max_epochs={MAXEP}"
if mc == "cpa" and ACC == "cpu":
    import configs.cpa_config as C2; C2.TRAIN_ARGS["use_gpu"] = False; device_note = "TRAIN_ARGS use_gpu=False (CPU; scvi 0.20.3 + torch 2.0 on MPS died with SIGBUS at epoch 1)"
    if MAXEP: C2.TRAIN_ARGS["max_epochs"] = int(MAXEP); device_note += f"; DIAGNOSTIC max_epochs={MAXEP}"
if os.environ.get("SF_THREADS"): torch.set_num_threads(int(os.environ["SF_THREADS"])); device_note += f"; torch threads {os.environ['SF_THREADS']}"
if mc == "cellina" and ACC == "cpu": device_note += " (MPS: vae_loss_validation NaN in epoch 1; CPU finite)"
if mc == "cellina" and ACC == "mps":
    # MPS workaround (P0 diagnosis, runs/p0_cellina/mps_diag2_cellina.json): torch 2.4.1 MPS lgamma returns non-finite values on
    # expanded (stride-0) views; scvi's NB broadcasts theta to mu's shape, so lgamma(theta) -> NaN loss. Contiguous copies are exact.
    import scvi.distributions._negative_binomial as nbm
    _lnb = nbm.log_nb_positive
    nbm.log_nb_positive = lambda x, mu, theta, eps=1e-8, log_fn=torch.log, lgamma_fn=torch.lgamma: _lnb(x, mu.contiguous(), theta.contiguous(), eps, log_fn, lgamma_fn)
    _lz = nbm.log_zinb_positive
    nbm.log_zinb_positive = lambda x, mu, theta, pi, eps=1e-8: _lz(x, mu.contiguous(), theta.contiguous(), pi.contiguous(), eps)
    device_note += "; MPS lgamma workaround: contiguous mu/theta in scvi NB/ZINB log-likelihood"
if mc == "cpa" and ACC == "mps":
    # MPS workaround (logging only): CPA logs numpy float64 epoch metrics; Lightning turns them into float64 MPS tensors (unsupported).
    import cpa._task as ct, numpy as _np
    _log = ct.CPATrainingPlan.log
    ct.CPATrainingPlan.log = lambda self, name, value, *a, **k: _log(self, name, float(value) if isinstance(value, (_np.floating, float)) else value, *a, **k)
    device_note = f"MPS (torch {torch.__version__}; torch 2.0.0 MPS nll_loss bus-errors on Darwin 25); float64 metric logging cast to float"
import train_loo, eval_loo
train_loo.DEFAULT_SEED = eval_loo.DEFAULT_SEED = SEED  # read at call time in main() / split_indices(seed=DEFAULT_SEED)
MER = not SID.isdigit(); DS = "merfish" if MER else "crc"  # MERFISH: SID = section label, file from f4_merfish_prep.py
if MER: holdouts = [h.replace("_", " ") for h in holdouts]  # job lines use GABAergic_neuron (no spaces)
HD = "Fiber_tracts" if MER else os.environ.get("SF_HOLDOUT", "CRC")  # 222: TVA (= adata_crc_holdout_config, which differs from adata_crc_config only here)
if not MER: train_loo.ADATA_CRC_ARGS["holdout_domains"] = [HD]; assert eval_loo.ADATA_CRC_ARGS["holdout_domains"] == [HD]  # same dict object
src = RUN / f"raw/{SID}.h5ad" if MER else next(f for f in (RUN / f"raw/cosmx_wt_{SID}.h5ad", Path(f"/Users/Hemato/project_3/data/raw/GSE1/cosmx_wt_{SID}.h5ad")) if f.exists())  # crop (f4_crop.py) first
lnk = RUN / (f"datasets/merfish/raw/{SID}.h5ad" if MER else f"datasets/crc/raw/crc_{SID}.h5ad"); lnk.parent.mkdir(parents=True, exist_ok=True)
if not lnk.exists(): lnk.symlink_to(src)
adata = str(lnk); log = RUN / f"timing_crc_{SID}_{NAME}{'_diag' if MAXEP else ''}.jsonl"
for h in holdouts:  # SF_PHASE=train|eval runs one phase per process (memory does not accumulate across tasks; P0 queue10/11 kills)
    done = RUN / (f"datasets/{SID}/{h}/" if MER else f"datasets/crc/crc_{SID}/{h}/") / f"{NAME}_counterfactual_x_{HD}.h5ad"; PH = os.environ.get("SF_PHASE", "both")
    t0 = time.time(); t_train = t_eval = None
    if PH in ("train", "both") and not (done.exists() and os.environ.get("SF_SKIP_DONE")):
        sys.argv = ["train_loo.py", "--dataset_name", DS, "--adata_path", adata, "--holdout_celltype", h, "--model_class", mc, "--model_name", NAME]
        train_loo.main(); t_train = round(time.time() - t0, 1)
    if PH in ("eval", "both") and not TRAIN_ONLY:
        assert done.exists(), f"no counterfactual for {h}: train first"; t1 = time.time()
        sys.argv = ["eval_loo.py", "--dataset_name", DS, "--adata_path", adata, "--holdout_celltype", h, "--model_class", mc, "--model_name", NAME, "--use_cf"]
        eval_loo.main(); t_eval = round(time.time() - t1, 1)
    rec = dict(sid=SID, seed=SEED, model=NAME, holdout=h, train_s=t_train, eval_s=t_eval, phase=PH, peak_rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30, 2),
               torch=torch.__version__, device=device_note, mps=torch.backends.mps.is_available())
    with open(log, "a") as f: f.write(json.dumps(rec) + "\n")
    print("TIMING", rec, flush=True)
    try: torch.mps.empty_cache()
    except Exception: pass
