"""F5 label transfer on one crop (primary = held-out-marker AUROC, metrics.f5.heldout_auroc; secondary = IF-gated
macro-F1, GSE1 only); every config of the model's grid, or only SF_CONFIG. Inputs = panel - F5 held-out markers,
intersected with the reference (Pelka for CRC, the unit's own reference for new cancers, floors.f5.for_unit).
  floors   InSituType ML, InSituType + spatial prior (floors.f5.spatial_prior, the locked floor), logistic regression
           C {.1, 1} x scaling {own, ref} (secondary)                                                        (env floors, R)
  scanvi   n_latent {10, 30} x likelihood {nb, zinb}; reference raw counts (stratified)                    (env scvi, MPS)
  tangram  density prior {uniform, rna_count_based}                                                         (env scvi, MPS)
Usage: python scripts/p3_f5.py <unit> <model> -> $SF_OUT/f5_<unit>_<model>.json (default runs/p3)"""
import json, sys, time
from pathlib import Path
import anndata as ad, numpy as np, scipy.sparse as sp
from spatialfloors.floors import f5
from spatialfloors.metrics.f5 import unit_inputs, score_classes
from spatialfloors.provenance import stamp

import os
u, model = sys.argv[1], sys.argv[2]; OUT = Path(os.environ.get("SF_OUT", "runs/p3")); OUT.mkdir(parents=True, exist_ok=True); ONLY = os.environ.get("SF_CONFIG")  # P4: one locked config
want = lambda c: ONLY in (None, c)
MK = f5.for_unit(u, json.loads(Path("manifests/panels/index.json").read_text())["f5_heldout_markers"]); held = {x for v in MK.values() for x in v}  # new cancers: own reference (PLAN §17)
a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad"); p = json.loads(Path(f"manifests/panels/{u}.json").read_text())
pg = set(f5.pelka_genes()); genes = [x for x in p["genes"] if x not in held and x in pg]; assert not held & set(genes); Xq = sp.csr_matrix(a[:, genes].X)
Xf, vf, gate, cls = unit_inputs(u, a)
LAB = {}


def score(fine, name):
    LAB[name] = np.asarray(fine).astype(str); return score_classes(f5.coarse(fine), Xf, vf, MK, gate, cls)


rep = dict(unit=u, model=model, n_cells=int(a.n_obs), n_genes=len(genes), grid={})
if model == "floors":
    negc = [x for x in a.var_names if x.startswith(("NegPrb", "Negative"))]; neg = np.asarray(a[:, negc].X.mean(1)).ravel() if negc else None  # Bruker WTx names them Negative<n> (SPEC F5: use negative probes when present; audit 2026-10-30)
    if want("insitutype") or want("insitutype_spatial"):  # v3 (PLAN §19): + spatial prior on the same run's log-likelihoods
        t = time.time(); prof, names = f5.reference(genes, rows="all"); lab, _ = f5.insitutype(Xq, prof, list(names), neg, workdir=OUT / f"ist_{u}")
        rep["grid"]["insitutype"] = dict(score(lab, "insitutype"), s=round(time.time() - t, 1))
        import pandas as pd
        LL = pd.read_csv(OUT / f"ist_{u}" / "logliks.csv"); lab_s = f5.spatial_prior(LL.values, list(LL.columns), a.obs[["x_um", "y_um"]].values.astype(float))
        rep["grid"]["insitutype_spatial"] = dict(score(lab_s, "insitutype_spatial"), s=round(time.time() - t, 1))
    Xr, yr = f5.reference(genes)
    for C in (.1, 1.):
        for sc in ("own", "ref"):
            if not want(f"logreg_C{C}_{sc}"): continue
            t = time.time(); rep["grid"][f"logreg_C{C}_{sc}"] = dict(score(f5.logreg(Xr, yr, Xq, C, sc), f"logreg_C{C}_{sc}"), s=round(time.time() - t, 1))
else:
    from spatialfloors.models import transfer
    import torch
    Xr, yr = f5.reference(genes, counts=True)
    cfgs = {f"nl{nl}_{lk}": dict(n_latent=nl, likelihood=lk) for nl in (10, 30) for lk in ("nb", "zinb")} if model == "scanvi" else {d: dict(density=d) for d in ("uniform", "rna_count_based")}
    for c, kw in cfgs.items():
        if not want(c): continue
        t = time.time(); lab = transfer.scanvi(Xr, yr, Xq, **kw) if model == "scanvi" else transfer.tangram(Xr, yr, Xq, **kw)
        rep["grid"][c] = dict(score(lab, c), s=round(time.time() - t, 1)); torch.mps.empty_cache()
for c, v in rep["grid"].items(): print(u, model, c, round(v["heldout_auroc"], 3), v.get("gated_macro_f1"), v["s"], flush=True)
(OUT / "labels").mkdir(exist_ok=True); np.savez_compressed(OUT / "labels" / f"f5_{u}_{model}.npz", **LAB)  # fine predictions (stress panel)
(OUT / f"f5_{u}_{model}.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
