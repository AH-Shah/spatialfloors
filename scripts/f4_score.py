"""F4 scores of one GSE1 slide or MERFISH section (SPEC F4, co-primary): official pearson (Cellina eval_loo JSON for Cellina / CPA; floors.l6.score
for L6) and spike-free pearson (pseudocount .1, own top 50) for all three, the learned models' counterfactual pseudobulk
(per-cell CP10k of the saved cf, which is on the count scale, then the mean) scored against L6's control / target pseudobulks on the genes both share.
Usage (env floors): python scripts/f4_score.py <sid> -> runs/f4/f4_<sid>.csv"""
import glob, json, os, re, sys
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scanpy as sc
from spatialfloors.floors import l6
from spatialfloors.floors import _labels_to_coarse as lm
from spatialfloors.provenance import commit

sid = sys.argv[1]; MER = not sid.isdigit(); R = Path(os.environ.get("SF_F4_RUN", "runs/f4")); tag = sid if MER else f"crc_{sid}"
HDS = ["Isocortex", "Fiber_tracts"] if MER else [os.environ.get("SF_HOLDOUT", "CRC")]; L = pd.read_csv(R / f"l6/{tag}.csv"); Z = np.load(R / f"l6/{tag}_lfc.npz")


def l6_genes():  # the gene order of l6.load_crc / load_merfish + prep (same deterministic steps, names kept)
    if MER:
        A = ad.read_h5ad(R / f"raw/{sid}.h5ad"); a = ad.AnnData(X=A.raw.X, obs=pd.DataFrame(dict(ct=A.obs.cell_type.values, dom=A.obs.major_brain_region.values), index=A.obs_names), var=pd.DataFrame(index=A.raw.var_names))
        a.obs_names_make_unique(); a = a[(a.obs.ct.notna() & a.obs.dom.notna()).values].copy(); n = 1120
    else:
        f = R / f"raw/cosmx_wt_{sid}.h5ad"; a = ad.read_h5ad(f if f.exists() else l6.RAW / f"cosmx_wt_{sid}.h5ad"); a.obs_names_make_unique()
        ct = a.obs["ist"].astype(str).map(lm.LABEL_TO_COARSE); dom = a.obs["typ"].astype(str).str.extract(r"(REF|TVA|CRC)", expand=False)
        a._inplace_subset_obs((ct.notna() & dom.notna()).values); n = 2000
    sc.pp.filter_cells(a, min_counts=3); sc.pp.filter_genes(a, min_counts=3)
    sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=n, subset=True); return list(a.var_names)


G = l6_genes(); rows = [dict(sid=sid, model="L6", seed=0, holdout=r.holdout_celltype, pearson=r.pearson, pearson_clean=r.pearson_clean, n_shared=len(G))
                    for r in L[L.model_name == "L6_twoway_pseudobulk"].itertuples()]
base = R / (f"datasets/{sid}" if MER else f"datasets/crc/crc_{sid}"); cor = R / f"datasets/{'merfish' if MER else 'crc'}/correlations"
for hd in HDS:
    for f in sorted(glob.glob(str(base / f"*/*_s*_counterfactual_x_{hd}.h5ad"))):
        h = Path(f).parent.name; m, s = re.match(r"(\w+)_s(\d+)_counterfactual", Path(f).name).groups(); task = f"{h}_{hd}"
        a = ad.read_h5ad(f); gi = pd.Index(G).get_indexer(a.var_names); ok = gi >= 0
        X = np.asarray(a.X.toarray() if hasattr(a.X, "toarray") else a.X, np.float64); X = (X / (X.sum(1, keepdims=True) + 1e-8) * 1e4)[:, ok]  # cf is saved on the count scale: per-cell CP10k over its HVGs, as L6's nmean (audit 2026-10-30)
        mf = np.full(len(G), np.nan); mf[gi[ok]] = np.asarray(X.mean(0)).ravel(); keep = ~np.isnan(mf)
        sc_ = l6.score(Z[f"{task}|mc"][keep], Z[f"{task}|mt"][keep], mf[keep], Z[f"{task}|bleed"][keep])
        off = json.loads((cor / f"{tag}_{m}_s{s}-cf_{task}.json").read_text())["pearson"]
        rows.append(dict(sid=sid, model=m, seed=int(s), holdout=task, pearson=off, pearson_clean=sc_["pearson_clean"], pearson_l6genes=sc_["pearson"], n_shared=int(keep.sum())))
        assert keep.sum() >= .95 * len(G), f"{f}: cf genes do not match L6's HVGs ({keep.sum()} / {len(G)})"
d = pd.DataFrame(rows).assign(spatialfloors_commit=commit()); d.to_csv(R / f"f4_{sid}.csv", index=False)
if not os.environ.get("SF_BLIND"): print(d.groupby("model")[["pearson", "pearson_clean"]].mean().round(3))
