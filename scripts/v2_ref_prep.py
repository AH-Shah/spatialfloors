"""F5 references for the replication round (registration v2): one CELLxGENE dataset per new tissue, chosen by a fixed rule
from a saved CELLxGENE Discover dataset index, standardised as scripts/ref_prep.py does (layout read by
floors.f5.reference: X = log1p CP10k csr, raw.X = integer counts, var.feature_name, obs.ClusterMidway = fine label,
obs.coarse), with held-out markers by the ref_prep.py rule on the tissue's platform gene list.
Selection: Homo sapiens only; is_primary_data all true; suspension type cell only; every assay a 10x Chromium assay
(label starts with "10x", multiome excluded); not spatial; at least one tissue label containing the organ word (ORGAN);
every disease "normal" or in the tissue's carcinoma list (CARCINOMA). Largest cell_count; ties: earliest published_at.
Cells kept: tissue label containing the organ word and disease in the carcinoma list. Raw counts must be integers
(raw.X, else X); a dataset that fails is replaced by the next one in the ranking.
Classes (Cell Ontology, CL basic): epithelial = descendant of CL:0000066, or CL:0001064 / CL:0001061; immune non-myeloid =
descendant of CL:0000542; myeloid = descendant of CL:0000766 except mast cells CL:0000097; everything else other.
Held-out markers per class = top 5 by one-vs-rest AUROC (log1p CP10k, stratified <= 80k / <= 4k per fine type, seed 0)
among the 200 genes with the largest mean log-fold change, detected in >= 10% of the class's cells, on the gene list,
not in another class's top-20 list. QC: >= 500 cells per class.
Gene list (--genes): "wtx" = CosMx WTx list (GSE1 panel, the round-1 rule); "units:<u1,u2,...>" = genes measured on every
listed round-2 unit (Gene Expression features of its cache, else of its raw cell_feature_matrix.h5).
Usage (env floors):
  python scripts/v2_ref_prep.py select <index.json>              -> registration_v2/refs/selection.json
  python scripts/v2_ref_prep.py prep <tissue> <h5ad> --genes G [--dataset ID] [--out <json>]
                                                                 -> <h5ad dir>/ref.h5ad + registration_v2/refs/<tissue>.json
Check: prep breast <Wu 2021 h5ad> --genes wtx --out <tmp json> --ref-out <tmp h5ad> gives the held-out markers and
marker AUROCs of manifests/refs/breast.json."""

import argparse, hashlib, json
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from spatialfloors.paths import RAW

ROOT = Path(__file__).resolve().parents[1]
ORGAN = dict(lung="lung", ovary="ovar", prostate="prostat", breast="breast")
CARCINOMA = dict(
    lung=("lung adenocarcinoma", "squamous cell lung carcinoma", "non-small cell lung carcinoma", "lung large cell carcinoma"),
    ovary=(
        "malignant ovarian serous tumor",
        "ovarian serous carcinoma",
        "high grade ovarian serous adenocarcinoma",
        "ovarian carcinoma",
        "ovarian cancer",
    ),
    prostate=("prostate adenocarcinoma", "prostatic acinar adenocarcinoma", "prostate carcinoma", "prostate cancer"),
    breast=(
        "breast carcinoma",
        "breast cancer",
        "invasive ductal breast carcinoma",
        "invasive lobular breast carcinoma",
        "triple-negative breast carcinoma",
        "estrogen-receptor positive breast cancer",
        "HER2 positive breast carcinoma",
    ),
)
def eligible(x, t):
    lab = lambda k: [o["label"] for o in x[k]]
    return (
        lab("organism") == ["Homo sapiens"]
        and all(x["is_primary_data"])
        and x["suspension_type"] == ["cell"]
        and all(a.startswith("10x") and "multiome" not in a for a in lab("assay"))
        and not x.get("spatial")
        and any(ORGAN[t] in s.lower() for s in lab("tissue"))
        and all(s == "normal" or s in CARCINOMA[t] for s in lab("disease"))
        and any(s in CARCINOMA[t] for s in lab("disease"))
    )


def select(index):
    raw = Path(index).read_bytes()
    d = json.loads(raw)
    out = dict(
        index=str(index),
        index_sha256=hashlib.sha256(raw).hexdigest(),
        rule=__doc__.split("Classes")[0].strip(),
        tissues={},
    )
    for t in ORGAN:
        c = sorted((x for x in d if eligible(x, t)), key=lambda x: (-x["cell_count"], x["published_at"]))
        out["tissues"][t] = [
            dict(
                dataset_id=x["dataset_id"],
                title=x["title"],
                collection=x["collection_name"],
                citation=x["citation"].split(" Dataset Version")[0],
                cell_count=x["cell_count"],
                published_at=x["published_at"],
                diseases=[o["label"] for o in x["disease"]],
                tissues=[o["label"] for o in x["tissue"]],
                h5ad=next(a["url"] for a in x["assets"] if a["filetype"] == "H5AD"),
                h5ad_bytes=next(a["filesize"] for a in x["assets"] if a["filetype"] == "H5AD"),
            )
            for x in c
        ]
        print(t, [(x["dataset_id"][:8], x["cell_count"], x["title"][:50]) for x in out["tissues"][t]])
    p = ROOT / "registration_v2/refs/selection.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=1))
    print("->", p)


def gene_list(spec, root):
    import anndata as ad, h5py
    from spatialfloors.data.units_v2 import SECTIONS

    if spec == "wtx":
        return set(ad.read_h5ad(RAW / "GSE1/cosmx_wt_232.h5ad", backed="r").var_names.astype(str))
    sets = []
    for u in spec.removeprefix("units:").split(","):
        c = Path(root) / f"runs/p2/cache/{u}.h5ad"
        if c.exists():
            a = ad.read_h5ad(c, backed="r")
            sets.append(set(a.var_names[a.var.feature_type == "Gene Expression"].astype(str)))
            a.file.close()
        else:
            with h5py.File(SECTIONS[u][1] / "cell_feature_matrix.h5") as f:
                ft = f["matrix/features/feature_type"][:].astype(str)
                sets.append(set(f["matrix/features/name"][:].astype(str)[ft == "Gene Expression"]))
    return set.intersection(*sets)


def cl_ancestors(obo):
    par, cur = {}, None
    for line in open(obo):
        line = line.strip()
        if line.startswith("id: "):
            cur = line[4:]
        elif line.startswith("is_a: ") and cur:
            par.setdefault(cur, set()).add(line[6:].split(" ")[0])

    def anc(t, seen=None):
        seen = set() if seen is None else seen
        for p in par.get(t, ()):
            if p not in seen:
                seen.add(p)
                anc(p, seen)
        return seen

    return anc


def prep(t, h5ad, genes_spec, dataset, out_json, ref_out, root, chunk=50000):
    import anndata as ad, h5py
    from sklearn.metrics import roc_auc_score
    from spatialfloors.floors.f5 import stratified_rows
    from spatialfloors.metrics.f5 import CLASSES
    from spatialfloors.provenance import stamp

    f = h5py.File(h5ad, "r")
    obs = ad.io.read_elem(f["obs"])
    has_raw = "raw" in f
    var = ad.io.read_elem(f["raw/var"] if has_raw else f["var"])
    M = ad.io.sparse_dataset(f["raw/X"] if has_raw else f["X"])
    keep = obs.tissue.astype(str).str.lower().str.contains(ORGAN[t]).values & obs.disease.astype(str).isin(
        CARCINOMA[t]
    ).values
    rows = np.flatnonzero(keep)
    genes = var["feature_name"].astype(str).values
    anc = cl_ancestors(RAW / "refs/cl-basic.obo")

    def cls(term):
        a = anc(term) | {term}
        if "CL:0000097" in a:
            return "other"
        if "CL:0000066" in a or term in ("CL:0001064", "CL:0001061"):
            return "epithelial"
        if "CL:0000542" in a:
            return "immune_nonmyeloid"
        if "CL:0000766" in a:
            return "myeloid"
        return "other"

    fine = obs.cell_type.astype(str).values[rows]
    terms = obs.cell_type_ontology_term_id.astype(str).values[rows]
    cmap = {fl: cls(tm) for fl, tm in zip(fine, terms)}
    coarse = np.array([cmap[x] for x in fine])
    n_cls = {c: int((coarse == c).sum()) for c in (*CLASSES, "other")}
    r = stratified_rows(fine)  # positions within the kept cells
    ref_out = Path(ref_out) if ref_out else Path(h5ad).parent / "ref.h5ad"
    gw = h5py.File(ref_out, "w")
    Lr, first, o = [], True, 0  # o = kept cells before this block
    for s in range(0, M.shape[0], chunk):
        k = keep[s : s + chunk]
        if not k.any():
            continue
        X = sp.csr_matrix(M[s : s + len(k)][np.flatnonzero(k)], dtype=np.float32)
        if first:
            assert np.allclose(X.data[:100000], np.round(X.data[:100000])), "reference counts are not integers"
        lib = np.asarray(X.sum(1)).ravel()
        L = sp.csr_matrix(sp.diags(1e4 / np.maximum(lib, 1)) @ X)
        L.data = np.log1p(L.data)
        Lr.append(L[r[(r >= o) & (r < o + X.shape[0])] - o])
        o += X.shape[0]
        kw = dict(dataset_kwargs=dict(compression="gzip"))
        if first:
            ad.io.write_elem(gw, "X", L, **kw)
            ad.io.write_elem(gw, "raw/X", X, **kw)
            first = False
        else:
            ad.io.sparse_dataset(gw["X"]).append(L)
            ad.io.sparse_dataset(gw["raw/X"]).append(X)
    f.close()
    vdf = pd.DataFrame(dict(feature_name=genes), index=pd.Index(genes).astype(str))
    out = ad.AnnData(obs=pd.DataFrame(dict(ClusterMidway=fine, coarse=coarse)), var=vdf)
    out.var_names_make_unique()
    out.obs_names = out.obs_names.astype(str)
    out.strings_to_categoricals()  # as write_h5ad
    ad.io.write_elem(gw, "obs", out.obs)
    ad.io.write_elem(gw, "var", out.var)
    ad.io.write_elem(gw, "raw/var", out.var.copy())
    for k in ("obsm", "varm", "obsp", "varp", "layers", "uns", "raw/varm"):
        ad.io.write_elem(gw, k, {})
    gw.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})
    gw["raw"].attrs.update({"encoding-type": "raw", "encoding-version": "0.1.0"})
    gw.close()
    assert o == len(rows)
    # held-out markers
    gl = gene_list(genes_spec, root)
    Lr = sp.vstack(Lr).tocsc()
    cs = coarse[r]
    ok_g = np.array([g in gl for g in genes])
    mk = {}
    col = lambda j: Lr[:, j].toarray().ravel()
    for c in CLASSES:
        pos = cs == c
        P_, N_ = Lr[pos], Lr[~pos]
        det = np.asarray((P_ > 0).mean(0)).ravel()
        lfc = np.asarray(P_.mean(0)).ravel() - np.asarray(N_.mean(0)).ravel()
        lfc[~ok_g | (det < 0.1)] = -np.inf
        cand = np.argsort(-lfc)[:200]
        au = [(roc_auc_score(pos, col(j)), genes[j]) for j in cand if np.isfinite(lfc[j])]
        mk[c] = [g for _, g in sorted(au, reverse=True)]
    MK = {c: [g for g in mk[c] if g not in set().union(*(set(mk[o][:20]) for o in CLASSES if o != c))][:5] for c in CLASSES}
    gix = {g: i for i, g in enumerate(genes)}
    auc = {c: {g: float(roc_auc_score(cs == c, col(gix[g]))) for g in MK[c]} for c in CLASSES}
    rec = dict(
        tissue=t,
        source=f"CELLxGENE {dataset}" if dataset else str(h5ad),
        ref=str(ref_out),
        h5ad=str(h5ad),
        n_cells=int(len(fine)),
        class_counts=n_cls,
        qc_pass=bool(min(n_cls[c] for c in CLASSES) >= 500),
        coarse_map=cmap,
        gene_list=genes_spec,
        n_genes_in_list=int(ok_g.sum()),
        heldout_markers=MK,
        marker_auroc_in_reference=auc,
        rule=__doc__.split("Usage")[0].strip(),
    )
    p = Path(out_json) if out_json else ROOT / f"registration_v2/refs/{t}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(stamp(rec), indent=1, default=float))
    print(t, n_cls, "qc", rec["qc_pass"], MK, "->", p)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("select", "prep"))
    ap.add_argument("args", nargs="+")
    ap.add_argument("--genes")
    ap.add_argument("--dataset")
    ap.add_argument("--out")
    ap.add_argument("--ref-out")
    ap.add_argument("--root", default=".", help="round-2 runs root (runs/p2/cache) for --genes units:")
    A = ap.parse_args()
    if A.cmd == "select":
        select(A.args[0])
    else:
        prep(A.args[0], A.args[1], A.genes, A.dataset, A.out, A.ref_out, A.root)
