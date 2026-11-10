"""F5 references for the new cancers, standardised to the Pelka layout read by floors.f5.reference:
X = log1p CP10k (csr), raw.X = integer counts, var.feature_name = gene symbol, obs.ClusterMidway = fine label,
obs.coarse.
Classes (Cell Ontology, CL basic): epithelial = descendant of epithelial cell CL:0000066, or malignant cell CL:0001064 /
abnormal cell CL:0001061 (carcinoma: tumour cells are epithelial); immune non-myeloid = descendant of lymphocyte
CL:0000542 (T, B, plasma, NK, ILC); myeloid = descendant of myeloid leukocyte CL:0000766 except mast cells CL:0000097
(as the locked Pelka map: Macro / Mono / DC / Granulo myeloid, Mast other); everything else other. Ji 2020 (GEO, no CL
terms): explicit level1 map below, level2 = fine label, Multiplet dropped.
Held-out markers per class = top 5 by one-vs-rest AUROC (log1p CP10k, stratified <= 80k / <= 4k per fine type, seed 0)
among the 200 genes with the largest mean log-fold change, detected in >= 10% of the class's cells, on the CosMx WTx
gene list (GSE1 panel), not in another class's list. QC: >= 500 cells per class (else the reference fails G5).
Usage (env floors): python scripts/ref_prep.py <breast|kidney|skin> -> <refdir>/ref.h5ad +
manifests/refs/<tissue>.json"""

import json, sys
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
from sklearn.metrics import roc_auc_score
from spatialfloors.floors.f5 import stratified_rows
from spatialfloors.metrics.f5 import CLASSES
from spatialfloors.provenance import stamp

T = sys.argv[1]
RD = __import__("spatialfloors.paths", fromlist=["RAW"]).RAW / "refs"
SRC = dict(
    breast=dict(
        dir="breast_wu2021",
        file="wu2021.h5ad",
        source="CELLxGENE 9fddb063-056d-4202-8b8a-4b0ee531d3ce (Wu et al. 2021, Nat Genet)",
    ),
    kidney=dict(
        dir="ccrcc_krishna2021",
        file="krishna2021.h5ad",
        source="CELLxGENE bd65a70f-b274-4133-b9dd-0d1431b6af34 (Krishna et al. 2021, Cancer Cell)",
    ),
    skin=dict(
        dir="cscc_ji2020",
        file="counts.txt.gz",
        source="GEO GSE144236 (Ji et al. 2020, Cell), level1 / level2 annotation",
    ),
)[T]
JI = {
    "Epithelial": "epithelial",
    "B Cell": "immune_nonmyeloid",
    "NK": "immune_nonmyeloid",
    "Tcell": "immune_nonmyeloid",
    "CD1C": "myeloid",
    "CLEC9A": "myeloid",
    "LC": "myeloid",
    "Mac": "myeloid",
    "MDSC": "myeloid",
    "ASDC": "myeloid",
    "PDC": "myeloid",
    "Endothelial Cell": "other",
    "Fibroblast": "other",
    "Melanocyte": "other",
}
# Krishna 2021: CELLxGENE codes the TAM clusters as 'abnormal
# cell' (80% of the CL 'epithelial' class) and PAX8+ renal epithelium as myofibroblast, so classes come from
# author_cluster_name.
# Ambiguous clusters are dropped (as Ji 2020 Multiplet).
KR_DROP = ("Ambiguous", "Ambiguous/Dead", "TAM/TCR (Ambiguos)")


def kr(n):
    if n in ("CD45- ccRCC CA9+", "CD45- PAX8+ renal epithelium"):
        return "epithelial"
    if n.startswith(("CD4+", "CD8A+", "NK ")) or n in ("B cell", "Conventional NK"):
        return "immune_nonmyeloid"
    if n.startswith(("TAM", "CD14+")) or n in ("cDC1", "cDC2", "pDC"):
        return "myeloid"
    return "other"  # Myofibroblast, Vascular Endothelium, Mast, Megakaryocyte


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


d = RD / SRC["dir"]
if T == "skin":
    m = pd.read_csv(d / "meta.txt.gz", sep="\t")
    m = m[m.level1_celltype != "Multiplet"]
    blk, gl, cols = [], [], None  # genes x cells text, read in gene blocks into sparse (dense whole file > 14 GB)
    for ch in pd.read_csv(d / "counts.txt.gz", sep="\t", index_col=0, chunksize=2000):
        cols = ch.columns if cols is None else cols
        blk.append(sp.csr_matrix(ch.values.astype(np.float32)))
        gl += list(ch.index.astype(str))
    keep_c = np.flatnonzero(np.isin(cols, m.index))
    X = sp.vstack(blk).T.tocsr()[keep_c]
    m = m.loc[cols[keep_c]]
    genes = np.array(gl)
    fine = m.level2_celltype.astype(str).values
    coarse = m.level1_celltype.map(JI).values
    cmap = {f: JI[l1] for f, l1 in zip(fine, m.level1_celltype)}
else:
    A = ad.read_h5ad(d / SRC["file"])
    X = sp.csr_matrix(A.raw.X if A.raw is not None else A.X, dtype=np.float32)
    genes = (A.raw.var if A.raw is not None else A.var)["feature_name"].astype(str).values
    obo = RD / "cl-basic.obo"
    anc = cl_ancestors(obo)
    terms = A.obs.cell_type_ontology_term_id.astype(str).values

    def cls(t):
        a = anc(t) | {t}
        if "CL:0000097" in a:
            return "other"
        if "CL:0000066" in a or t in ("CL:0001064", "CL:0001061"):
            return "epithelial"
        if "CL:0000542" in a:
            return "immune_nonmyeloid"
        if "CL:0000766" in a:
            return "myeloid"
        return "other"

    fine = A.obs.cell_type.astype(str).values
    cmap = {f: cls(t) for f, t in zip(fine, terms)}
    coarse = np.array([cmap[f] for f in fine])
    if T == "kidney":
        an = A.obs.author_cluster_name.astype(str).values
        keep_c = ~np.isin(an, KR_DROP)
        X, fine = X[keep_c], an[keep_c]
        cmap = {f: kr(f) for f in set(fine)}
        coarse = np.array([cmap[f] for f in fine])
assert np.allclose(X.data[:100000], np.round(X.data[:100000])), "reference counts are not integers"
n_cls = {c: int((coarse == c).sum()) for c in (*CLASSES, "other")}
lib = np.asarray(X.sum(1)).ravel()
L = sp.csr_matrix(sp.diags(1e4 / np.maximum(lib, 1)) @ X)
L.data = np.log1p(L.data)
# held-out markers
wtx = set(ad.read_h5ad(RD.parent / "GSE1/cosmx_wt_232.h5ad", backed="r").var_names.astype(str))
r = stratified_rows(fine)
Lr = L[r].tocsc()
cs = coarse[r]
ok_g = np.array([g in wtx for g in genes])
mk = {}  # sparse until the 200 candidates (memory)
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
MK = {}
for c in CLASSES:
    MK[c] = [g for g in mk[c] if g not in set().union(*(set(mk[o][:20]) for o in CLASSES if o != c))][:5]
gix = {g: i for i, g in enumerate(genes)}
auc = {c: {g: float(roc_auc_score(cs == c, col(gix[g]))) for g in MK[c]} for c in CLASSES}
out = ad.AnnData(
    L,
    obs=pd.DataFrame(dict(ClusterMidway=fine, coarse=coarse)),
    var=pd.DataFrame(dict(feature_name=genes), index=pd.Index(genes).astype(str)),
)
out.var_names_make_unique()
out.raw = ad.AnnData(X, var=out.var.copy())
out.write_h5ad(d / "ref.h5ad", compression="gzip")
rec = dict(
    tissue=T,
    source=SRC["source"],
    ref=str(d / "ref.h5ad"),
    n_cells=int(len(fine)),
    class_counts=n_cls,
    qc_pass=bool(min(n_cls[c] for c in CLASSES) >= 500),
    coarse_map=cmap,
    heldout_markers=MK,
    marker_auroc_in_reference=auc,
    rule=__doc__.split("Usage")[0].strip(),
)
Path("manifests/refs").mkdir(parents=True, exist_ok=True)
Path(f"manifests/refs/{T}.json").write_text(json.dumps(stamp(rec), indent=1, default=float))
print(T, n_cls, "qc", rec["qc_pass"], MK)
