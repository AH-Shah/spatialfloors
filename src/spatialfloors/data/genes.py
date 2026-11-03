"""Gene panel (SPEC §Gene panel): filter, then seurat_v3 HVG on X_a for whole-transcriptome panels; targeted panels keep
all filtered genes. Annotation: HGNC complete set 2026-10-06 (sha256 in HGNC_SHA), human; mouse panels use symbol rules."""
import re
from pathlib import Path
import numpy as np
import pandas as pd

HGNC = Path("/Users/Hemato/project_3/data/raw/hgnc_2026-10-06/hgnc_complete_set_2026-10-06.txt")
HGNC_SHA = "8bf6f686e640dbc864a346484626953621b0500d900d26cd085a35c2782d725d"
WTX_MIN = 10000  # panels with more genes than this are whole-transcriptome -> HVG 2000; else targeted (keep all filtered)
N_HVG = 2000

PROBE = re.compile(r"^(NEGATIVE|NegPrb|Neg_?Control|SYSTEMCONTROL|SystemControl|BLANK|Blank|FalseCode|falsecode|Unassigned|Deprecated|Intergenic)", re.I)
CLONE = re.compile(r"^((A[CFLP]|B[X]|C[RU]|FP|U|Z)\d{5,}\.\d+|ENS[A-Z]*G\d+(\.\d+)?|LOC\d+)$")  # unannotated loci without an HGNC symbol
UNANNOT_LNC = re.compile(r"^LINC\d+|-(AS|DT|IT|OT)\d*$")
PSEUDO = {"pseudogene", "immunoglobulin pseudogene", "T cell receptor pseudogene"}


def hgnc_table(path=HGNC):
    """symbol -> (locus_type, gene_group); previous / alias symbols resolve only when unique and not an approved symbol."""
    h = pd.read_csv(path, sep="\t", usecols=["symbol", "locus_type", "gene_group", "prev_symbol", "alias_symbol"], dtype=str).fillna("")
    main = {s: (lt, gg) for s, lt, gg in zip(h.symbol, h.locus_type, h.gene_group)}
    alt = {}
    for col in ("prev_symbol", "alias_symbol"):
        for sym, xs in zip(h.symbol, h[col]):
            for x in filter(None, xs.split("|")): alt.setdefault(x, set()).add(sym)
    for x, ss in alt.items():
        if x not in main and len(ss) == 1: main[x] = main[next(iter(ss))]
    return main


def filter_reason(symbols, species="human", table=None):
    """Per gene: '' (kept) or the removal reason. Human: HGNC locus type / gene group; mouse: symbol rules (targeted panels)."""
    if species == "human": table = table if table is not None else hgnc_table()
    out = []
    for g in map(str, symbols):
        if PROBE.match(g): out.append("probe"); continue
        if g.upper().startswith("MT-"): out.append("mito"); continue
        if species == "mouse":
            out.append("ribosomal" if re.match(r"^Rp[ls]\d+[a-z]?\d*$|^Rplp\d$|^Rpsa$", g) else ""); continue
        lt, gg = table.get(g, (None, ""))
        if lt is None: out.append("unannotated" if CLONE.match(g) or UNANNOT_LNC.search(g) else ""); continue
        groups = gg.split("|")
        if "L ribosomal proteins" in groups or "S ribosomal proteins" in groups: out.append("ribosomal")
        elif lt in PSEUDO: out.append("pseudogene")
        elif lt == "RNA, long non-coding" and UNANNOT_LNC.search(g): out.append("unannotated_lncRNA")
        else: out.append("")
    return np.array(out, dtype=object)


def panel(Xa, symbols, species="human", n_hvg=N_HVG):
    """Xa: cells x genes counts (X_a of the count split, or the counts when no split applies). Returns dict with the kept
    gene indices (filtered; HVG seurat_v3 on Xa for whole-transcriptome panels), the removal counts, and the mode."""
    import anndata as ad, scanpy as sc
    reason = filter_reason(symbols, species); keep = np.flatnonzero(reason == "")
    mode = "wtx_hvg" if len(symbols) > WTX_MIN else "targeted_all"
    if mode == "wtx_hvg":
        a = ad.AnnData(Xa[:, keep].astype(np.float32)); a.var_names = np.asarray(symbols)[keep].astype(str)
        sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=n_hvg, subset=False)
        keep = keep[np.flatnonzero(a.var.highly_variable.values)]
    return dict(mode=mode, genes=np.asarray(symbols)[keep].astype(str).tolist(), idx=keep.tolist(),
                removed={k: int(v) for k, v in pd.Series(reason[reason != ""]).value_counts().items()}, n_input=len(symbols))
