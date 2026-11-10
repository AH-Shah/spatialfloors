"""Nicheformer inputs per platform (authors' tokenization notebooks, external/nicheformer/notebooks/tokenization):
human CosMx -> symbols to Ensembl via the Pelka CELLxGENE var, technology mean cosmx; SPATCH Xenium -> the file's own
Ensembl ids, technology mean xenium, assay token "Xenium"; Allen MERFISH (Zhuang, Yao) -> mouse Ensembl ids mapped to
human orthologs ("Gene stable ID" -> "Human gene stable ID", first row per mouse gene as in merfish_mouse_brain.ipynb;
Ensembl release 116 BioMart export, sha256 below), technology mean merfish, species token mouse.
Review 2026-10-28: the P0/P3 adapter was human-CosMx-only; mouse genes did not map at all."""
from pathlib import Path
import numpy as np
import pandas as pd

from spatialfloors.paths import RAW
ORTHO = RAW / "ensembl_orthologs/mouse_human_orthologs.tsv"
ORTHO_SHA = "38d7c752d021d5088ffb2ac4d4dd263fdfda65cfb429522c25ee71c49d022d6f"


def _mouse_ids(name, symbols):
    import anndata as ad
    from ..data.units import ALLEN, YAO
    if name.startswith("allen_yao"):
        v = ad.read_h5ad(YAO / "C57BL6J-638850-raw.h5ad", backed="r").var; m = dict(zip(v.gene_symbol.astype(str), v.index.astype(str)))
    else:
        an = "animal2" if name.startswith("allen2_") else name.split("_")[1]; v = ad.read_h5ad(ALLEN[an], backed="r").raw.var
        m = dict(zip(v.feature_name.astype(str), v.index.astype(str)))
    return [m.get(str(s)) for s in symbols]


def inputs(name, symbols):
    """dict(ensg=list or None, specie, assay, tech) for run_nicheformer on crop `name` with var names `symbols`."""
    if name.startswith(("gse1_", "gse0_", "gse303070_", "bruker_")): return dict(ensg=None, specie="human", assay="cosmx", tech="cosmx")  # Bruker WTx = human CosMx
    if name.startswith("spatch"):
        import anndata as ad
        from ..data.units import RAW
        v = ad.read_h5ad(RAW / "spatch_xenium_5k/spatch_xenium_5k.h5ad", backed="r").var; m = dict(zip(v.index.astype(str), v.gene_ids.astype(str)))
        return dict(ensg=[m.get(str(s)) for s in symbols], specie="human", assay="Xenium", tech="xenium")
    o = pd.read_csv(ORTHO, sep="\t").drop_duplicates("Gene stable ID").set_index("Gene stable ID")["Human gene stable ID"]
    mid = _mouse_ids(name, symbols); ensg = [o.get(x) if x is not None and isinstance(o.get(x), str) else None for x in mid]
    return dict(ensg=ensg, specie="mouse", assay="merfish", tech="merfish")
