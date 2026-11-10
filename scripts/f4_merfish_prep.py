"""F4 MERFISH input: one h5ad per F4 section (manifests/allen.json f4_sections) in Cellina's Allen schema:
raw.X = counts, obs cell_type / major_brain_region / brain_section_label, obsm X_spatial_coords (um), as the authors'
animal 2
file. Zhuang animals: rows of the CELLxGENE file. Yao 638850: types and regions mapped as for the section choice: class
Glut / GABA -> glutamatergic / GABAergic neuron, subclass Astro / Endo / Oligo NN -> astrocyte / endothelial cell /
oligodendrocyte, other cells keep their class name; TH / Isocortex / 'fiber tracts' -> Thalamus / Isocortex /
Fiber_tracts.
Usage (env floors): SF_F4_RUN=runs/p4/f4 python scripts/f4_merfish_prep.py <animal1|animal3|yao638850>"""

import json, os, sys
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
from spatialfloors.data.units import ALLEN, YAO

an = sys.argv[1]
secs = json.loads(Path("manifests/allen.json").read_text())[an]["f4_sections"]
out = Path(os.environ.get("SF_F4_RUN", "runs/f4")) / "raw"
out.mkdir(parents=True, exist_ok=True)
if an.startswith("animal"):
    A = ad.read_h5ad(ALLEN[an], backed="r")
    for s in secs:
        m = (A.obs.brain_section_label == s).values
        rows = np.flatnonzero(m)
        b = ad.AnnData(
            sp.csr_matrix(A.raw.X[rows]),
            obs=A.obs.iloc[rows][["cell_type", "major_brain_region", "brain_section_label"]].copy(),
            var=A.raw.var.copy(),
        )
        b.obsm["X_spatial_coords"] = np.asarray(A.obsm["X_spatial_coords"][rows])
        b.raw = b.copy()
        b.write_h5ad(out / f"{s}.h5ad")
        print(s, b.shape, flush=True)
else:
    T = ["GABAergic neuron", "astrocyte", "endothelial cell", "glutamatergic neuron", "oligodendrocyte"]
    m = pd.read_csv(
        YAO / "cell_metadata_with_parcellation_annotation.csv",
        usecols=[
            "cell_label",
            "brain_section_label",
            "class",
            "subclass",
            "parcellation_division",
            "parcellation_category",
            "x_section",
            "y_section",
        ],
    )
    m = m[m.brain_section_label.isin(secs)].set_index("cell_label")
    c, s_ = m["class"].astype(str), m["subclass"].astype(str)
    ct = np.select(
        [
            c.str.contains("Glut"),
            c.str.contains("GABA"),
            s_.str.contains(r"\bAstro"),
            s_.str.contains(r"\bEndo\b"),
            s_.str.contains(r"\bOligo NN"),
        ],
        T[3:4] + T[0:1] + T[1:3] + T[4:],
        default=c.values,
    )
    reg = np.select(
        [
            m.parcellation_division == "TH",
            m.parcellation_division == "Isocortex",
            m.parcellation_category == "fiber tracts",
        ],
        ["Thalamus", "Isocortex", "Fiber_tracts"],
        default=m.parcellation_division.astype(str),
    )
    m = m.assign(cell_type=ct, major_brain_region=reg)
    A = ad.read_h5ad(YAO / "C57BL6J-638850-raw.h5ad", backed="r")
    pos = pd.Index(A.obs_names.astype(str)).get_indexer(m.index.astype(str))
    assert (pos >= 0).all()
    o = np.argsort(pos)
    X = sp.csr_matrix(A.X[pos[o]])
    inv = np.empty_like(o)
    inv[o] = np.arange(len(o))
    X = X[inv]
    var = A.var.copy()
    var.index = (var["gene_symbol"] if "gene_symbol" in var else var.index).astype(str).values
    A.file.close()
    for s in secs:
        k = (m.brain_section_label == s).values
        b = ad.AnnData(
            X[np.flatnonzero(k)],
            obs=m.loc[k, ["cell_type", "major_brain_region", "brain_section_label"]].copy(),
            var=var.copy(),
        )
        b.var_names_make_unique()
        b.obsm["X_spatial_coords"] = m.loc[k, ["x_section", "y_section"]].values * 1000
        b.raw = b.copy()
        b.write_h5ad(out / f"{s}.h5ad")
        print(s, b.shape, flush=True)
