"""Pathologist annotations -> per-cell classes (registration v2, ANNOTATION_PROTOCOL.md; data.annotations).
  cells: one annotator's QuPath GeoJSON (H&E pixels) on one 10x section -> <out>.tsv.gz (cell_id, class) + <out>.json
         (class counts, files and sha256)
  agree: two annotators' cell tables of one section -> <out>.tsv.gz (consensus class) + <out>.json (Cohen's kappa,
         ARI, balanced ARI, consensus shares, regions >= 2%, whether the section enters F2)
Usage (env floors):
  python scripts/v2_annot_to_cells.py cells <section> <geojson> <out stem> [--extra "Class,Class"]
  python scripts/v2_annot_to_cells.py agree <cells A .tsv.gz> <cells B .tsv.gz> <out stem>"""

import argparse, hashlib, json
from pathlib import Path
import numpy as np, pandas as pd
from spatialfloors.data import annotations as an
from spatialfloors.data.units_v2 import SECTIONS
from spatialfloors.provenance import stamp

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
ap = argparse.ArgumentParser()
ap.add_argument("cmd", choices=("cells", "agree"))
ap.add_argument("args", nargs=3)
ap.add_argument("--extra", default="")
A = ap.parse_args()
out = Path(A.args[2])
out.parent.mkdir(parents=True, exist_ok=True)
if A.cmd == "cells":
    sec, gj = A.args[0], A.args[1]
    d = SECTIONS[sec][1]
    extra = tuple(x.strip() for x in A.extra.split(",") if x.strip())
    polys = an.read_geojson(gj, extra)
    M = np.loadtxt(d / "he_imagealignment.csv", delimiter=",")
    c = pd.read_parquet(d / "cells.parquet", columns=["cell_id", "x_centroid", "y_centroid"])
    lab = an.label_cells(c[["x_centroid", "y_centroid"]].values, M, polys)
    pd.DataFrame({"cell_id": c.cell_id.astype(str).values, "class": lab}).to_csv(
        out.with_suffix(".tsv.gz"), sep="\t", index=False
    )
    rec = dict(
        section=sec,
        geojson=str(gj),
        geojson_sha256=sha(gj),
        alignment_sha256=sha(d / "he_imagealignment.csv"),
        n_polygons=len(polys),
        polygons_per_class=pd.Series([p[0] for p in polys]).value_counts().to_dict(),
        n_cells=int(len(lab)),
        cells_per_class=pd.Series(lab).replace("", "(none)").value_counts().to_dict(),
        extra_classes=list(extra),
    )
else:
    a, b = (pd.read_csv(p, sep="\t", keep_default_na=False) for p in A.args[:2])
    assert (a.cell_id.values == b.cell_id.values).all(), "cell tables differ in cells or order"
    rec, cons = an.agreement(a["class"].values, b["class"].values)
    pd.DataFrame({"cell_id": a.cell_id.values, "class": cons}).to_csv(out.with_suffix(".tsv.gz"), sep="\t", index=False)
    rec.update(cells_a=A.args[0], cells_a_sha256=sha(A.args[0]), cells_b=A.args[1], cells_b_sha256=sha(A.args[1]))
out.with_suffix(".json").write_text(json.dumps(stamp(rec), indent=1, default=float))
print(json.dumps({k: v for k, v in rec.items() if not k.endswith("sha256")}, default=float))
