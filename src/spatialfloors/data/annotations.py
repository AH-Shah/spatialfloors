"""Pathologist region annotations (registration v2, ANNOTATION_PROTOCOL.md): QuPath GeoJSON polygons drawn on the H&E
image -> one class per Xenium cell; agreement between two annotators and their consensus.
Coordinates: he_imagealignment.csv (10x) is the 3 x 3 affine from H&E pixels to Xenium morphology-image pixels
(0.2125 um); a cell centroid in um is mapped to H&E pixels with its inverse. A cell takes the class of the smallest
polygon containing its centroid (an immune aggregate drawn inside stroma keeps its own class); "Exclude" overrides every
other class; cells in no polygon are unlabelled ("")."""

import json
import numpy as np
from matplotlib.path import Path as MPath

XENIUM_UM_PER_PX = 0.2125
CLASSES = ("Tumour", "Tumour stroma", "Normal tissue", "Immune aggregate", "Necrosis", "Other", "Exclude")
UNLABELLED = ("", "Exclude")
_ALIAS = {"tumor": "Tumour", "tumor stroma": "Tumour stroma", "stroma": "Tumour stroma", "normal": "Normal tissue"}


def canonical(name, extra=()):
    """Protocol class for a QuPath class name (case-insensitive; US spellings accepted); extra = tissue-specific
    classes fixed after the pilot."""
    allowed = {c.lower(): c for c in (*CLASSES, *extra)}
    n = str(name).strip().lower()
    if n in allowed:
        return allowed[n]
    if n in _ALIAS:
        return _ALIAS[n]
    raise ValueError(f"class {name!r} is not in the protocol ({', '.join((*CLASSES, *extra))})")


def read_geojson(path, extra=()):
    """[(class, [exterior ring, hole rings ...]), ...] from a QuPath 'Export objects as GeoJSON' file (pixel
    coordinates of the full-resolution H&E image). Polygon and MultiPolygon only; an unclassified polygon is an error."""
    d = json.loads(open(path).read())
    feats = d["features"] if isinstance(d, dict) and "features" in d else d if isinstance(d, list) else [d]
    out = []
    for i, f in enumerate(feats):
        g = f.get("geometry") or {}
        cl = (f.get("properties") or {}).get("classification")
        name = cl.get("name") if isinstance(cl, dict) else cl
        if name is None:
            raise ValueError(f"feature {i}: polygon without a class")
        if g.get("type") == "Polygon":
            polys = [g["coordinates"]]
        elif g.get("type") == "MultiPolygon":
            polys = g["coordinates"]
        else:
            raise ValueError(f"feature {i}: geometry {g.get('type')!r}, expected Polygon or MultiPolygon")
        c = canonical(name, extra)
        out += [(c, [np.asarray(r, float)[:, :2] for r in rings]) for rings in polys]
    return out


def um_to_he_px(xy_um, M):
    """Cell centroids (n x 2, um) -> H&E pixels, M = he_imagealignment.csv (H&E px -> Xenium px)."""
    p = np.linalg.inv(np.asarray(M, float)) @ np.vstack([np.asarray(xy_um, float).T / XENIUM_UM_PER_PX, np.ones(len(xy_um))])
    return (p[:2] / p[2]).T


def _area(r):
    x, y = r[:, 0], r[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def label_cells(xy_um, M, polys):
    """One class per cell (object array; "" = in no polygon)."""
    xy = um_to_he_px(xy_um, M)
    lab = np.full(len(xy), "", dtype=object)
    best = np.full(len(xy), np.inf)
    excl = np.zeros(len(xy), bool)
    for c, rings in polys:
        ext = rings[0]
        lo, hi = ext.min(0), ext.max(0)
        idx = np.flatnonzero(np.all((xy >= lo) & (xy <= hi), 1))
        if not len(idx):
            continue
        inside = MPath(ext).contains_points(xy[idx])
        for h in rings[1:]:
            inside &= ~MPath(h).contains_points(xy[idx])
        idx = idx[inside]
        if c == "Exclude":
            excl[idx] = True
            continue
        a = _area(ext) - sum(_area(h) for h in rings[1:])
        upd = idx[a < best[idx]]
        lab[upd], best[upd] = c, a
    lab[excl] = "Exclude"
    return lab


def agreement(a, b, min_share=0.02):
    """Two annotators' per-cell classes (same cell order) -> (summary dict, consensus labels). Scored cells: labelled by
    both (neither "" nor "Exclude"). Cohen's kappa, ARI and balanced ARI (the F2 metric, annotator a as the reference).
    Consensus: the class where both agree, else "". Regions: consensus classes with >= min_share of the consensus
    cells; the section enters F2 only with >= 2 regions."""
    from sklearn.metrics import adjusted_rand_score, cohen_kappa_score
    from spatialfloors.metrics.f2 import ari_bal

    a, b = np.asarray(a, dtype=object), np.asarray(b, dtype=object)
    both = ~np.isin(a, UNLABELLED) & ~np.isin(b, UNLABELLED)
    cons = np.where(both & (a == b), a, "").astype(object)
    lab = cons[cons != ""]
    share = {c: float((lab == c).mean()) for c in sorted(set(lab))} if len(lab) else {}
    regions = sorted(c for c, s in share.items() if s >= min_share)
    out = dict(
        n_cells=int(len(a)),
        n_labelled_both=int(both.sum()),
        n_agree=int(len(lab)),
        cohen_kappa=float(cohen_kappa_score(a[both], b[both])) if both.any() else None,
        ari=float(adjusted_rand_score(a[both], b[both])) if both.any() else None,
        ari_bal=float(ari_bal(a[both].astype(str), b[both].astype(str))) if both.any() else None,
        consensus_share=share,
        regions=regions,
        enters_f2=len(regions) >= 2,
    )
    return out, cons
