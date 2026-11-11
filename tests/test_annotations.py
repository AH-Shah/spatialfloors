import json
import numpy as np, pytest
from spatialfloors.data import annotations as an

# H&E px -> Xenium px: x2 scale and a 90-degree rotation, as in the 10x files (H&E portrait, Xenium landscape)
M = np.array([[0.0, 2.0, 10.0], [-2.0, 0.0, 1000.0], [0, 0, 1]])


def he_to_um(p):
    q = M @ np.vstack([np.asarray(p, float).T, np.ones(len(p))])
    return q[:2].T * an.XENIUM_UM_PER_PX


def sq(x0, y0, x1, y1):
    return [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]


def test_round_trip():
    p = np.array([[3.0, 4.0], [100.0, 7.5]])
    assert np.allclose(an.um_to_he_px(he_to_um(p), M), p)


def test_label_rules(tmp_path):
    feats = [
        dict(type="Feature", geometry=dict(type="Polygon", coordinates=[sq(0, 0, 100, 100), sq(40, 40, 60, 60)]),
             properties=dict(classification=dict(name="Tumor"))),
        dict(type="Feature", geometry=dict(type="Polygon", coordinates=[sq(80, 80, 95, 95)]),
             properties=dict(classification=dict(name="Immune aggregate"))),
        dict(type="Feature", geometry=dict(type="MultiPolygon", coordinates=[[sq(0, 0, 10, 10)]]),
             properties=dict(classification="Exclude")),
    ]
    f = tmp_path / "a.geojson"
    f.write_text(json.dumps(dict(type="FeatureCollection", features=feats)))
    polys = an.read_geojson(f)
    cells_he = np.array([[20, 20], [50, 50], [90, 90], [5, 5], [200, 200]])
    lab = an.label_cells(he_to_um(cells_he), M, polys)
    assert list(lab) == ["Tumour", "", "Immune aggregate", "Exclude", ""]  # hole, smallest wins, Exclude wins


def test_unknown_class(tmp_path):
    f = tmp_path / "b.geojson"
    f.write_text(json.dumps([dict(geometry=dict(type="Polygon", coordinates=[sq(0, 0, 1, 1)]),
                                  properties=dict(classification=dict(name="Gleason 4")))]))
    with pytest.raises(ValueError):
        an.read_geojson(f)
    assert an.read_geojson(f, extra=("Gleason 4",))[0][0] == "Gleason 4"


def test_agreement():
    a = np.array(["Tumour"] * 300 + ["Tumour stroma"] * 300 + ["Exclude"] * 10 + [""] * 5, dtype=object)
    s, cons = an.agreement(a, a.copy())
    assert s["cohen_kappa"] == 1 and s["ari"] == 1 and s["n_labelled_both"] == 600 and s["enters_f2"]
    b = a.copy()
    b[:300] = "Tumour stroma"
    s, cons = an.agreement(a, b)
    assert s["n_agree"] == 300 and s["regions"] == ["Tumour stroma"] and not s["enters_f2"]
    assert (cons[:300] == "").all()
