"""Round-2 loaders (data.units_v2) on tiny synthetic sections in each producer's layout."""

import json
import numpy as np, pandas as pd, scipy.sparse as sp
import pytest

h5py = pytest.importorskip("h5py")
from spatialfloors.data import units_v2  # noqa: E402


def _xenium(d):
    d.mkdir()
    # 3 features (gene, negative probe, deprecated codeword) x 4 cells
    M = np.array([[1, 0, 2, 0], [0, 1, 0, 0], [5, 5, 5, 5]])
    C = sp.csc_matrix(M)
    with h5py.File(d / "cell_feature_matrix.h5", "w") as f:
        m = f.create_group("matrix")
        m["data"], m["indices"], m["indptr"], m["shape"] = C.data, C.indices, C.indptr, np.array(M.shape)
        m["barcodes"] = np.array([b"c1", b"c2", b"c3", b"c4"])
        g = m.create_group("features")
        g["name"] = np.array([b"KRT8", b"NegControlProbe_00001", b"DeprecatedCodeword_0"])
        g["id"] = np.array([b"ENSG00000170421", b"NegControlProbe_00001", b"DeprecatedCodeword_0"])
        g["feature_type"] = np.array([b"Gene Expression", b"Negative Control Probe", b"Deprecated Codeword"])
    pd.DataFrame(dict(cell_id=["c4", "c3", "c2", "c1"], x_centroid=[4.0, 3, 2, 1], y_centroid=[40.0, 30, 20, 10])).to_parquet(d / "cells.parquet")  # fmt: skip
    (d / "gene_panel.json").write_text(json.dumps({}))


def _nsclc(d):
    d.mkdir()
    pd.DataFrame(dict(fov=[1, 1, 1, 2], cell_ID=[0, 1, 2, 1], EPCAM=[9, 3, 0, 1], NegPrb1=[9, 0, 1, 0])).to_csv(d / "S_exprMat_file.csv", index=False)  # fmt: skip
    pd.DataFrame({"fov": [2, 1, 1], "cell_ID": [1, 1, 2], "CenterX_global_px": [100.0, 0, 10], "CenterY_global_px": [0.0, 0, 50], "Mean.PanCK": [1.0, 2, 3], "Mean.CD45": [0.0, 0, 0], "Mean.CD3": [0.0, 0, 0]}).to_csv(d / "S_metadata_file.csv", index=False)  # fmt: skip


def test_xenium_and_nsclc(tmp_path, monkeypatch):
    _xenium(tmp_path / "x")
    _nsclc(tmp_path / "n")
    monkeypatch.setitem(units_v2.SECTIONS, "t_x", ("xenium", tmp_path / "x", "t_x", "s", "breast", "carcinoma"))
    monkeypatch.setitem(units_v2.SECTIONS, "t_n", ("nsclc", tmp_path / "n", "t_n", "s", "lung", "carcinoma"))

    a = units_v2.load("t_x")
    assert list(a.obs_names) == ["c4", "c3", "c2", "c1"]  # cells.parquet order
    assert list(a.var_names) == ["KRT8", "NegControlProbe_00001"]  # deprecated codewords dropped
    assert a.X.toarray().tolist() == [[0, 0], [2, 0], [0, 1], [1, 0]]
    assert a.obs.x_um.tolist() == [4, 3, 2, 1] and set(a.obs.fov) == {"0"}
    assert a.var.gene_id.iloc[0] == "ENSG00000170421"
    b = units_v2.load("t_x", ids=["c1", "c3"])
    assert b.X.toarray().tolist() == [[1, 0], [2, 0]]

    n = units_v2.load("t_n")
    assert list(n.obs_names) == ["t_n_2_1", "t_n_1_1", "t_n_1_2"]  # metadata order; cell_ID 0 dropped
    assert n.X.toarray().tolist() == [[1, 0], [3, 0], [0, 1]]
    assert np.allclose(n.obs.x_um, [18.0, 0, 1.8]) and np.allclose(n.obs.y_um, [0, 0, 9.0])
    assert n.obs.IF_PanCK.tolist() == [1.0, 2.0, 3.0]
    assert n.var.feature_type.tolist() == ["Gene Expression", "Negative Control Probe"]
    o = units_v2.obs("t_n")
    assert list(o.index) == list(n.obs_names) and o.fov.tolist() == ["2", "1", "1"]
