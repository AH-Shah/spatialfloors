import numpy as np
import pytest
from spatialfloors.data.genes import HGNC, filter_reason


@pytest.mark.skipif(not HGNC.exists(), reason="HGNC complete set not under SF_DATA (DATA.md)")
def test_filter_rules():
    g = ["NEGR1", "NEGATIVE1", "SystemControl3", "RPS6KA1", "RPL13", "RPSA", "MRPL12", "LINC00511", "MALAT1", "ACTBP2", "MT-CO1",
         "AC012345.1", "ENSG00000123", "KRT20", "HAND2-AS1"]
    r = dict(zip(g, filter_reason(g)))
    assert r["NEGR1"] == "" and r["RPS6KA1"] == "" and r["MRPL12"] == "" and r["MALAT1"] == "" and r["KRT20"] == ""
    assert r["NEGATIVE1"] == r["SystemControl3"] == "probe" and r["MT-CO1"] == "mito"
    assert r["RPL13"] == r["RPSA"] == "ribosomal" and r["ACTBP2"] == "pseudogene"
    assert r["LINC00511"] == r["HAND2-AS1"] == "unannotated_lncRNA" and r["AC012345.1"] == r["ENSG00000123"] == "unannotated"
    m = dict(zip(["Rpl13", "Rps6ka1", "Blank-3", "Gad1"], filter_reason(["Rpl13", "Rps6ka1", "Blank-3", "Gad1"], "mouse")))
    assert m == {"Rpl13": "ribosomal", "Rps6ka1": "", "Blank-3": "probe", "Gad1": ""}
