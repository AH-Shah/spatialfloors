"""Bit-identity of the vendored L6 ladder vs row 323 (crc_232). ~minutes on CPU; needs GSE1 raw data."""
import os, pytest
from pathlib import Path

RAW = Path("/Users/Hemato/project_3/data/raw/GSE1/cosmx_wt_232.h5ad")


@pytest.mark.slow
@pytest.mark.data
@pytest.mark.skipif(not RAW.exists() or not os.environ.get("SF_SLOW"), reason="set SF_SLOW=1 with GSE1 on disk")
def test_row323_bit_identical(tmp_path):
    from spatialfloors.cli import main
    assert main(["check", "--out", str(tmp_path)]) == 0
