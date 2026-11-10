"""One data root (2026-11-10, external review R1). Every input file is read under SF_DATA (default: <repo>/data, a
directory or a symlink to one). Layout and sources: DATA.md.
  raw/        downloaded datasets as published (one sub-directory per dataset)
  processed/  labelled h5ads made from raw/ by this repository's scripts
  derived/    inputs made by the investigators' earlier project (DATA.md lists each file, its SHA-256 and origin)
  pelka.h5ad  the CRC single-cell reference (Pelka et al. 2021)"""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("SF_DATA", REPO / "data"))
RAW = DATA / "raw"
PROC = DATA / "processed" / "p5" / "labeled_h5ad"
DERIVED = DATA / "derived"
PELKA = DATA / "pelka.h5ad"
