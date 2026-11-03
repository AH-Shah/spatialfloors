"""Commit stamp written into every output."""
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def commit():
    try:
        h = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--untracked-files=no"], capture_output=True, text=True).stdout.strip()
        return h + ("-dirty" if dirty else "")
    except Exception:
        return "nogit"


def stamp(d=None):
    d = dict(d or {}); d["spatialfloors_commit"] = commit(); return d
