"""Analysis checks added after the external review (2026-11-10, R7).
1. Every tier-1 test in the committed report recomputes from its per-unit D (same strata, same margin).
2. A documented property of the F5 primary metric: calling fewer, more confident cells of a class can raise its
   held-out-marker AUROC (the reason for the matched-call-rate analysis, scripts/x_f5_matched.py)."""
import json
from pathlib import Path
import numpy as np
import pytest
import scipy.sparse as sp
from spatialfloors import stats
from spatialfloors.metrics.f5 import heldout_auroc

REPORT = Path(__file__).resolve().parents[1] / "results/reports/P5_corrected.json"
unit_stratum = lambda u: stats.stratum("gse303070_" if u.startswith("Donor") else u)


@pytest.mark.skipif(not REPORT.exists(), reason="tier-1 report not in this checkout")
def test_tier1_recomputes_from_unit_values():
    pairs = json.loads(REPORT.read_text())["pairs"]
    n = 0
    for k, v in pairs.items():
        if "|" not in k or "per_unit" not in v: continue
        units = sorted(v["per_unit"]); D = np.array([v["per_unit"][u] for u in units])
        t = stats.family_test(D, [unit_stratum(u) for u in units], margin=v["margin"])
        for f in ("est", "lower", "upper"): assert abs(t[f] - v[f]) < 1e-12, (k, f, t[f], v[f])
        assert t["noninferior"] == v["noninferior"]; n += 1
    assert n == 12


def test_f5_auroc_can_rise_when_fewer_cells_are_called():
    rng = np.random.default_rng(0); n = 2000
    marker = rng.gamma(2., 1., n)  # held-out marker signal of class "epithelial"
    X = sp.csr_matrix(np.column_stack([rng.poisson(marker * 3), rng.poisson(20., n)]))
    order = np.argsort(-marker)
    wide = np.full(n, "other", dtype=object); wide[order[:1000]] = "epithelial"
    narrow = np.full(n, "other", dtype=object); narrow[order[:300]] = "epithelial"
    mk = {"epithelial": ["M1"]}
    a_wide = heldout_auroc(wide, X, ["M1", "BG"], mk, classes=("epithelial",))[1]["epithelial"]["auroc"]
    a_narrow = heldout_auroc(narrow, X, ["M1", "BG"], mk, classes=("epithelial",))[1]["epithelial"]["auroc"]
    assert a_narrow > a_wide  # the metric rewards conservative calling; reported, not corrected (RESULTS.md, F5)
