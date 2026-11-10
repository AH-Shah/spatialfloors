"""Per-unit paired differences: D ~ 1 + S_centred, exact t, Holm, non-inferiority (PROTOCOL §4)."""

import numpy as np
from scipy import stats as st

MARGINS = dict(F1=0.03, F2=0.05, F3=0.03, F4=0.03, F5=0.05)  # caps; locked margins = margin_from_dev
REL_MARGIN = 0.10
STRATUM = {
    "gse1": "CosMx WTx",
    "gse0": "CosMx WTx",
    "gse303070": "CosMx 1k",
    "spatch": "Xenium 5K",
    "allen": "MERFISH mouse",
    "bruker": "CosMx WTx",
}  # bruker = Bruker WTx pan-cancer units


def stratum(crop):
    """Platform stratum of a crop name (strata = platform; both Allen datasets are MERFISH mouse brain)."""
    return next(v for k, v in STRATUM.items() if crop.startswith(k))


def aggregate_units(D_by_crop, unit_of_crop):
    """Per-crop paired differences -> one D per unit (patient / donor / animal): mean over that unit's crops
    (e.g. GSE303070 donor 280 = samples 280N + 280T; Allen animal 2 F4 = three sections). Returns (units, D, strata)."""
    import pandas as pd

    df = pd.DataFrame(dict(crop=list(D_by_crop), D=list(D_by_crop.values())))
    df["unit"] = df.crop.map(unit_of_crop)
    df["stratum"] = df.crop.map(stratum)
    assert df.unit.notna().all(), "crop without a unit"
    assert (df.groupby("unit").stratum.nunique() == 1).all(), "unit spans strata"
    g = df.groupby("unit").agg(D=("D", "mean"), stratum=("stratum", "first"))
    return g.index.tolist(), g.D.values, g.stratum.values


def margin_from_dev(family, floor_dev, chance_dev):
    """Locked margin = min(cap, REL_MARGIN x median over development units of the floor's above-chance score)
    (PROTOCOL §4). Chance: F1 permutation chance, F2 0, F3 0, F5 .5; F4 has no development units and keeps its cap."""
    if family == "F4":
        return MARGINS["F4"]
    above = np.median(np.asarray(floor_dev, float) - np.asarray(chance_dev, float))
    return float(min(MARGINS[family], REL_MARGIN * above))


def centred_intercept(D, strata=None):
    """OLS of D on 1 + mean-centred stratum indicators (one column per non-reference stratum).
    Returns intercept (= mean D), its SE, df (N - number of strata; N - 1 with one stratum)."""
    D = np.asarray(D, float)
    n = len(D)
    X = np.ones((n, 1))
    if strata is not None and len(np.unique(strata)) > 1:
        lv = np.unique(strata)
        S = (np.asarray(strata)[:, None] == lv[1:][None, :]).astype(float)
        X = np.hstack([X, S - S.mean(0)])
    b, *_ = np.linalg.lstsq(X, D, rcond=None)
    df = n - X.shape[1]
    r = D - X @ b
    s2 = r @ r / df
    se = np.sqrt(s2 * np.linalg.inv(X.T @ X)[0, 0])
    return float(b[0]), float(se), int(df)


def wilcoxon(x, alternative):
    x = np.asarray(x, float)
    return float(st.wilcoxon(x, alternative=alternative).pvalue) if len(x) > 1 and np.any(x != 0) else np.nan


def family_test(D, strata=None, margin=0.03, alpha=0.05):
    est, se, df = centred_intercept(D, strata)
    tq = st.t.ppf(1 - alpha, df)
    p_sup = float(st.t.sf(est / se, df)) if se > 0 else float(est <= 0)
    return dict(
        est=est,
        se=se,
        df=df,
        upper=est + tq * se,
        lower=est - tq * se,
        noninferior=bool(est + tq * se < margin),
        equivalent=bool(est + tq * se < margin and est - tq * se > -margin),
        floor_superior=bool(est + tq * se < 0),  # TOST = 90% CI in (-m, m)
        p_superior=p_sup,
        wilcoxon_superior_p=wilcoxon(D, "greater"),  # sensitivity: signed-rank, learned superior (D > 0)
        wilcoxon_noninferior_p=wilcoxon(np.asarray(D, float) - margin, "less"),
    )  # and floor non-inferior (D - margin < 0)


def holm(p, alpha=0.05):
    p = np.asarray(p, float)
    o = np.argsort(p)
    rej = np.zeros(len(p), bool)
    for i, j in enumerate(o):
        if p[j] <= alpha / (len(p) - i):
            rej[j] = True
        else:
            break
    return rej
