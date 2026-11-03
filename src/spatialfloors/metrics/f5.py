"""F5 label transfer (SPEC §F5): primary = held-out-marker AUROC (heldout_auroc, all units; decision 2026-10-28);
secondary = macro-F1 on IF-gated cells (GSE1 2-class gate)."""
import numpy as np
from sklearn.metrics import f1_score

CLASSES = ("epithelial", "immune_nonmyeloid", "myeloid")


def tails(v, lo_q, hi_q):
    """(negative, positive) masks from confident tails: <= lo_q quantile negative, >= hi_q quantile positive."""
    lo, hi = np.nanquantile(v, lo_q), np.nanquantile(v, hi_q); return v <= lo, v >= hi, (float(lo), float(hi))


def if_gate(panck, cd45, myelo, q=(.5, .8), myelo_q=(.5, .8)):
    """3-class gate on log1p IF means. epithelial = PanCK+ & CD45-; myeloid = CD45+ & myelo+ & PanCK-;
    immune non-myeloid = CD45+ & myelo- & PanCK-. Everything else ungated ('' label). Returns labels, thresholds."""
    pn, pp, tp = tails(np.log1p(panck), *q); cn, cp, tc = tails(np.log1p(cd45), *q); mn, mp, tm = tails(np.log1p(myelo), *myelo_q)
    lab = np.full(len(panck), "", dtype=object)
    lab[pp & cn] = "epithelial"; lab[cp & mp & pn] = "myeloid"; lab[cp & mn & pn] = "immune_nonmyeloid"
    return lab, dict(panck=tp, cd45=tc, myelo=tm, q=q, myelo_q=myelo_q)


def if_gate2(panck, myelo, q_panck=(.40, .60), q_myelo=.85):
    """GSE1 2-class fallback (manifests/if_gate_gse1_232.json decision): epithelial = PanCK >= q.60; myeloid =
    CD68_CK8_18 >= q.85 and PanCK <= q.40. Quantiles of the unit being scored (rule fixed, thresholds not re-tuned)."""
    pn, pp, tp = tails(np.log1p(panck), *q_panck); _, mp, tm = tails(np.log1p(myelo), .5, q_myelo)
    lab = np.full(len(panck), "", dtype=object); lab[pp] = "epithelial"; lab[mp & pn] = "myeloid"
    return lab, dict(panck=tp, myelo=tm, q_panck=q_panck, q_myelo=q_myelo)


def macro_f1_gated(gate, pred, classes=CLASSES):
    m = np.isin(gate, classes); return float(f1_score(gate[m], pred[m], labels=list(classes), average="macro", zero_division=0)), int(m.sum())


def unit_inputs(u, a):
    """Scoring inputs of unit u (AnnData a): counts over its filtered genes (library size of the held-out-marker score)
    and, on GSE1, the 2-class IF gate with its classes (else None, None)."""
    import scipy.sparse as sp
    from ..data.genes import filter_reason
    keep = filter_reason(a.var_names, "human") == ""
    gate, cls = (if_gate2(a.obs.IF_PanCK.values, a.obs.IF_CD68_CK8_18.values)[0], ("epithelial", "myeloid")) if u.startswith("gse1_") else (None, None)
    return sp.csr_matrix(a.X[:, keep]), a.var_names[keep], gate, cls


def score_classes(c, X, var_names, markers, gate=None, gate_classes=None):
    """Coarse class calls c -> held-out-marker AUROC (primary), per class, predicted shares, and macro-F1 on the IF-gated
    cells when a gate is given."""
    h, per = heldout_auroc(c, X, var_names, markers); r = dict(heldout_auroc=h, per_class=per, pred_share={k: round(float((c == k).mean()), 3) for k in (*CLASSES, "other")})
    if gate is not None: r["gated_macro_f1"], r["n_gated"] = macro_f1_gated(gate, c, gate_classes)
    return r


MIN_PRED = 50


def heldout_auroc(pred, X, var_names, markers, classes=CLASSES, min_pred=MIN_PRED):
    """Per class c: score_c = log1p CP10k of the summed held-out markers of c present on the platform (never model
    inputs); AUROC of score_c for cells predicted c vs every other cell (another F5 class or "other": a model cannot
    raise its AUROC by dumping hard cells into "other", review 2026-10-28). A class with < min_pred cells on either side
    scores .5 (no information). Returns (mean over classes, {c: dict(auroc, n_pred)}).
    X = counts over the unit's filtered genes (library size = row sum)."""
    import scipy.sparse as sp
    from sklearn.metrics import roc_auc_score
    X = sp.csr_matrix(X); lib = np.maximum(np.asarray(X.sum(1)).ravel(), 1); vn = list(map(str, var_names)); pred = np.asarray(pred).astype(str)
    out = {}
    for c in classes:
        gi = [vn.index(g) for g in markers[c] if g in vn]; sc = np.log1p(np.asarray(X[:, gi].sum(1)).ravel() / lib * 1e4)
        pos = pred == c; neg = ~pos
        a = float(roc_auc_score(np.r_[np.ones(pos.sum()), np.zeros(neg.sum())], np.r_[sc[pos], sc[neg]])) if min(pos.sum(), neg.sum()) >= min_pred and gi else .5
        out[c] = dict(auroc=a, n_pred=int(pos.sum()), n_markers=len(gi))
    return float(np.mean([v["auroc"] for v in out.values()])), out
