"""F5 floors (PLAN §3/§5): label transfer from the Pelka CRC atlas (ClusterMidway) to a crop, scored on IF-gated cells.
Input genes = the unit's panel minus every F5 held-out marker (P1, scripts/p1_panels.py), intersected with Pelka.
Pelka X is log1p-normalised (CELLxGENE): linear values = expm1(X).
- InSituType ML (insitutypeML, R env; /Users/crocodilus/code/label_transfer.R): reference profiles = mean linear
  expression per ClusterMidway over all annotated Pelka cells on the shared genes; neg = mean negative-probe count per
  cell, or 0.05 when the panel has none (that script's default).
- Logistic regression: multinomial, on the same stratified Pelka subsample as scANVI (<= 80k cells, <= 4k per
  ClusterMidway, seed 0); log1p-CP10k renormalised on the shared genes in both datasets, genes z-scored per dataset
  ("own") or with reference statistics ("ref"); C in {.1, 1} (grid = 4 configs, dev crops only).
Fine labels -> F5 classes by COARSE (locked); everything else is "other" (counts as a miss on gated cells)."""
import numpy as np
import scipy.sparse as sp

PELKA = "/Users/Hemato/project_3/data/pelka.h5ad"
COARSE = {**{k: "epithelial" for k in ("EpiT", "Epi")},
          **{k: "immune_nonmyeloid" for k in ("Plasma", "TCD4", "TCD8", "B", "Tgd", "TZBTB16", "NK", "ILC")},
          **{k: "myeloid" for k in ("Macro", "Mono", "DC", "Granulo")}}  # Mast, Endo, Fibro, Peri, SmoothMuscle, Schwann -> other


UNIT_REF = {"bruker_breast": "breast", "bruker_mobreast": "breast", "bruker_skin": "skin", "bruker_kidney": "kidney"}  # PLAN §17


def for_unit(u, mk_locked):
    """Reference and held-out markers of unit u: Pelka + the locked markers for CRC units (unchanged); for a new-cancer
    unit (PLAN §17) point reference() / pelka_genes() at its standardised reference (scripts/ref_prep.py) and extend
    COARSE with its fine -> class map. Returns the held-out marker dict."""
    global PELKA
    t = UNIT_REF.get(u)
    if t is None: return mk_locked
    import json
    from pathlib import Path
    r = json.loads((Path(__file__).resolve().parents[3] / f"manifests/refs/{t}.json").read_text()); assert r["qc_pass"], f"reference {t} failed QC"
    PELKA = r["ref"]; COARSE.update(r["coarse_map"]); return r["heldout_markers"]


def coarse(labels):
    return np.array([COARSE.get(str(l), "other") for l in labels], dtype=object)


def pelka_genes():
    import anndata as ad
    P = ad.read_h5ad(PELKA, backed="r"); g = P.var.feature_name.astype(str).values; P.file.close(); return g


def stratified_rows(labels, n_max=80000, per_class=4000, seed=0):
    rng = np.random.default_rng(seed); rows = []
    for c in np.unique(labels):
        i = np.flatnonzero(labels == c); rows.append(rng.choice(i, min(len(i), per_class), replace=False))
    rows = np.sort(np.concatenate(rows)); return np.sort(rng.choice(rows, n_max, replace=False)) if len(rows) > n_max else rows


def reference(genes, rows="stratified", chunk=20000, seed=0, counts=False):
    """(linear expression csr on `genes` in that order, ClusterMidway labels) for the stratified subsample, or
    rows="all" -> (per-label mean profile genes x labels, label names) over all annotated cells (chunked, backed).
    counts=True (stratified only): raw integer counts from Pelka .raw.X (scANVI / Tangram inputs; review 2026-10-28:
    X is log1p-normalised, so expm1(X) is not counts). Logistic regression is unaffected by the choice (it renormalises
    on the shared genes, which cancels any per-cell factor)."""
    import anndata as ad
    P = ad.read_h5ad(PELKA, backed="r"); lab = P.obs.ClusterMidway.astype(str).values; ok = ~np.isin(lab, ["nan", "", "undefined"])
    pg = P.var.feature_name.astype(str).values; first = {}
    for i, g in enumerate(pg): first.setdefault(g, i)
    col = np.array([first[g] for g in genes])
    if rows == "stratified":
        r = stratified_rows(np.where(ok, lab, "nan"), seed=seed); r = r[ok[r]]
        if counts:
            rg = P.raw.var["feature_name"].astype(str).values if "feature_name" in P.raw.var else P.var.feature_name.astype(str).values
            rf = {}
            for i, g in enumerate(rg): rf.setdefault(g, i)
            X = sp.csr_matrix(P.raw.X[r])[:, np.array([rf[g] for g in genes])]; P.file.close()
            assert np.allclose(X.data, np.round(X.data)), "Pelka raw is not integer counts"; return X.astype(np.float32), lab[r]
        X = sp.csr_matrix(P.X[r])[:, col]; X.data = np.expm1(X.data); P.file.close(); return X.astype(np.float32), lab[r]
    names = np.unique(lab[ok]); S = np.zeros((len(names), len(col))); n = np.zeros(len(names)); li = {c: i for i, c in enumerate(names)}
    for a in range(0, P.n_obs, chunk):
        X = sp.csr_matrix(P.X[a:a + chunk])[:, col]; X.data = np.expm1(X.data); m = ok[a:a + chunk]; idx = np.array([li.get(c, -1) for c in lab[a:a + chunk]])
        M = sp.csr_matrix((np.ones(m.sum()), (idx[m], np.flatnonzero(m))), shape=(len(names), X.shape[0])); S += (M @ X).toarray(); n += np.asarray(M.sum(1)).ravel()
    P.file.close(); return (S / n[:, None]).T, names


def _lognorm(X):
    X = sp.csr_matrix(X, dtype=np.float32); s = np.asarray(X.sum(1)).ravel(); Y = sp.csr_matrix(sp.diags(1e4 / np.maximum(s, 1e-9)) @ X); Y.data = np.log1p(Y.data); return Y.toarray()


def logreg(Xref, yref, Xq, C=1., scaling="own", seed=0):
    """Multinomial LR on reference fine labels; returns predicted fine labels for the query."""
    from sklearn.linear_model import LogisticRegression
    R, Q = _lognorm(Xref), _lognorm(Xq); mr, sr = R.mean(0), R.std(0) + 1e-6
    mq, sq = (Q.mean(0), Q.std(0) + 1e-6) if scaling == "own" else (mr, sr)
    m = LogisticRegression(C=C, max_iter=2000, random_state=seed).fit((R - mr) / sr, yref)
    return m.predict((Q - mq) / sq)


def spatial_prior(logliks, names, xy_um, k=15, beta=1., sweeps=5):
    """F5 floor v3 (PLAN §19): mean-field neighbour prior on InSituType's log-likelihoods (cells x types). The
    log-likelihoods are divided by their median top-1 - top-2 margin (so beta counts margins), then for `sweeps` rounds
    z = ll + beta * (mean neighbour posterior over the k nearest cells); labels = argmax. No parameter is fitted."""
    from scipy.spatial import cKDTree
    ll = np.asarray(logliks, float); s2 = np.sort(ll, 1); z0 = ll / max(float(np.median(s2[:, -1] - s2[:, -2])), 1e-6)
    nb = cKDTree(xy_um).query(xy_um, k + 1)[1][:, 1:]; R = np.exp(z0 - z0.max(1, keepdims=True)); R /= R.sum(1, keepdims=True)
    for _ in range(sweeps):
        z = z0 + beta * R[nb].mean(1); z -= z.max(1, keepdims=True); R = np.exp(z); R /= R.sum(1, keepdims=True)
    return np.asarray(names)[R.argmax(1)].astype(str)


def insitutype(Xq, profiles, names, neg=None, workdir=None):
    """insitutypeML via Rscript (scripts/insitutype_ml.R, pixi env r) -> (fine labels, posterior prob)."""
    import subprocess, tempfile, scipy.io as sio, pandas as pd
    from pathlib import Path
    d = Path(workdir or tempfile.mkdtemp()); d.mkdir(parents=True, exist_ok=True)
    from . import lean
    if lean.on():  # whole-slide arm (SF_LEAN=1): the same matrix as raw CSC arrays (scripts/insitutype_ml.R reads them), not Matrix Market text
        C = sp.csc_matrix(Xq, dtype=np.float64); C.sort_indices(); assert C.nnz < 2 ** 31
        C.indices.astype(np.int32).tofile(d / "x_i.bin"); C.indptr.astype(np.int32).tofile(d / "x_p.bin"); C.data.tofile(d / "x_x.bin"); (d / "x_dim.txt").write_text(f"{C.shape[0]} {C.shape[1]} {C.nnz}")
    else: sio.mmwrite(str(d / "x.mtx"), sp.csr_matrix(Xq).astype(np.float64))
    pd.DataFrame(profiles, columns=names).to_csv(d / "profiles.csv", index=False)
    np.savetxt(d / "neg.txt", np.full(Xq.shape[0], .05) if neg is None else np.asarray(neg, float))
    root = Path(__file__).resolve().parents[3]
    subprocess.run(["pixi", "run", "-e", "r", "Rscript", str(root / "scripts/insitutype_ml.R"), str(d)], check=True, cwd=root)
    o = pd.read_csv(d / "out.csv"); return o.clust.astype(str).values, o.prob.values
