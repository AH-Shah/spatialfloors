"""F3 learned model: NCEM (Fischer et al. 2022, Nat Biotechnol; theislab/ncem 0.1.5, TF 2.15 in external/ncem_env).
The authors' ED-NCEM (EstimatorEDncem / ModelEDncem) and linear interactions NCEM at their MIBI-TOF settings (grid()),
inside the shared F3 design: per receiver and fold, one image per quadrant (graph within a block), node labels = one-hot
receiver ids (the label-free ids every F3 model uses), radius graph restricted by the F3 bleed rule (other lineages, same
FOV, same quadrant), target = the receiver's bleed-safe state PCs from block_cv (zeros elsewhere; those nodes are never
trained or scored). The encoder sees the neighbours' labels only, never the target. Training = the authors'
Estimator.train. The only data-driven choice is the radius (19 / 46 um), selected per species on development units.
score() also accepts the earlier tuned variants (type / type+ctx labels, higher capacity), reported on development only."""
import numpy as np
import scipy.sparse as sp

N_EVAL = 32
# The authors' MIBI-TOF setup (ncem_benchmarks grid_searches/mibitof_cancer; its val-best runs reproduced on their data in
# runs/p3/ncem_repro, HYPOTHESES row 375): Adam lr .05, all training images per step ('S'), 10 nodes per image, <= 2000
# epochs x 1 step, early stopping after 100 and LR halving after 50 epochs without validation-loss gain, 10 validation steps.
AUTH_TRAIN = dict(lr=.05, batch_size="S", n_eval=10, epochs=2000, patience=100, lr_patience=50,
                  train_kw=dict(lr_schedule_factor=.5, lr_schedule_min_lr=1e-10, max_validation_steps=10, shuffle_buffer_size=None))
# ED-NCEM GCN (ed_ncem_gcn.sh): latent 4, encoder 8 x {0, 1}, decoder 4 x 1, one conditional GCN layer of 8 with bias,
# dropout .05, l2 1e-6, own type as encoder condition; use_domain (one-hot image) is off: F3 test images are held-out
# quadrants, so their domain weights would never be trained.
AUTH_ED = dict(latent_dim=4, dropout_rate=.05, l2_coef=1e-6, l1_coef=0., enc_intermediate_dim=8, enc_depth=1, dec_intermediate_dim=4, dec_depth=1,
               cond_depth=1, cond_dim=8, cond_dropout_rate=0., cond_activation="relu", cond_l2_reg=1e-6, cond_use_bias=True)
# Our own type+ctx tuning (previous F3 recipe: lr 5e-3, 512 nodes per image, cond 64, latent 30) is not in the grid: it
# is not the authors' model (expression labels, our training) and lost to their setup on 278N (19 um .042 vs .044,
# 46 um .055 vs .058; runs/p3v3/f3_gse303070_148191A_278N_ncemA.json edT_*). Toys: scripts/diag_ncem_auth_tune.py.
# F3 port of the authors' sampling (scripts/diag_ncem_auth_fit.py, 278N r19): their step = 10 cells from each of 58
# similar-sized fields, every cell a target (~580 equally weighted cells); F3's FOV x quadrant images with receiver
# cells only weight cells by 1 / image count. One image per quadrant, 200 nodes per image, all validation nodes:
# linear interactions .033 -> .044 (input ridge .051), ED-NCEM .033 -> .042 (input ridge .050). 278N FOV-image runs
# kept as fovimg_* in runs/p3v3/f3_gse303070_148191A_278N_ncemA.json.
BAL = dict(image="block", n_eval=200, pad_eval=True)
RADII = (19., 46.)  # 278N, FOV images: 13 um lowest for every model and input ridge (CosMx: ~1-3 other-lineage cells)


def grid(radii=RADII):
    """Configs for score(): per radius the authors' ED-NCEM (their val-best MIBI config, enc depth 1; depth 0 never beat
    it on MIBI) and their linear interactions NCEM (l2 0, their MIBI val-best; l2 .1 tied on 278N, l2 1 lost .03 on MIBI).
    The only data-driven choice is the radius (the authors' free parameter, selected per species on development units)."""
    g = {}
    for R in radii:
        g[f"edA_e1_r{R:g}"] = dict(AUTH_TRAIN, **BAL, type_cond=True, radius_um=R, init_kw=dict(AUTH_ED))
        g[f"linA_l20_r{R:g}"] = dict(AUTH_TRAIN, **BAL, kind="interactions", radius_um=R)
    return g


def input_ridge(D, radius_um, mode="mean", min_cells=None, receivers=None):
    """Engineering floor of each NCEM variant: ridge (f3.ridge_fit_predict) on exactly its neighbour inputs on the same
    radius graph and folds: 'mean' = degree-scaled neighbour type labels (ED-NCEM's GCN input), 'presence' = binarised
    neighbour types (the interactions design for a fixed receiver type). Same output format as score()."""
    from ..floors import f3
    from ..metrics import f3 as m3
    A = radius_graph(D, radius_um); oh = sp.csr_matrix(np.eye(D["K"], dtype=np.float32)[D["ids"]]); C = (A @ oh).toarray()
    X = C / np.maximum(C.sum(1, keepdims=True), 1) if mode == "mean" else (C > 0).astype(np.float32); out = {}
    for r in (receivers if receivers is not None else f3.receivers(D, min_cells or f3.MIN_CELLS)):
        rows, Y, dep = f3.target(D, r); det = {}
        s, fo = m3.block_cv(X[rows], Y, D["xy"][rows], f3.ridge_fit_predict, f3.RF_UM, blocks=D["blocks"][rows], split=D["split"], covars=dep, detail=det)
        out[r] = (s, fo, np.nanmean(det.get("per_pc", [[np.nan]]), 0).tolist())
    return out


def _graph(D, nb, w, top=15):
    """Symmetrised binary adjacency of the floor kernel's top-15 other-lineage neighbours (w > 0 only)."""
    o = np.argsort(-w, 1)[:, :top]; ww = np.take_along_axis(w, o, 1); j = np.take_along_axis(nb, o, 1); i = np.repeat(np.arange(len(nb)), top)
    keep = ww.ravel() > 0; n = len(nb); A = sp.csr_matrix((np.ones(keep.sum(), np.float32), (i[keep], j.ravel()[keep])), shape=(n, n))
    A = ((A + A.T) > 0).astype(np.float32); return sp.csr_matrix(A)


def radius_graph(D, radius_um):
    """Binary radius graph (NCEM's get_data: squidpy spatial_neighbors, generic coordinates, no self loops) under the F3
    bleed rule: other-lineage cells of the same FOV and quadrant, within the 200 nearest (f3.kernel)."""
    from ..floors.f3 import kernel
    nb, w = kernel(D, radius_um); keep = w >= np.exp(-.5)  # w = exp(-d^2 / 2 R^2) on allowed pairs -> d <= R
    n = len(nb); i = np.repeat(np.arange(n), nb.shape[1])[keep.ravel()]
    A = sp.csr_matrix((np.ones(len(i), np.float32), (i, nb[keep])), shape=(n, n)); return sp.csr_matrix(((A + A.T) > 0).astype(np.float32))


def ctx_features(D, r, q, n_pc, seed=0):
    """Bleed-safe neighbour-state labels (as the GAT node features, models.gat.make): SVD of the genes receiver r does not
    dominate, fitted outside quadrant q, z-scored on those cells, zeroed for the receiver's lineage."""
    import scipy.sparse as sp_
    from sklearn.decomposition import TruncatedSVD
    from ..floors.bleedsafe import renorm_log
    g = np.flatnonzero(~D["mask"][r]); xg = renorm_log(D["x"][:, g]); k = min(n_pc, len(g) - 1)
    fit = np.flatnonzero(D["blocks"] != q); fit = np.sort(np.random.default_rng(seed).choice(fit, min(len(fit), 100000), replace=False))
    m1 = np.asarray(xg[fit].mean(0)).ravel(); sd = np.sqrt(np.maximum(np.asarray(xg[fit].multiply(xg[fit]).mean(0)).ravel() - m1 ** 2, 1e-6))
    xs = sp_.csr_matrix(xg @ sp_.diags(1 / sd)); E = TruncatedSVD(k, random_state=seed).fit(xs[fit]).transform(xs).astype(np.float32)
    E = (E - E[fit].mean(0)) / (E[fit].std(0) + 1e-6); E[D["lin"] == D["lin"][np.flatnonzero(D["ids"] == r)[0]]] = 0; return E


def make(D, r, rows, nb, w, cond_depth=1, latent_dim=10, seed=0, epochs=1000, lr=1e-4, batch_size=16, h0="type", n_pc=10, patience=20, cond_dim=8, l2_coef=0., dropout_rate=.1, n_eval=N_EVAL, pad_eval=False, lr_patience=5, init_kw=None,
         kind="ed_ncem", radius_um=None, type_cond=False, train_kw=None, image="fov"):
    """fit_predict for metrics.f3.block_cv (wants_rows): (Xtr, Ytr, Xte, tr, te, q) -> predicted PCs of rows[te].
    h0: node labels of the neighbours: 'type' (the authors' cell-type labels) or 'type+ctx' (plus the bleed-safe neighbour
    state, the information the floor and LightGBM get; NCEM's node_label_space 'standard' uses expression).
    kind: 'ed_ncem' (EstimatorEDncem, 'scaled' adjacency) or 'interactions' (EstimatorInteractions, the authors' linear
    interactions NCEM: own type x binarised neighbour types, 'full' adjacency; h0 'type' only). radius_um: NCEM's radius
    graph under the bleed rule (None: the floor kernel's top-15). type_cond: the encoder also sees the cell's own type
    (ED-NCEM use_type_cond; node covariate = own type, as the authors' use_covar_node_label). batch_size 'S' = all training
    images per step (the authors' 'S'). init_kw / train_kw override init_model / train arguments."""
    import tensorflow as tf
    from ncem.estimators import EstimatorEDncem, EstimatorInteractions
    from ncem.data import DataLoader
    from ..floors.f3 import groups
    assert kind == "ed_ncem" or h0 == "type"
    A = _graph(D, nb, w) if radius_um is None else radius_graph(D, radius_um); B = D["blocks"]; K = D["K"]; oh0 = np.eye(K, dtype=np.float32)[D["ids"]]
    tra = (lambda a: a) if kind == "interactions" else DataLoader._transform_a
    # images = FOV (or 500 um tile on single-FOV MERFISH sections) x quadrant: NCEM's own unit is a FOV-sized image, and an
    # epoch is len(train images) // batch_size steps (3 quadrant images gave 1 step per epoch and stopped after ~40 steps)
    # image='block': one image per quadrant. NCEM samples n_eval nodes per image per step, so FOV x quadrant images weight
    # each receiver cell by 1 / its image's receiver count (Kish n_eff 5-37% of n on 278N: 1-cell slivers weigh as much
    # as 460-cell FOVs); the authors' images are similar-sized fields in which every cell is a target. Same graph.
    img = B.astype(str) if image == "block" else np.char.add(np.char.add(B.astype(str), "_"), groups(D["fov"], D["xy"])); keys = list(np.unique(img))
    cells = {k: np.flatnonzero(img == k) for k in keys}; loc = np.full(len(B), -1); [loc.__setitem__(c, np.arange(len(c))) for c in cells.values()]
    Ab = {k: sp.csr_matrix(tra(sp.csr_matrix(A[c][:, c])) if len(c) > 1 else sp.csr_matrix((len(c), len(c)), dtype=np.float32))
          for k, c in cells.items()}  # one-cell image (FOV x quadrant sliver): NCEM's transform fails on a 1 x 1 matrix; no edges anyway
    kb = {k: int(k.split("_")[0]) for k in keys}

    def fp(Xtr, Ytr, Xte, tr, te, q):
        tf.keras.backend.clear_session(); tf.random.set_seed(seed); rng = np.random.default_rng(seed); np.random.seed(seed)
        e = EstimatorInteractions(log_transform=False) if kind == "interactions" else EstimatorEDncem(cond_type="gcn", use_type_cond=type_cond, log_transform=False)
        g_tr, g_te = rows[tr], rows[te]; npc = Ytr.shape[1]
        oh = oh0 if h0 == "type" else np.c_[oh0, ctx_features(D, r, q, n_pc, seed)]; K = oh.shape[1]
        h1 = {k: np.zeros((len(c), npc), np.float32) for k, c in cells.items()}
        for gi, y in zip(g_tr, Ytr): h1[img[gi]][loc[gi]] = y
        v = rng.random(len(g_tr)) < .1
        idx = lambda g, k: np.sort(loc[g[img[g] == k]]).astype(np.int64)
        e.a = Ab; e.h_0 = {k: oh[c] for k, c in cells.items()}; e.h_1 = h1; e.node_types = e.h_0
        e.node_type_names = {i: str(i) for i in range(K)}; e.n_features_type = K; e.n_features_standard = npc; e.node_feature_names = [f"pc{i}" for i in range(npc)]
        e.size_factors = {k: np.ones(len(c), np.float32) for k, c in cells.items()}; e.graph_covar_names = []; e.graph_covar = {k: np.array([], ndmin=1) for k in keys}
        e.node_covar = {k: oh0[c] if type_cond else np.empty((len(c), 0), np.float32) for k, c in cells.items()}; e.n_features_0, e.n_features_1 = K, npc
        e.n_graph_covariates, e.n_node_covariates = 0, (D["K"] if type_cond else 0)
        e.max_nodes = max(len(c) for c in cells.values()); e.complete_img_keys = keys; e.img_to_patient_dict = {k: "unit" for k in keys}; e.undefined_node_types = None
        e.domains = {k: i for i, k in enumerate(keys)}; e.n_domains = len(keys); e.cond_depth = cond_depth
        trk = [k for k in keys if kb[k] != q and (img[g_tr] == k).any()]; tek = [k for k in keys if (img[g_te] == k).any()]
        e.img_keys_train = trk; e.img_keys_eval = [k for k in trk if (img[g_tr[v]] == k).any()]; e.img_keys_test = tek
        e.nodes_idx_train = {k: idx(g_tr[~v], k) for k in trk}; e.nodes_idx_eval = {k: idx(g_tr[v], k) for k in e.img_keys_eval}
        if pad_eval:  # NCEM's eval generator takes whole chunks of n_eval nodes and drops the rest (an image with < n_eval validation nodes gives none)
            e.nodes_idx_eval = {k: np.resize(l, int(np.ceil(len(l) / n_eval) * n_eval)).astype(np.int64) for k, l in e.nodes_idx_eval.items()}
        lt = {k: idx(g_te, k) for k in tek}; padn = {k: int(np.ceil(len(l) / n_eval) * n_eval) for k, l in lt.items()}
        e.nodes_idx_test = {k: np.r_[l, np.resize(l, padn[k] - len(l))].astype(np.int64) for k, l in lt.items()}
        mk = dict(learning_rate=lr, n_eval_nodes_per_graph=n_eval, use_domain=False, scale_node_size=False, l2_coef=l2_coef)  # l2 0 / dropout .1: package defaults
        mk.update(dict(use_interactions=True, output_layer="linear") if kind == "interactions" else dict(latent_dim=latent_dim, cond_depth=cond_depth, cond_dim=cond_dim, output_layer="gaussian", dropout_rate=dropout_rate))
        e.init_model(**{**mk, **(init_kw or {})})
        tk = dict(epochs=epochs, batch_size=batch_size, validation_batch_size=16, max_steps_per_epoch=20, patience=patience, lr_schedule_patience=lr_patience, early_stopping=True, reduce_lr_plateau=True, log_dir=None)
        tk.update(train_kw or {}); tk["batch_size"] = len(trk) if tk["batch_size"] == "S" else tk["batch_size"]
        e.train(**tk)
        R = np.asarray(e.predict(batch_size=1)); R = R.reshape(-1, R.shape[-1])[:, :npc]  # gaussian layer: [mean, var]; chunks in image order
        pos, P, o = {}, np.zeros((len(g_te), npc), np.float32), 0
        for k in tek: pos.update({(k, int(i)): o + j for j, i in enumerate(lt[k])}); o += padn[k]
        for j, gi in enumerate(g_te): P[j] = R[pos[(img[gi], int(loc[gi]))]]
        h = e.history if isinstance(getattr(e, "history", None), dict) else {}
        fp.epochs.append(int(len(h["loss"])) if "loss" in h else -1)
        vl = h.get("val_loss", []); fp.val.append([float(vl[0]), float(np.min(vl))] if len(vl) else [])  # degeneracy check: validation loss must fall
        return P
    fp.wants_rows = True; fp.epochs = []; fp.val = []
    return fp


def score(D, ell_um=25., cond_depth=1, latent_dim=10, min_cells=None, seed=0, receivers=None, lr=1e-4, batch_size=16, h0="type", n_pc=10, **kw):
    """Per receiver block_cv score (same output format as floors.f3.score_floor)."""
    from ..floors import f3
    from ..metrics import f3 as m3
    nb, w = f3.kernel(D, ell_um); T = f3.id_mix(D, nb, w); out, ep = {}, {}
    for r in (receivers if receivers is not None else f3.receivers(D, min_cells or f3.MIN_CELLS)):
        rows, Y, dep = f3.target(D, r); fp = make(D, r, rows, nb, w, cond_depth, latent_dim, seed, lr=lr, batch_size=batch_size, h0=h0, n_pc=n_pc, **kw); det = {}
        s, fo = m3.block_cv(T[rows], Y, D["xy"][rows], fp, f3.RF_UM, blocks=D["blocks"][rows], split=D["split"], covars=dep, detail=det)  # X unused by NCEM
        out[r] = (s, fo, np.nanmean(det.get("per_pc", [[np.nan]]), 0).tolist()); pp = np.asarray(det.get("per_pc", [[np.nan]]), float)
        ep[r] = dict(epochs=fp.epochs, val=fp.val, nan_frac=float(np.isnan(pp).mean()), nan0=float(np.nan_to_num(pp).mean(1).mean())); print("ncem receiver", r, round(s, 4), fp.epochs, flush=True)
    return out, ep
