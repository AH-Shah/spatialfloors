"""F3 learned model: GAT (PyG GATConv) predicting a receiver's state PCs from its niche, inside the shared F3 design
(floors.f3: same receivers, targets, quadrant CV, buffer, kernel graph). Inputs are bleed-safe by the floor's rules:
node j carries [one-hot id_j, the receiver id's context SVD of j (genes the receiver does not dominate; SVD fitted
outside the held-out quadrant)], zeroed for cells of the receiver's lineage (its own and same-lineage expression never
enter); edges = the floor kernel's top-15 neighbours (same FOV, same quadrant, other lineages), into every cell.
Head: MLP on [GAT output, the floor's features T + ctx] plus a linear map of the floor features (the model contains
its floor; engineering rule: a learned model must contain and beat its trivial baseline).
Training (per receiver, per fold): full graph of the training quadrants, MSE on the receivers' training PCs, Adam lr 1e-3,
weight decay 1e-4, <= 300 epochs, MLP head (hidden, ELU) after message passing, early stopping (patience 20) on a random 15% of training receivers; seed 0.
Grid (PLAN §6, amended 2026-10-28): attention {GATConv, GATv2Conv} x layers {1, 2}, hidden 128."""
import json, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import torch
from ..floors import f3
import torch.utils.checkpoint
CHUNK_ABOVE = 1_000_000  # edges; every crop (<= 51.5k cells, <= ~0.8M edges) keeps the plain path, bit-identical to development
CHUNK_T = 40_000  # targets per padded chunk
CKPT_EDGES = 250_000  # edges per chunk; 2026-10-31 crash fix: w3g_231 (298k-cell whole slide, ~3.4M edges) hit 19-21 GB on MPS


def hops(ei, idx, layers, n):
    """Whole slides only (2026-11-01 speed fix: w3g_231 ran ~15 s per epoch over the full 3.4M-edge graph): the outputs at
    idx depend only on their `layers`-hop in-neighbourhood, so run on that subgraph (exact). Crops keep the full graph."""
    if ei.shape[1] <= CHUNK_ABOVE: return None, ei, idx
    from torch_geometric.utils import k_hop_subgraph
    sub, e, inv, _ = k_hop_subgraph(idx, layers, ei, relabel_nodes=True, num_nodes=n, flow="source_to_target"); return sub, e, inv


def edges(D, nb, w, top=15):
    """Directed edges j -> i for the top-`top` kernel weights of each i (w > 0 only)."""
    o = np.argsort(-w, 1)[:, :top]; ww = np.take_along_axis(w, o, 1); j = np.take_along_axis(nb, o, 1)
    i = np.repeat(np.arange(len(nb)), top); keep = ww.ravel() > 0
    return np.vstack([j.ravel()[keep], i[keep]])


class _Sync(torch.autograd.Function):
    """Identity; its backward synchronises MPS (frees the chunk's temporaries)."""
    @staticmethod
    def forward(ctx, x): return x.view_as(x)
    @staticmethod
    def backward(ctx, g):
        if g.is_mps: torch.mps.synchronize()
        return g


class Net(torch.nn.Module):
    """GAT message passing (no self loops: a receiver's h is built from its neighbours only) -> MLP head on
    [h_i, floor features of i]; the linear part (lin) on the floor features alone means the model contains its floor."""
    def __init__(s, d_in, d_floor, hidden, layers, d_out, conv="gat"):
        from torch_geometric.nn import GATConv, GATv2Conv
        super().__init__(); d = [d_in] + [hidden] * layers; Conv = GATv2Conv if conv == "gatv2" else GATConv
        s.c = torch.nn.ModuleList([Conv(d[k], hidden // 4, heads=4, add_self_loops=False) for k in range(layers)])
        s.o = torch.nn.Sequential(torch.nn.Linear(hidden + d_floor, hidden), torch.nn.ELU(), torch.nn.Linear(hidden, d_out)); s.lin = torch.nn.Linear(d_floor, d_out)

    def forward(s, x, ei, idx, f):
        h = x
        for c in s.c: h = torch.nn.functional.elu(c(h, ei) if ei.shape[1] <= CHUNK_ABOVE and not getattr(s, "force_chunk", False) else s.chunked(c, h, ei))
        return s.o(torch.cat([h[idx], f], 1)) + s.lin(f)

    @staticmethod
    def chunked(c, h, ei):
        """Same layer on whole slides: targets in chunks (softmax is per target, so exact), each chunk one bipartite call
        (sources -> targets, shared lin) recomputed in backward (checkpoint). Every chunk is padded to one shape
        (CKPT_EDGES edges, CHUNK_T targets; padding edges go to a dummy target that is dropped): MPS caches a graph per
        shape. MPS frees a temporary only once the queued command buffer has run, so one unsynchronised backward over all
        chunks held ~0.35 GB per chunk (2026-10-31 profile: 4.5M edges, tensors < 1.1 GB, 12.8 GB driver): sync per chunk."""
        o = torch.argsort(ei[1]); ei = ei[:, o]; n = h.shape[0]; cnt = torch.bincount(ei[1], minlength=n).cumsum(0).cpu()
        out, lo, E, T = [], 0, CKPT_EDGES, CHUNK_T
        while lo < n:
            e0 = int(cnt[lo - 1]) if lo else 0; hi = min(int(torch.searchsorted(cnt, e0 + E, right=True)), lo + T, n); hi = max(hi, lo + 1)
            e1 = int(cnt[hi - 1]); m = e1 - e0; assert m <= E, "a target with more than CKPT_EDGES in-edges"
            src = torch.zeros(E, dtype=ei.dtype, device=ei.device); src[:m] = ei[0, e0:e1]
            dst = torch.full((E,), T, dtype=ei.dtype, device=ei.device); dst[:m] = ei[1, e0:e1] - lo
            loc = torch.stack([torch.arange(E, device=ei.device), dst]); k = hi - lo
            def f_(h_, src_=src, loc_=loc, a=lo, k_=k):
                xd = torch.cat([h_[a:a + k_], h_.new_zeros(T + 1 - k_, h_.shape[1])])
                return c((h_[src_], xd), loc_, size=(E, T + 1))[:k_]
            hs = _Sync.apply(h) if h.requires_grad else h  # backward: wait for this chunk's GPU work before the next one
            out.append(_Sync.apply(torch.utils.checkpoint.checkpoint(f_, hs, use_reentrant=False)) if torch.is_grad_enabled() else f_(hs)); lo = hi
            if h.is_mps: torch.mps.synchronize()
        return torch.cat(out)


def make(D, r, rows, nb, w, layers=2, hidden=128, n_pc=50, device="mps", seed=0, max_epochs=300, patience=20, lr=1e-3, conv="gat"):
    """(X callable, fit_predict) for f3 block_cv on receiver r."""
    ctx_full = {}
    ei_all = edges(D, nb, w); onehot = np.eye(D["K"], dtype=np.float32)[D["ids"]]; own = D["lin"] == D["lin"][rows[0]]
    from ..floors.bleedsafe import renorm_log
    import scipy.sparse as sp
    from sklearn.decomposition import TruncatedSVD
    from ..floors.bleedsafe import kernel_mean, group_center
    T = f3.id_mix(D, nb, w); g = np.flatnonzero(~D["mask"][r]); xg = renorm_log(D["x"][:, g]); k = min(n_pc, len(g) - 1)

    def X(q, tr_local=None):
        f_ = pre.pop(q, None); return f_.result() if f_ is not None else X_(q)

    def X_(q):
        fit = np.flatnonzero(D["blocks"] != q); fit = np.sort(np.random.default_rng(seed).choice(fit, min(len(fit), 100000), replace=False))
        m1 = np.asarray(xg[fit].mean(0)).ravel(); sd = np.sqrt(np.maximum(np.asarray(xg[fit].multiply(xg[fit]).mean(0)).ravel() - m1 ** 2, 1e-6))
        xs = sp.csr_matrix(xg @ sp.diags(1 / sd)); E = TruncatedSVD(k, random_state=seed).fit(xs[fit]).transform(xs).astype(np.float32)
        floor = np.c_[T[rows], kernel_mean(nb, w, group_center(E, D["fov"]).astype(np.float32), rows)]  # = floors.f3 [T, ctx]
        E = (E - E[fit].mean(0)) / (E[fit].std(0) + 1e-6); E[own] = 0; ctx_full[q] = np.c_[onehot, E]
        return floor  # block_cv scales it on training rows; the head reads it next to the graph

    def fp(Xtr, Ytr, Xte, tr, te, q):
        torch.manual_seed(seed); rng = np.random.default_rng(seed); feats = torch.as_tensor(ctx_full[q], device=device)
        keep = (D["blocks"][ei_all[0]] != q) & (D["blocks"][ei_all[1]] != q); ei_tr = torch.as_tensor(ei_all[:, keep], device=device)
        ei_te = torch.as_tensor(ei_all[:, (D["blocks"][ei_all[1]] == q) & (D["blocks"][ei_all[0]] == q)], device=device)
        ftr, fte = torch.as_tensor(Xtr, dtype=torch.float32, device=device), torch.as_tensor(Xte, dtype=torch.float32, device=device)
        net = Net(feats.shape[1], Xtr.shape[1], hidden, layers, Ytr.shape[1], conv).to(device); opt = torch.optim.Adam(net.parameters(), lr, weight_decay=1e-4)
        y = torch.as_tensor(Ytr, dtype=torch.float32, device=device); gi = torch.as_tensor(rows[tr], device=device)
        v = rng.random(len(tr)) < .15; vt, vv = torch.as_tensor(np.flatnonzero(~v), device=device), torch.as_tensor(np.flatnonzero(v), device=device)
        sub, ei_tr, gi = hops(ei_tr, gi, layers, len(feats)); ftrain = feats if sub is None else feats[sub]
        net.force_chunk = sub is not None  # 2026-11-01 crash fix: plain path on <= 1M-edge subgraphs spiked to 12-17 GB (MPS); chunked 1.5 G
        best, bad, state, t0 = np.inf, 0, None, time.time()
        for ep in range(max_epochs):
            net.train(); opt.zero_grad(); out = net(ftrain, ei_tr, gi, ftr); loss = ((out[vt] - y[vt]) ** 2).mean(); loss.backward(); opt.step()
            net.eval()
            with torch.no_grad(): lv = float(((net(ftrain, ei_tr, gi, ftr)[vv] - y[vv]) ** 2).mean())
            if lv < best - 1e-5: best, bad, state = lv, 0, {kk: vv_.detach().clone() for kk, vv_ in net.state_dict().items()}
            else:
                bad += 1
                if bad >= patience: break
        if state is not None: net.load_state_dict(state)  # None only if validation loss was never finite
        net.eval(); fp.epochs.append(ep + 1); print(f"gat receiver {r} fold {q}: {ep + 1} epochs, {time.time() - t0:.0f} s, {ei_tr.shape[1]} train edges", flush=True)
        sub, ei_te, gt = hops(ei_te, torch.as_tensor(rows[te], device=device), layers, len(feats))
        with torch.no_grad(): pred = net(feats if sub is None else feats[sub], ei_te, gt, fte).cpu().numpy()
        torch.mps.empty_cache(); return pred
    fp.wants_rows = True; fp.epochs = []
    # 2026-11-01 (user: no GPU idle during CPU prep): every fold's features (SVD, kernel means) are built on a CPU thread
    # while the GPU trains the earlier folds; same function, same seeds, so identical features.
    ex = ThreadPoolExecutor(1); pre = {q: ex.submit(X_, q) for q in np.unique(D["blocks"][rows])}; ex.shutdown(wait=False)
    return X, fp


def score(D, ell_um=25., layers=2, hidden=128, n_pc=50, min_cells=f3.MIN_CELLS, device="mps", seed=0, receivers=None, full=False, ckpt=None, **train):  # train: lr, patience, max_epochs, conv
    """Per receiver block_cv score of the GAT (same output format as floors.f3.score_floor). ckpt: json path; finished
    receivers are saved there and skipped on a rerun (2026-11-01: w3g_231 lost 6.4 h to a swap kill with nothing saved)."""
    from ..metrics import f3 as m3
    nb, w = f3.kernel(D, ell_um); out = {}; ep = {}
    done = json.loads(Path(ckpt).read_text()) if ckpt and not full and Path(ckpt).exists() else {}
    todo = [r for r in (receivers if receivers is not None else f3.receivers(D, min_cells)) if str(r) not in done]
    for r in done: print(f"gat receiver {r}: resumed from {ckpt}", flush=True)
    out.update({int(r): tuple(v[:3]) for r, v in done.items()}); ep.update({int(r): v[3] for r, v in done.items()})
    prep = lambda r: (lambda t: (t, make(D, r, t[0], nb, w, layers, hidden, n_pc, device, seed, **train)))(f3.target(D, r))
    ex = ThreadPoolExecutor(1); nxt = ex.submit(prep, todo[0]) if todo else None
    for i, r in enumerate(todo):
        (rows, Y, dep), (X, fp) = nxt.result(); det = {"full": full}
        nxt = ex.submit(prep, todo[i + 1]) if i + 1 < len(todo) else None  # next receiver's CPU prep overlaps this one's GPU training
        s, fo = m3.block_cv(X, Y, D["xy"][rows], fp, f3.RF_UM, blocks=D["blocks"][rows], split=D["split"], covars=dep, detail=det)
        out[r] = (s, fo, np.nanmean(det.get("per_pc", [[np.nan]]), 0).tolist()) + ((det,) if full else ()); ep[r] = fp.epochs
        if ckpt and not full: done[str(r)] = [*out[r], ep[r]]; Path(ckpt).write_text(json.dumps(done, default=float))
    return out, ep
