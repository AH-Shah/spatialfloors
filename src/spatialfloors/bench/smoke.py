"""P0 MPS smoke tests on the GSE1 232 dev crop (50k cells x 5 epochs): s/epoch, peak MPS memory, finite gradients, NaN.
Usage: python -m spatialfloors.bench.smoke {scvi,scanvi,cellcharter,novae,nicheformer,gat} [--epochs 5]
Writes runs/p0_smoke/<model>.json. One model per process; the queue (scripts/p0_smoke_queue.sh) runs them in sequence."""
import argparse, json, os, time, threading
from pathlib import Path
import numpy as np

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import torch

OUT = Path("runs/p0_smoke"); CROP = OUT / "gse1_232_crop.h5ad"; PELKA = "/Users/Hemato/project_3/data/pelka.h5ad"


class MemPeak:
    """Polls torch.mps.driver_allocated_memory every 0.2 s (driver view: includes the cache) and process RSS."""
    def __init__(self): self.peak = 0; self.rss = 0; self._stop = False; self.t = threading.Thread(target=self._run, daemon=True)
    def _run(self):
        import resource
        while not self._stop:
            try: self.peak = max(self.peak, torch.mps.driver_allocated_memory())
            except Exception: pass
            self.rss = max(self.rss, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss); time.sleep(.2)
    def __enter__(self): self.t.start(); return self
    def __exit__(self, *a): self._stop = True; self.t.join()


def lightning_telemetry():
    """Lightning callback: per-epoch wall time, finite-gradient check every 20 steps, loss finiteness."""
    try: from lightning.pytorch.callbacks import Callback
    except ImportError: from pytorch_lightning.callbacks import Callback

    class Tele(Callback):
        def __init__(self): self.ep, self.t0, self.grad_ok, self.grad_checks, self.step = [], None, True, 0, 0
        def on_train_epoch_start(self, trainer, pl_module): self.t0 = time.time()
        def on_train_epoch_end(self, trainer, pl_module): self.ep.append(time.time() - self.t0)
        def on_after_backward(self, trainer, pl_module):
            self.step += 1
            if self.step % 20 == 1:
                self.grad_checks += 1
                for p in pl_module.parameters():
                    if p.grad is not None and not torch.isfinite(p.grad).all(): self.grad_ok = False; break
    return Tele()


def finite(x): x = np.asarray(x, dtype=float); return bool(np.isfinite(x).all())


def load_crop():
    import anndata as ad
    return ad.read_h5ad(CROP)


def run_scvi(a, epochs, n_latent=10, tag="scvi"):
    import scvi
    b = a[:, a.var.hvg.values].copy(); scvi.model.SCVI.setup_anndata(b)
    m = scvi.model.SCVI(b, n_latent=n_latent, gene_likelihood="nb"); tele = lightning_telemetry()
    m.train(max_epochs=epochs, accelerator="mps", devices=1, callbacks=[tele], early_stopping=False, enable_progress_bar=False)
    Z = m.get_latent_representation(); h = m.history.get("elbo_train")
    return dict(s_per_epoch=tele.ep, grad_finite=tele.grad_ok, grad_checks=tele.grad_checks, latent_finite=finite(Z), elbo_finite=finite(h.values) if h is not None else None,
                n_cells=b.n_obs, n_genes=b.n_vars), Z


def run_scanvi(a, epochs):
    """Joint scVI on a stratified Pelka subsample (<= 80k, <= 4k per ClusterMidway) + the crop, batch_key=modality, then SCANVI.from_scvi_model."""
    import scvi, anndata as ad, scipy.sparse as sp, pandas as pd
    P = ad.read_h5ad(PELKA, backed="r"); lab = P.obs.ClusterMidway.astype(str).values; rng = np.random.default_rng(0)
    take = np.sort(np.concatenate([rng.choice(np.flatnonzero(lab == c), min(4000, (lab == c).sum()), replace=False) for c in np.unique(lab)]))[:80000]
    sym = P.var.feature_name.astype(str).values; hv = a.var_names[a.var.hvg.values]; common = [g for g in hv if g in set(sym)]
    first = {}; [first.setdefault(g, i) for i, g in enumerate(sym)]; col = np.array([first[g] for g in common]); Xp = P.raw.X[take] if P.raw is not None else P.X[take]; Xp = sp.csr_matrix(Xp)[:, col]
    assert np.allclose(Xp.data[:1000], np.round(Xp.data[:1000])), "Pelka raw is not integer counts"
    p = ad.AnnData(Xp, obs=pd.DataFrame(dict(modality="dissociated", label=lab[take]), index=P.obs_names[take]), var=pd.DataFrame(index=common))
    s = a[:, common].copy(); s.obs = pd.DataFrame(dict(modality="spatial", label="Unknown"), index=a.obs_names)
    j = ad.concat([p, s]); j.obs_names_make_unique(); scvi.model.SCVI.setup_anndata(j, batch_key="modality", labels_key="label")
    m = scvi.model.SCVI(j, n_latent=10); t1 = lightning_telemetry()
    m.train(max_epochs=epochs, accelerator="mps", devices=1, callbacks=[t1], early_stopping=False, enable_progress_bar=False)
    m2 = scvi.model.SCANVI.from_scvi_model(m, unlabeled_category="Unknown"); t2 = lightning_telemetry()
    m2.train(max_epochs=epochs, accelerator="mps", devices=1, callbacks=[t2], early_stopping=False, enable_progress_bar=False)
    pred = m2.predict(j[j.obs.modality == "spatial"]); Z = m2.get_latent_representation()
    return dict(n_pelka=int(len(take)), n_spatial=s.n_obs, n_genes_common=len(common), n_hvg=len(hv), s_per_epoch_scvi=t1.ep, s_per_epoch_scanvi=t2.ep,
                grad_finite=t1.grad_ok and t2.grad_ok, latent_finite=finite(Z), pred_top=pd.Series(pred).value_counts().head(5).to_dict()), None


def run_cellcharter(a, epochs):
    """scVI representation on counts -> Delaunay + remove_long_links -> aggregate_neighbors(n_layers=3) -> cc.tl.Cluster (GMM) at k=2."""
    import cellcharter as cc, squidpy as sq, pandas as pd
    info, Z = run_scvi(a, epochs); a.obsm["X_scVI"] = Z; t = time.time()
    sq.gr.spatial_neighbors(a, coord_type="generic", delaunay=True, spatial_key="spatial"); cc.gr.remove_long_links(a)
    cc.gr.aggregate_neighbors(a, n_layers=3, use_rep="X_scVI", out_key="X_cellcharter"); t_agg = time.time() - t; t = time.time()
    k = int(a.obs.region.nunique()); g = cc.tl.Cluster(n_clusters=k, random_state=0, trainer_params=dict(accelerator="cpu", enable_progress_bar=False))
    g.fit(a, use_rep="X_cellcharter"); lab = np.asarray(g.predict(a, use_rep="X_cellcharter")); t_gmm = time.time() - t
    from sklearn.metrics import adjusted_rand_score
    return dict(info, aggregate_s=t_agg, gmm_s=t_gmm, k=k, dim=int(a.obsm["X_cellcharter"].shape[1]), cluster_sizes=pd.Series(lab).value_counts().tolist(),
                ari_smoke=float(adjusted_rand_score(a.obs.region.astype(str), lab))), None


def run_novae(a, epochs, steps=40, device="cpu"):
    """Raw counts, µm coordinates, locked 2000-HVG panel. MPS capped (PYTORCH_MPS_HIGH/LOW_WATERMARK_RATIO .5/.4 from the
    queue; set_per_process_memory_fraction rejects on this build) so it fails loudly instead of swapping. Zero-shot representations, then ONE fine-tune epoch capped at `steps` batches;
    s/epoch extrapolated from the median step x full batches/epoch. Known cost: serial dataloader (data tensors live on MPS, workers crash)."""
    import novae
    from sklearn.metrics import adjusted_rand_score
    a = a[:, a.var.hvg.values].copy(); a.var = a.var.drop(columns=["hvg"])
    m = novae.Novae.from_pretrained("prism-oncology/novae-human-0", revision="b8c0a5d7612bac6bc719ab57ed3cd16ad814728c"); voc = set(m.cell_embedder.gene_to_index)
    share = float(np.mean([g.lower() in voc for g in a.var_names]))
    novae.utils.spatial_neighbors(a, radius=80); t = time.time(); m.compute_representations(a, zero_shot=True, accelerator=device); t_zs = time.time() - t
    k = int(a.obs.region.nunique()); key0 = m.assign_domains(a, n_domains=k); ari0 = float(adjusted_rand_score(a.obs.region.astype(str), a.obs[key0].astype(str)))
    Z0 = np.asarray(a.obsm["novae_latent"]); tele = lightning_telemetry(); steps_rec = []
    try: from lightning.pytorch.callbacks import Callback
    except ImportError: from pytorch_lightning.callbacks import Callback

    class StepT(Callback):
        def on_train_start(s, tr, pl): s.n_full = len(tr.train_dataloader) if tr.train_dataloader is not None else None; s.t = None
        def on_train_batch_start(s, tr, pl, b, i): s.t = time.time()
        def on_train_batch_end(s, tr, pl, o, b, i):
            if device == "mps": torch.mps.synchronize()
            steps_rec.append(time.time() - s.t)
    st = StepT(); t = time.time()
    m.fine_tune(a, max_epochs=1, accelerator=device, callbacks=[tele, st], limit_train_batches=steps); t_ft = time.time() - t
    n_full = getattr(st, "n_full", None); step_med = float(np.median(steps_rec[3:] or steps_rec))
    m.compute_representations(a, accelerator=device); key1 = m.assign_domains(a, n_domains=k); Z = np.asarray(a.obsm["novae_latent"])
    sv = np.linalg.svd(Z - Z.mean(0), compute_uv=False); er = float(sv.sum() ** 2 / (sv ** 2).sum())
    sizes = a.obs[key1].value_counts(normalize=True)
    return dict(n_genes=int(a.n_vars), vocab_share_hvg=share, zero_shot_s=t_zs, ari_zero_shot=ari0, steps_run=len(steps_rec), batches_per_epoch=n_full,
                median_step_s=step_med, wall_per_step_s=t_ft / max(len(steps_rec), 1), s_per_epoch_extrapolated=(t_ft / max(len(steps_rec), 1)) * n_full if n_full else None,
                grad_finite=tele.grad_ok, latent_finite=finite(Z0) and finite(Z), effective_rank=er, max_domain_share=float(sizes.max()),
                ari_after_partial_epoch=float(adjusted_rand_score(a.obs.region.astype(str), a.obs[key1].astype(str))), device=device, device_reason="MPS footprint grows ~0.1 GB/5 s even with empty_cache and HIGH_WATERMARK .5 (12 GB guard kill twice); CPU flat 1.4-3.7 GB, same step time (profile runs/p0_smoke/novae_profile_*)"), None


def _nf_masked(core, ii, am, layer, cap=64, budget=4 * 1500 ** 2):
    """Mean over non-pad gene tokens (context tokens 0-2 dropped). Each sequence is cut to its own non-pad length: cells are
    length-sorted and run in microbatches of <= cap cells with cells x T^2 <= budget (half the 8 x 1,500 attention that fits).
    Non-pad outputs equal the full-length forward because padding keys are masked; only the pooling differs from the authors'
    get_embeddings, which averages the 1,497 positions including padding. layer: int, or a tuple -> {layer: embedding}
    from one forward pass (every requested layer read on the way up)."""
    Ls = tuple(layer) if isinstance(layer, (tuple, list)) else (layer,); want = {L % core.config.nlayers: L for L in Ls}
    n_tok = am.sum(1); n = len(n_tok); order = torch.argsort(n_tok); outs = {L: torch.empty(n, core.config.dim_model, device=ii.device) for L in Ls}
    assert bool((am == (torch.arange(am.shape[1], device=am.device)[None] < n_tok[:, None])).all()), "padding is not a suffix"
    i = 0
    while i < n:
        T = int(n_tok[order[min(n, i + cap) - 1]]); mb = max(1, min(cap, budget // T ** 2)); r = order[i:i + mb]
        T = min(am.shape[1], -(-int(n_tok[r].max()) // 128) * 128)  # 128-token buckets: few distinct shapes for the MPS allocator
        x = core.embeddings(ii[r, :T]) + core.positional_embedding(core.pos[:T].to(ii.device)); m = am[r, :T].bool()
        w = m[:, 3:, None].to(x.dtype)
        for li in range(max(want) + 1):
            x = core.encoder.layers[li](x, src_key_padding_mask=~m, is_causal=False)
            if li in want: outs[want[li]][r] = (x[:, 3:] * w).sum(1) / w.sum(1).clamp(min=1)
        i += mb; del x, w, m
        torch.mps.empty_cache()
    return outs if isinstance(layer, (tuple, list)) else outs[layer]


def run_nicheformer(a, epochs, batch=256, mbatch=8, layers=(-1,), n_cells=None, half=False, save_emb=None, pool="authors", ensg=None, specie="human", assay="cosmx", tech="cosmx", grad_check=True):
    """HF theislab/Nicheformer (pinned); CosMx technology mean. Genes mapped symbol -> Ensembl (Pelka CELLxGENE var) into the
    model's 20,310-gene order as CSR; each 256-cell batch is tokenized from CSR (the HF tokenizer densifies its input, so only
    256 x 20,310 is ever dense); one layer per pass. One MLM backward step checks finite gradients.
    pool: "authors" (get_embeddings, mean incl. padding), "masked" (_nf_masked, trimmed), "both" (keys "<L>" and "<L>_masked").
    ensg / specie / assay / tech: per-platform inputs (models.nicheformer.inputs); defaults = human CosMx."""
    import sys, anndata as ad, pandas as pd, scipy.sparse as sp
    hf = Path("external/hf_nicheformer").resolve(); sys.path.insert(0, str(hf))
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    model = AutoModelForMaskedLM.from_pretrained(str(hf), trust_remote_code=True).to("mps").eval()
    if half: model = model.half()  # P0 speed profile: fp16 inference (embeddings compared with fp32 on the same cells)
    tok = AutoTokenizer.from_pretrained(str(hf), trust_remote_code=True)
    tm = np.nan_to_num(np.load(f"external/nicheformer/data/model_means/{tech}_mean_script.npy"))  # authors' notebooks (tokenization/{cosmx,xenium,merfish}_*):
    tm = np.where(tm == 0, 1, np.where((tm % 1) >= .5, np.ceil(tm), np.floor(tm))); tok._load_technology_mean(tm)  # nan->0, round, 0->1 (raw file: 18,453 NaN -> every gene a token)
    ref = ad.read_h5ad(hf / "model.h5ad"); refg = pd.Index(ref.var_names.astype(str))
    if n_cells: a = a[:n_cells]
    if ensg is None:  # human symbols -> Ensembl via the Pelka CELLxGENE var; else models.nicheformer.gene_ids (orthologs for mouse)
        pv = ad.read_h5ad(PELKA, backed="r").var; m = dict(zip(pv.feature_name.astype(str), pv.index.astype(str))); ensg = [m.get(g) for g in a.var_names]
    col = np.array([refg.get_loc(e) if e in refg else -1 for e in ensg]); ok = col >= 0
    X = sp.csr_matrix(a.X, dtype=np.float32); lib = np.asarray(X.sum(1)).ravel(); cnt_share = float(X[:, ok].sum() / lib.sum())
    P = sp.csr_matrix((np.ones(ok.sum(), np.float32), (np.flatnonzero(ok), col[ok])), shape=(a.n_vars, len(refg)))
    B = sp.csr_matrix(X @ P)  # cells x model genes, sparse
    tok._load_reference_model = lambda: None  # we already aligned to model.h5ad gene order
    obs = pd.DataFrame(dict(modality="spatial", specie=specie, assay=assay), index=a.obs_names)
    out = {L: [] for L in layers}; outm = {L: [] for L in layers}; t_layer = {}; t_m = {L: 0. for L in layers}
    with torch.no_grad():
        for L in ([tuple(layers)] if pool == "masked" else layers):  # masked: all layers from one pass
            t = time.time()
            for s in range(0, B.shape[0], batch):
                inp = tok(ad.AnnData(B[s:s + batch], obs=obs.iloc[s:s + batch]))
                if pool in ("masked", "both"):
                    tm_ = time.time(); ii = torch.as_tensor(inp["input_ids"]).to("mps"); am = torch.as_tensor(inp["attention_mask"]).to("mps")
                    o_ = _nf_masked(model.nicheformer, ii, am, L)
                    for L_, v in (o_.items() if isinstance(L, tuple) else [(L, o_)]): outm[L_].append(v.float().cpu().numpy())
                    del ii, am, o_; torch.mps.empty_cache(); t_m[L if not isinstance(L, tuple) else L[0]] += time.time() - tm_
                for r in range(0, inp["input_ids"].shape[0] if pool != "masked" else 0, mbatch):  # 8 cells x 1,500 tokens per forward (256 needed a 34 GB attention buffer; 16 held 11 GB)
                    ii = torch.as_tensor(inp["input_ids"][r:r + mbatch]).to("mps"); am = torch.as_tensor(inp["attention_mask"][r:r + mbatch]).to("mps")
                    out[L].append(model.get_embeddings(input_ids=ii, attention_mask=am, layer=L, with_context=False).float().cpu().numpy())
                    del ii, am
                if (s // batch) % 20 == 0: torch.mps.empty_cache()
            torch.mps.synchronize(); t_layer[L if not isinstance(L, tuple) else L[0]] = time.time() - t
    E = {str(L): np.vstack(v) for L, v in out.items() if v} | {f"{L}_masked": np.vstack(v) for L, v in outm.items() if v}; torch.mps.empty_cache()
    if save_emb: np.savez(save_emb, **E)
    if half: model = model.float()
    gfin, mlm = None, float("nan")  # P3 zero-shot: no backward (it pushed gse0_S0 past the 14 GB cap after the embeddings were saved)
    if grad_check:
        model.train(); ii = torch.as_tensor(inp["input_ids"][:2]).to("mps"); am = torch.as_tensor(inp["attention_mask"][:2]).to("mps")
        o = model(input_ids=ii, attention_mask=am, labels=ii); o.loss.backward(); mlm = float(o.loss.item())
        gfin = bool(all(torch.isfinite(p.grad).all().item() for p in model.parameters() if p.grad is not None))
    return dict(n_cells=int(B.shape[0]), tokens_per_cell=float(inp["attention_mask"].sum(1).float().mean()), vocab_share_genes=float(ok.mean()), vocab_share_counts=cnt_share, n_genes_mapped=int(ok.sum()),
                tokens_nonpad_median=float(np.median(np.asarray(inp["attention_mask"]).sum(1))), half=half, mbatch=mbatch,
                s_per_layer_pass={str(L): v for L, v in t_layer.items()}, cells_per_s={str(L): B.shape[0] / v for L, v in t_layer.items()}, pool=pool,
                masked_s={str(L): v for L, v in t_m.items()}, masked_cells_per_s={str(L): B.shape[0] / v for L, v in t_m.items() if v},
                emb_dim=int(next(iter(E.values())).shape[1]), emb_finite=all(finite(v) for v in E.values()), mlm_loss=mlm,
                grad_finite=gfin, vocab_gate_passed=bool(ok.mean() >= .8)), None


def run_gat(a, epochs, hidden=128, layers=2):
    """F3-shaped GAT smoke: input = own Expr50 (log1p CP10k, TruncatedSVD-50 on the crop), graph = 15-NN within 4 quadrant
    blocks (no self loops), target = top-10 PCs of the receiver's own state; full batch on MPS, Adam, MSE."""
    import scipy.sparse as sp
    from torch_geometric.nn import GATConv
    from sklearn.decomposition import TruncatedSVD, PCA
    from sklearn.neighbors import NearestNeighbors
    from spatialfloors.metrics.f3 import quadrant_blocks
    X = a.X.astype(np.float32); lib = np.asarray(X.sum(1)).ravel(); X = sp.csr_matrix(sp.diags(1e4 / np.maximum(lib, 1)) @ X); X.data = np.log1p(X.data)
    F = TruncatedSVD(50, random_state=0).fit_transform(X).astype(np.float32); Y = PCA(10, random_state=0).fit_transform(X[:, a.var.hvg.values].toarray()).astype(np.float32)
    xy = a.obsm["spatial"]; bl, _ = quadrant_blocks(xy); nb = NearestNeighbors(n_neighbors=16).fit(xy).kneighbors(xy, return_distance=False)[:, 1:]
    src = nb.ravel(); dst = np.repeat(np.arange(len(xy)), 15); keep = bl[src] == bl[dst]
    ei = torch.as_tensor(np.vstack([src[keep], dst[keep]]), dtype=torch.long, device="mps")
    x = torch.as_tensor((F - F.mean(0)) / F.std(0), device="mps"); y = torch.as_tensor((Y - Y.mean(0)) / Y.std(0), device="mps")

    class Net(torch.nn.Module):
        def __init__(s):
            super().__init__(); d = [50] + [hidden] * layers
            s.c = torch.nn.ModuleList([GATConv(d[i], hidden // 4, heads=4, add_self_loops=False) for i in range(layers)]); s.o = torch.nn.Linear(hidden, 10)
        def forward(s, x, ei):
            for c in s.c: x = torch.nn.functional.elu(c(x, ei))
            return s.o(x)
    net = Net().to("mps"); opt = torch.optim.Adam(net.parameters(), 1e-3); ep, gok, losses = [], True, []
    for e in range(epochs):
        t = time.time(); opt.zero_grad(); loss = torch.nn.functional.mse_loss(net(x, ei), y); loss.backward()
        gok &= all(torch.isfinite(p.grad).all().item() for p in net.parameters() if p.grad is not None); opt.step(); torch.mps.synchronize()
        ep.append(time.time() - t); losses.append(float(loss.item()))
    return dict(s_per_epoch=ep, grad_finite=gok, loss=losses, loss_finite=finite(losses), n_edges=int(ei.shape[1]), cross_block_edges_dropped=int((~keep).sum())), None


RUN = dict(scvi=run_scvi, scanvi=run_scanvi, cellcharter=run_cellcharter, novae=run_novae, nicheformer=run_nicheformer, gat=run_gat)


def main():
    p = argparse.ArgumentParser(); p.add_argument("model", choices=list(RUN)); p.add_argument("--epochs", type=int, default=5)
    p.add_argument("--n-cells", type=int, default=None); p.add_argument("--layers", default="-1"); p.add_argument("--tag", default=None)
    p.add_argument("--half", action="store_true"); p.add_argument("--mbatch", type=int, default=8)
    p.add_argument("--pool", default="authors", choices=["authors", "masked", "both"]); a_ = p.parse_args()
    from spatialfloors.provenance import stamp
    torch.manual_seed(0); np.random.seed(0); a = load_crop(); t = time.time(); err = None
    with MemPeak() as mp:
        try:
            kw = dict(n_cells=a_.n_cells, layers=tuple(int(x) for x in a_.layers.split(",")), half=a_.half, mbatch=a_.mbatch, pool=a_.pool, save_emb=OUT / f"{a_.tag}_emb.npz" if a_.tag else None) if a_.model == "nicheformer" else {}
            info, _ = RUN[a_.model](a, a_.epochs, **kw)
        except Exception as e:
            import traceback; info, err = {}, traceback.format_exc()
    info.setdefault("n_cells", int(a.n_obs)); rep = stamp(dict(model=a_.model, epochs=a_.epochs, wall_s=round(time.time() - t, 1), peak_mps_driver_gb=round(mp.peak / 2**30, 2),
                     peak_rss_gb=round(mp.rss / 2**30, 2), torch=torch.__version__, ok=err is None, error=err, **info))
    (OUT / f"{a_.tag or a_.model}.json").write_text(json.dumps(rep, indent=1, default=float)); print(json.dumps(rep, indent=1, default=float)[:3000])
    try: torch.mps.empty_cache()
    except Exception: pass


if __name__ == "__main__":
    main()
