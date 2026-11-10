"""Nicheformer (F1, zero-shot): per-platform inputs and cell embeddings. Inputs follow the authors' tokenization
notebooks
(external/nicheformer/notebooks/tokenization): human CosMx -> gene symbols mapped to Ensembl ids via the Pelka CELLxGENE
var, technology mean "cosmx"; SPATCH Xenium -> the file's own Ensembl ids, technology mean "xenium", assay token
"Xenium"; Allen MERFISH -> mouse Ensembl ids mapped to human orthologs (Ensembl 116 BioMart; first row per mouse gene,
as
in the authors' merfish_mouse_brain notebook), technology mean "merfish", species token "mouse"."""

from pathlib import Path
import numpy as np
import pandas as pd

from spatialfloors.paths import RAW

ORTHO = RAW / "ensembl_orthologs/mouse_human_orthologs.tsv"
ORTHO_SHA = "38d7c752d021d5088ffb2ac4d4dd263fdfda65cfb429522c25ee71c49d022d6f"


def _mouse_ids(name, symbols):
    import anndata as ad
    from ..data.units import ALLEN, YAO

    if name.startswith("allen_yao"):
        v = ad.read_h5ad(YAO / "C57BL6J-638850-raw.h5ad", backed="r").var
        m = dict(zip(v.gene_symbol.astype(str), v.index.astype(str)))
    else:
        an = "animal2" if name.startswith("allen2_") else name.split("_")[1]
        v = ad.read_h5ad(ALLEN[an], backed="r").raw.var
        m = dict(zip(v.feature_name.astype(str), v.index.astype(str)))
    return [m.get(str(s)) for s in symbols]


def inputs(name, symbols):
    """dict(ensg=list or None, specie, assay, tech) for `embed` on crop `name` with var names `symbols`."""
    if name.startswith(("gse1_", "gse0_", "gse303070_", "bruker_")):
        return dict(ensg=None, specie="human", assay="cosmx", tech="cosmx")  # Bruker WTx = human CosMx
    if name.startswith("spatch"):
        import anndata as ad
        from ..data.units import RAW

        v = ad.read_h5ad(RAW / "spatch_xenium_5k/spatch_xenium_5k.h5ad", backed="r").var
        m = dict(zip(v.index.astype(str), v.gene_ids.astype(str)))
        return dict(ensg=[m.get(str(s)) for s in symbols], specie="human", assay="Xenium", tech="xenium")
    o = (
        pd.read_csv(ORTHO, sep="\t")
        .drop_duplicates("Gene stable ID")
        .set_index("Gene stable ID")["Human gene stable ID"]
    )
    mid = _mouse_ids(name, symbols)
    ensg = [o.get(x) if x is not None and isinstance(o.get(x), str) else None for x in mid]
    return dict(ensg=ensg, specie="mouse", assay="merfish", tech="merfish")


def _finite(x):
    return bool(np.isfinite(np.asarray(x, dtype=float)).all())


def _nf_masked(core, ii, am, layer, cap=64, budget=4 * 1500**2):
    """Mean over non-pad gene tokens (context tokens 0-2 dropped). Each sequence is cut to its own non-pad length: cells
    are
    length-sorted and run in microbatches of <= cap cells with cells x T^2 <= budget. Non-pad outputs equal the
    full-length
    forward because padding keys are masked; only the pooling differs from the authors' get_embeddings, which averages
    all
    1,497 positions including padding. layer: int, or a tuple -> {layer: embedding} from one forward pass."""
    import torch

    Ls = tuple(layer) if isinstance(layer, (tuple, list)) else (layer,)
    want = {L % core.config.nlayers: L for L in Ls}
    n_tok = am.sum(1)
    n = len(n_tok)
    order = torch.argsort(n_tok)
    outs = {L: torch.empty(n, core.config.dim_model, device=ii.device) for L in Ls}
    assert bool((am == (torch.arange(am.shape[1], device=am.device)[None] < n_tok[:, None])).all()), (
        "padding is not a suffix"
    )
    i = 0
    while i < n:
        T = int(n_tok[order[min(n, i + cap) - 1]])
        mb = max(1, min(cap, budget // T**2))
        r = order[i : i + mb]
        T = min(
            am.shape[1], -(-int(n_tok[r].max()) // 128) * 128
        )  # 128-token buckets: few distinct shapes for the MPS allocator
        x = core.embeddings(ii[r, :T]) + core.positional_embedding(core.pos[:T].to(ii.device))
        m = am[r, :T].bool()
        w = m[:, 3:, None].to(x.dtype)
        for li in range(max(want) + 1):
            x = core.encoder.layers[li](x, src_key_padding_mask=~m, is_causal=False)
            if li in want:
                outs[want[li]][r] = (x[:, 3:] * w).sum(1) / w.sum(1).clamp(min=1)
        i += mb
        del x, w, m
        torch.mps.empty_cache()
    return outs if isinstance(layer, (tuple, list)) else outs[layer]


def embed(
    a,
    epochs=0,
    batch=256,
    mbatch=8,
    layers=(-1,),
    n_cells=None,
    half=False,
    save_emb=None,
    pool="authors",
    ensg=None,
    specie="human",
    assay="cosmx",
    tech="cosmx",
    grad_check=True,
):
    """Zero-shot Nicheformer cell embeddings (Hugging Face theislab/Nicheformer, external/hf_nicheformer). Genes are
    mapped
    into the model's 20,310-gene order as a sparse matrix; each 256-cell batch is tokenized separately, so only
    256 x 20,310 values are ever dense. pool: "authors" (get_embeddings, mean including padding), "masked" (_nf_masked),
    "both". ensg / specie / assay / tech: per-platform inputs from `inputs`; defaults = human CosMx. grad_check: one
    masked
    language-model backward step to check finite gradients. Returns (dict of run statistics, None); embeddings are saved
    to save_emb (npz, one key per layer)."""
    import os, sys, time

    os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
    import anndata as ad, scipy.sparse as sp, torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    from spatialfloors.paths import PELKA

    hf = Path("external/hf_nicheformer").resolve()
    sys.path.insert(0, str(hf))
    model = AutoModelForMaskedLM.from_pretrained(str(hf), trust_remote_code=True).to("mps").eval()
    if half:
        model = model.half()
    tok = AutoTokenizer.from_pretrained(str(hf), trust_remote_code=True)
    tm = np.nan_to_num(
        np.load(f"external/nicheformer/data/model_means/{tech}_mean_script.npy")
    )  # technology means from the authors' tokenization notebooks
    tm = np.where(tm == 0, 1, np.where((tm % 1) >= 0.5, np.ceil(tm), np.floor(tm)))
    tok._load_technology_mean(tm)  # NaN -> 0, rounded, 0 -> 1, as the authors' notebooks
    ref = ad.read_h5ad(hf / "model.h5ad")
    refg = pd.Index(ref.var_names.astype(str))
    if n_cells:
        a = a[:n_cells]
    if ensg is None:  # human symbols -> Ensembl ids via the Pelka CELLxGENE var
        pv = ad.read_h5ad(str(PELKA), backed="r").var
        m = dict(zip(pv.feature_name.astype(str), pv.index.astype(str)))
        ensg = [m.get(g) for g in a.var_names]
    col = np.array([refg.get_loc(e) if e in refg else -1 for e in ensg])
    ok = col >= 0
    X = sp.csr_matrix(a.X, dtype=np.float32)
    lib = np.asarray(X.sum(1)).ravel()
    cnt_share = float(X[:, ok].sum() / lib.sum())
    P = sp.csr_matrix((np.ones(ok.sum(), np.float32), (np.flatnonzero(ok), col[ok])), shape=(a.n_vars, len(refg)))
    B = sp.csr_matrix(X @ P)  # cells x model genes, sparse
    tok._load_reference_model = lambda: None  # input already aligned to the model.h5ad gene order
    obs = pd.DataFrame(dict(modality="spatial", specie=specie, assay=assay), index=a.obs_names)
    out = {L: [] for L in layers}
    outm = {L: [] for L in layers}
    t_layer = {}
    t_m = {L: 0.0 for L in layers}
    with torch.no_grad():
        for L in [tuple(layers)] if pool == "masked" else layers:  # masked: every layer from one pass
            t = time.time()
            for s in range(0, B.shape[0], batch):
                inp = tok(ad.AnnData(B[s : s + batch], obs=obs.iloc[s : s + batch]))
                if pool in ("masked", "both"):
                    tm_ = time.time()
                    ii = torch.as_tensor(inp["input_ids"]).to("mps")
                    am = torch.as_tensor(inp["attention_mask"]).to("mps")
                    o_ = _nf_masked(model.nicheformer, ii, am, L)
                    for L_, v in o_.items() if isinstance(L, tuple) else [(L, o_)]:
                        outm[L_].append(v.float().cpu().numpy())
                    del ii, am, o_
                    torch.mps.empty_cache()
                    t_m[L if not isinstance(L, tuple) else L[0]] += time.time() - tm_
                for r in range(
                    0, inp["input_ids"].shape[0] if pool != "masked" else 0, mbatch
                ):  # 8 cells x 1,500 tokens per forward pass (memory)
                    ii = torch.as_tensor(inp["input_ids"][r : r + mbatch]).to("mps")
                    am = torch.as_tensor(inp["attention_mask"][r : r + mbatch]).to("mps")
                    out[L].append(
                        model.get_embeddings(input_ids=ii, attention_mask=am, layer=L, with_context=False)
                        .float()
                        .cpu()
                        .numpy()
                    )
                    del ii, am
                if (s // batch) % 20 == 0:
                    torch.mps.empty_cache()
            torch.mps.synchronize()
            t_layer[L if not isinstance(L, tuple) else L[0]] = time.time() - t
    E = {str(L): np.vstack(v) for L, v in out.items() if v} | {
        f"{L}_masked": np.vstack(v) for L, v in outm.items() if v
    }
    torch.mps.empty_cache()
    if save_emb:
        np.savez(save_emb, **E)
    if half:
        model = model.float()
    gfin, mlm = None, float("nan")
    if grad_check:
        model.train()
        ii = torch.as_tensor(inp["input_ids"][:2]).to("mps")
        am = torch.as_tensor(inp["attention_mask"][:2]).to("mps")
        o = model(input_ids=ii, attention_mask=am, labels=ii)
        o.loss.backward()
        mlm = float(o.loss.item())
        gfin = bool(all(torch.isfinite(p.grad).all().item() for p in model.parameters() if p.grad is not None))
    return dict(
        n_cells=int(B.shape[0]),
        tokens_per_cell=float(inp["attention_mask"].sum(1).float().mean()),
        vocab_share_genes=float(ok.mean()),
        vocab_share_counts=cnt_share,
        n_genes_mapped=int(ok.sum()),
        tokens_nonpad_median=float(np.median(np.asarray(inp["attention_mask"]).sum(1))),
        half=half,
        mbatch=mbatch,
        s_per_layer_pass={str(L): v for L, v in t_layer.items()},
        cells_per_s={str(L): B.shape[0] / v for L, v in t_layer.items()},
        pool=pool,
        masked_s={str(L): v for L, v in t_m.items()},
        masked_cells_per_s={str(L): B.shape[0] / v for L, v in t_m.items() if v},
        emb_dim=int(next(iter(E.values())).shape[1]),
        emb_finite=all(_finite(v) for v in E.values()),
        mlm_loss=mlm,
        grad_finite=gfin,
        vocab_gate_passed=bool(ok.mean() >= 0.8),
    ), None
