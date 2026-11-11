"""Exploratory: is the registered F3 LightGBM held back by its fixed 100 rounds or by its extra features?
LightGBM with early stopping (leaves 15, learning rate .05, up to 1,000 rounds, patience 30; validation = about 15% of
the training FOVs or 500 um tiles, floors.f3.groups; the recipe of the development residual variant) on
  floor_feats  exactly the v3 floor's features (floors.f3.v3_floor_parts)
  all_feats    every v3 feature group, as the registered LightGBM (floors.f3.v3_parts)
Same receivers, targets, folds and score as the registered runs (floors.f3.score_v3).
Usage (unit: env scvi, which has LightGBM; report: env floors; cwd = the runs root):
  python scripts/x_f3_lgbm_es.py unit <crop> <out dir>      -> <out dir>/f3es_<crop>.json
  python scripts/x_f3_lgbm_es.py report <out dir> <p4 dir>  -> <out dir>/X_f3_lgbm_es.{json,md}"""

import json, sys, time
from pathlib import Path
import numpy as np

PUB = Path(__file__).resolve().parents[1]


def fit_predict_es(grp, stops):
    from lightgbm import LGBMRegressor, early_stopping

    def fp(Xtr, Ytr, Xte, tr, te, q):
        fv = np.unique(grp[tr])
        v = np.isin(grp[tr], np.random.default_rng(0).choice(fv, max(1, round(0.15 * len(fv))), replace=False))
        es = 20 <= v.sum() <= len(v) - 20  # else too few groups for a validation part: the registered 100 rounds
        P = np.empty((len(Xte), Ytr.shape[1]))
        for j in range(Ytr.shape[1]):
            m = LGBMRegressor(
                num_leaves=15,
                learning_rate=0.05,
                n_estimators=1000 if es else 100,
                n_jobs=4,
                random_state=0,
                verbose=-1,
            )
            if es:
                m.fit(
                    Xtr[~v], Ytr[~v, j], eval_set=[(Xtr[v], Ytr[v, j])], callbacks=[early_stopping(30, verbose=False)]
                )
                stops.append(int(m.best_iteration_))
            else:
                m.fit(Xtr, Ytr[:, j])
                stops.append(100)
            P[:, j] = m.predict(Xte)
        return P

    fp.wants_rows = True
    return fp


if sys.argv[1] == "unit":
    import anndata as ad
    from spatialfloors.floors import f3
    from spatialfloors.metrics import f3 as m3
    from spatialfloors.provenance import stamp

    u, OUT = sys.argv[2], Path(sys.argv[3])
    OUT.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(f"runs/p2/cache/{u}.h5ad")
    g = a.var_names.get_indexer(json.loads((PUB / f"manifests/panels/{u}.json").read_text())["genes"])
    assert (g >= 0).all()
    D = f3.prepare(a.X[:, g], a.obs[["x_um", "y_um"]].values, a.obs.fov.values, ids=np.load(f"runs/p2/f3_ids/{u}.npy"))
    L = f3.V3_LOCK["mouse" if u.startswith("allen") else "human"]
    K = f3.kernel(D, L["ell"])
    T = f3.id_mix(D, *K)
    R = {k: f3.ring(D, *v) for k, v in f3.RINGS_V3.items()}
    rep = dict(unit=u, grid={})
    for arm in ("floor_feats", "all_feats"):
        t, per, stops = time.time(), {}, []
        for r in f3.receivers(D):
            rows, Y, dep = f3.target(D, r)
            Ec = f3.svd_fn(D, np.flatnonzero(~D["mask"][r]), L["npc"])
            Es = f3.svd_fn(D, np.flatnonzero(D["mask"][r]), f3.V3_NPC_OWN)
            F = (
                f3.v3_floor_parts(rows, L["feat"], K, T, R, Es, Ec)
                if arm == "floor_feats"
                else f3.v3_parts(D, r, rows, L["feat"], K, T, R, Ec, Es)
            )
            fp = fit_predict_es(f3.groups(D["fov"][rows], D["xy"][rows]), stops)
            s, fo = m3.block_cv(
                F, Y, D["xy"][rows], fp, f3.RF_UM, blocks=D["blocks"][rows], split=D["split"], covars=dep
            )
            per[r] = (s, fo, [np.nan])
        rep["grid"][arm] = dict(
            unit=f3.unit_score(per) if per else float("nan"),
            n_receivers=len(per),
            rounds_median=float(np.median(stops)) if stops else None,
            rounds_at_cap=float(np.mean(np.array(stops) >= 1000)) if stops else None,
            s=round(time.time() - t, 1),
        )
        print(u, arm, rep["grid"][arm], flush=True)
    (OUT / f"f3es_{u}.json").write_text(json.dumps(stamp(rep), indent=1, default=float))
    print("F3ES_UNIT_DONE", u, flush=True)
else:
    import pandas as pd
    from spatialfloors import stats

    OUT, P4 = Path(sys.argv[2]), Path(sys.argv[3])
    Lk = json.loads((PUB / "manifests/locked_v3.json").read_text())["locked"]["F3"]
    donors = {
        f"gse303070_{d['slide'].replace('Tumor', 'Tumour')}_{d['pref']}": d["donor_cells"]
        for d in json.loads((PUB / "manifests/gse303070_donors.json").read_text())["check_vs_fov_csv"]
    }
    rows = []
    for f in sorted(OUT.glob("f3es_*.json")):
        d = json.loads(f.read_text())
        u, sp = d["unit"], "mouse" if d["unit"].startswith("allen") else "human"
        fl = json.loads((P4 / f"f3_{u}_ridge.json").read_text())["grid"]["v3"]["unit"]
        lg = json.loads((P4 / f"f3_{u}_lightgbm.json").read_text())["grid"][Lk["lightgbm"][sp]["config"]]["unit"]
        if not np.isfinite(fl):
            continue
        rows.append(
            dict(
                crop=u,
                floor=fl,
                registered=lg,
                **{k: v["unit"] for k, v in d["grid"].items()},
                **{f"rounds_{k}": v["rounds_median"] for k, v in d["grid"].items()},
            )
        )
    df = pd.DataFrame(rows)
    unit_of = {c: donors.get(c, c) for c in df.crop}
    out, md = (
        dict(n_crops=len(df), pairs={}),
        [
            "# F3 LightGBM with early stopping (exploratory)",
            "",
            "D = LightGBM − v3 floor per unit; registered test (platform strata, margin of F3).",
            "",
            "| LightGBM | Features | Units | Mean score | D | One-sided 95% bounds | Non-inferior | Median rounds |",
            "|---|---|---|---|---|---|---|---|",
        ],
    )
    for key, name, feats in (
        ("registered", "100 rounds (registered)", "all v3 groups"),
        ("all_feats", "early stopping", "all v3 groups"),
        ("floor_feats", "early stopping", "the floor's features"),
    ):
        units, D, strata = stats.aggregate_units(dict(zip(df.crop, df[key] - df.floor)), unit_of)
        t = stats.family_test(D, strata, margin=Lk["margin"])
        out["pairs"][key] = {k: t[k] for k in ("est", "lower", "upper", "noninferior")} | dict(
            n_units=len(D),
            mean_score=float(stats.aggregate_units(dict(zip(df.crop, df[key])), unit_of)[1].mean()),
            per_unit=dict(zip(units, map(float, D))),
        )
        rd = "100" if key == "registered" else f"{df['rounds_' + key].median():.0f}"
        md.append(
            f"| {name} | {feats} | {len(D)} | {out['pairs'][key]['mean_score']:.4f} | {t['est']:+.4f} | "
            f"{t['lower']:+.4f}, {t['upper']:+.4f} | {'yes' if t['noninferior'] else 'no'} | {rd} |"
        )
    fm = float(stats.aggregate_units(dict(zip(df.crop, df.floor)), unit_of)[1].mean())
    md += ["", f"v3 floor mean score {fm:.4f} on the same units."]
    (OUT / "X_f3_lgbm_es.json").write_text(json.dumps(out, indent=1, default=float))
    (OUT / "X_f3_lgbm_es.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
