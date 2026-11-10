# X_posctl: positive controls of the learned models (run before confirmation; rules fixed before each run)

| Family | Model | Control | Pass rule | Observed | Pass |
|---|---|---|---|---|---|
| F1 | scVI | planted types (5 types, 200 genes) | F1 ≥ floor − .05 and ≥ chance + .20 | scVI 0.415, floor 0.416, chance -0.021 | yes |
| F1 | Nicheformer | vocabulary gate (cannot read synthetic genes) | ≥ 80% panel in vocabulary; F1 ≥ chance + .05 on every dev crop | see runs/p3/f1_* | see note |
| F2 | CellCharter | 2 domains differing only in neighbourhood | ARI ≥ .80 (k = 2, seeds 0–2) | ARI 0.967 | yes |
| F2 | Novae | authors' tutorial slide (Xenium colon P5) | ARI vs published domains ≥ .50; eff. rank ≥ 5; largest domain ≤ 50% | ARI 0.540, rank 7.2, largest 0.36 | yes |
| F3 | LightGBM | planted niche toys (quadratic, couple, bleed) | quadratic ≥ ridge + .08; couple ≥ .15; bleed null | {'quadratic': True, 'couple': True, 'bleed': True, 'pass': True} | yes |
| F3 | NCEM authors' setup edA | planted toys (compo, couple, bleed) | compo / couple recovered; bleed null | {'bleed': True, 'compo': True, 'pass': True} | yes |
| F3 | NCEM authors' setup linA (locked: linA) | planted toys (compo, couple, bleed) | compo / couple recovered; bleed null | {'bleed': True, 'compo': False, 'pass': False} | no |
| F3 | NCEM authors' setup edT_type | planted toys (compo, couple, bleed) | compo / couple recovered; bleed null | {'bleed': True, 'compo': False, 'pass': False} | no |
| F3 | NCEM authors' setup edT_ctx | planted toys (compo, couple, bleed) | compo / couple recovered; bleed null | {'bleed': True, 'couple': False, 'pass': False} | no |
| F3 | NCEM ncem_type (earlier port, superseded) | planted niche toys | as LightGBM | {'quadratic': False, 'couple': False, 'bleed': True, 'pass': False} | no |
| F3 | NCEM ncem_type+ctx (earlier port, superseded) | planted niche toys | as LightGBM | {'quadratic': False, 'couple': False, 'bleed': True, 'pass': False} | no |
| F5 | scANVI | planted-type toy (query depth ×0.3, platform factor) | coarse macro-F1 ≥ .90 | 1.000 | yes |
| F5 | Tangram | same toy | coarse macro-F1 ≥ .90 | 1.000 | yes |

F4 (Cellina, CPA): the control is the reproduction of the authors' published GSE1 results (runs/p0_cellina; PROTOCOL §3).
