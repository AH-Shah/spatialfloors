# X_f5_scale: scANVI on whole slides with a crop-size query (exploratory diagnostic)

scANVI (locked nl30_nb) trained on the reference + 50,000 random cells of the slide, then applied to every cell; floor = InSituType + spatial prior on the whole slide. On GSE1 221 a run with crop-equivalent epochs on all cells (69 / 10) scored .588, as the registered run.

| Slide | Cells | Query : reference | scANVI registered | scANVI subsampled | Floor | D registered | D subsampled | D crop |
|---|---|---|---|---|---|---|---|---|
| 110 | 694,553 | 10.6 | 0.659 | 0.691 | 0.686 | -0.027 | +0.005 | +0.020 |
| 120 | 741,155 | 11.3 | 0.734 | 0.765 | 0.790 | -0.056 | -0.024 | -0.025 |
| 210 | 617,498 | 9.5 | 0.593 | 0.691 | 0.702 | -0.109 | -0.010 | -0.015 |
| 221 | 370,074 | 5.7 | 0.588 | 0.688 | 0.716 | -0.128 | -0.028 | +0.000 |
| 222 | 140,542 | 2.2 | 0.628 | 0.632 | 0.615 | +0.013 | +0.017 | +0.022 |
| 242 | 474,471 | 7.3 | 0.627 | 0.681 | 0.703 | -0.076 | -0.022 | -0.012 |

| D (6 slides) | Mean | One-sided 95% bounds | Non-inferior (margin .0255) |
|---|---|---|---|
| D_registered | -0.0639 | -0.1069, -0.0208 | True |
| D_subsampled | -0.0106 | -0.0253, +0.0041 | True |
| D_crop | -0.0015 | -0.0174, +0.0144 | True |
