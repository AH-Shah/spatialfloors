# X_f1_euclid: F1 with Euclidean k-NN (exploratory, after results-v1; no decision weight)

commit 0714d4f3e838ef7ef7cc2cecf4f14294e1fb8965; margin 0.0194878990028946.

Crop means: cosine_floor 0.1624; cosine_scvi 0.1536; cosine_nicheformer 0.1245; euclid_floor 0.1489; euclid_scvi 0.1113; euclid_nicheformer 0.1235

| Comparison | Units | D | One-sided 95% bounds | Non-inferior | Floor superior |
|---|---|---|---|---|---|
| scvi euclid − floor euclid | 25 | -0.0376 | -0.0429, -0.0323 | True | True |
| scvi euclid − floor cosine | 25 | -0.0510 | -0.0566, -0.0454 | True | True |
| scvi cosine − floor cosine (registered) | 25 | -0.0087 | -0.0112, -0.0062 | True | True |
| scvi best − floor cosine | 25 | -0.0087 | -0.0112, -0.0062 | True | True |
| nicheformer euclid − floor euclid | 25 | -0.0256 | -0.0310, -0.0202 | True | True |
| nicheformer euclid − floor cosine | 25 | -0.0389 | -0.0443, -0.0336 | True | True |
| nicheformer cosine − floor cosine (registered) | 25 | -0.0379 | -0.0432, -0.0327 | True | True |
| nicheformer best − floor cosine | 25 | -0.0378 | -0.0430, -0.0325 | True | True |
