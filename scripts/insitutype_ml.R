# InSituType ML for spatialfloors F5 (floors/f5.py insitutype): <dir>/x.mtx (cells x genes counts), profiles.csv
# (genes x labels, linear mean expression, same gene order), neg.txt (per-cell background) -> <dir>/out.csv (clust, prob).
# Call as in /Users/crocodilus/code/label_transfer.R (insitutypeML(x, neg, reference_profiles)).
# v3 (PLAN §19, additive; out.csv unchanged): also <dir>/logliks.csv (the log-likelihoods, for the spatial prior).
suppressPackageStartupMessages({library(InSituType); library(Matrix)})
d <- commandArgs(trailingOnly = TRUE)[1]
x <- as(readMM(file.path(d, "x.mtx")), "CsparseMatrix"); P <- as.matrix(read.csv(file.path(d, "profiles.csv"), check.names = FALSE))
neg <- scan(file.path(d, "neg.txt"), quiet = TRUE); rownames(P) <- paste0("g", seq_len(nrow(P))); colnames(x) <- rownames(P); rownames(x) <- paste0("c", seq_len(nrow(x)))
r <- insitutypeML(x = x, neg = neg, reference_profiles = P)
write.csv(data.frame(clust = as.character(r$clust), prob = as.numeric(r$prob)), file.path(d, "out.csv"), row.names = FALSE)
write.csv(as.matrix(r$logliks), file.path(d, "logliks.csv"), row.names = FALSE)
