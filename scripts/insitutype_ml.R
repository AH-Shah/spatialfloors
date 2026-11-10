# InSituType ML for spatialfloors F5 (floors/f5.py insitutype): <dir>/x.mtx (cells x genes counts), profiles.csv
# (genes x labels, linear mean expression, same gene order), neg.txt (per-cell background) -> <dir>/out.csv (clust, prob).
# Call as in the InSituType vignette: insitutypeML(x, neg, reference_profiles).
# v3 (PLAN §19, additive; out.csv unchanged): also <dir>/logliks.csv (the log-likelihoods, for the spatial prior).
suppressPackageStartupMessages({library(InSituType); library(Matrix)})
d <- commandArgs(trailingOnly = TRUE)[1]
if (file.exists(file.path(d, "x_i.bin"))) {  # whole-slide arm (SF_LEAN=1, floors/f5.py): raw CSC arrays of the same matrix
  dm <- scan(file.path(d, "x_dim.txt"), quiet = TRUE)
  x <- new("dgCMatrix", i = readBin(file.path(d, "x_i.bin"), "integer", dm[3], size = 4), p = readBin(file.path(d, "x_p.bin"), "integer", dm[2] + 1, size = 4),
           x = readBin(file.path(d, "x_x.bin"), "double", dm[3]), Dim = as.integer(dm[1:2]))
} else x <- as(readMM(file.path(d, "x.mtx")), "CsparseMatrix")
P <- as.matrix(read.csv(file.path(d, "profiles.csv"), check.names = FALSE))
neg <- scan(file.path(d, "neg.txt"), quiet = TRUE); rownames(P) <- paste0("g", seq_len(nrow(P))); colnames(x) <- rownames(P); rownames(x) <- paste0("c", seq_len(nrow(x)))
r <- insitutypeML(x = x, neg = neg, reference_profiles = P)
write.csv(data.frame(clust = as.character(r$clust), prob = as.numeric(r$prob)), file.path(d, "out.csv"), row.names = FALSE)
write.csv(as.matrix(r$logliks), file.path(d, "logliks.csv"), row.names = FALSE)
