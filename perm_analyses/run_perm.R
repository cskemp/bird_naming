# Permutation tests for Analysis 1 (predicting name-sharing).
#
# Usage (from the repository root):
#   Rscript --vanilla perm_analyses/run_perm.R           # full run
#   Rscript --vanilla perm_analyses/run_perm.R --quick   # smoke test
#
# Results go to perm_analyses/output/ (or perm_analyses/output/quick/).

suppressPackageStartupMessages(library(here))
source(here("perm_analyses", "perm_functions.R"))

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args

B <- if (quick) 99 else 2000        # MRQAP permutations
B_boot <- if (quick) 99 else 2000   # category bootstrap resamples
R_cal <- if (quick) 20 else 100     # noise predictors per language
B_cal <- if (quick) 49 else 199     # permutations per noise predictor
seed <- 302
cores <- max(1, parallel::detectCores() - 1)
RNGkind("L'Ecuyer-CMRG")
out_dir <- if (quick) here("perm_analyses", "output", "quick") else
  here("perm_analyses", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) message(format(Sys.time(), "%H:%M:%S "), ...)
log_msg("B = ", B, ", B_boot = ", B_boot, ", R_cal = ", R_cal,
        ", B_cal = ", B_cal, ", seed = ", seed, ", cores = ", cores)

pairs <- load_pairs(languages)
Ls <- lapply(pairs, make_language)
bird_data <- read.csv(here("data", "bird_data.csv"))
cats <- lapply(languages, function(l) home_categories(bird_data, Ls[[l]], l))
names(cats) <- languages

# Checks ---------------------------------------------------------------------
log_msg("Checks")
set.seed(seed)
for (l in languages) {
  L <- Ls[[l]]
  d <- pairs[[l]]
  # pair vectors line up with the pair data
  i <- match(d$Species1, L$species)
  j <- match(d$Species2, L$species)
  stopifnot(all(L$y[match(paste(pmin(i, j), pmax(i, j)),
                          paste(L$row, L$col))] == d$SameName))
  # a permuted pair vector equals a direct lookup of the relabelled species
  p <- sample.int(L$n)
  M <- L$mats$DistPer
  stopifnot(isTRUE(all.equal(pair_vec(L, M, p), M[p, p][upper.tri(M)])))
  # identity permutation: the residual-based statistic equals the observed one
  id <- list(setNames(list(seq_len(L$n)), l))
  r <- mrqap(Ls[l], "DistPer", "Dist74", id)
  stopifnot(isTRUE(all.equal(r$null, r$s_obs)))
}

# Category bootstrap with every category drawn once reproduces the data
obs_74_23 <- observed_74_23(Ls, cats)

# Observed LL gains agree with the in-sample fits in cv_analyses
insample <- read.csv(here("cv_analyses", "output", "in_sample.csv"))
model_name <- function(cols) {
  short <- c(Dist74 = "74", Dist23 = "23", DistPhy = "Phy", DistPer = "Per")
  ord <- c("Dist23", "Dist74", "DistPhy", "DistPer")
  if (length(cols) == 0) "null" else
    paste(short[intersect(ord, cols)], collapse = "+")
}
ll_insample <- function(l, cols) {
  insample$ll[insample$language == l & insample$model == model_name(cols)]
}
for (l in languages) {
  stopifnot(isTRUE(all.equal(unname(obs_74_23[l]),
                             ll_insample(l, "Dist74") - ll_insample(l, "Dist23"))))
}
log_msg("Checks passed")

# MRQAP ----------------------------------------------------------------------
set.seed(seed)
perms <- draw_perms(Ls, B)
res <- run_mrqap(Ls, perm_tests, perms, cores = cores, log = log_msg)

# Observed gains for nested tests reproduce the in-sample LL differences
for (k in seq_len(nrow(res))) {
  if (res$language[k] == "All") next
  t <- perm_tests[[match(res$test[k], sapply(perm_tests, `[[`, "test"))]]
  full <- ll_insample(res$language[k], c(t$z, t$x))
  base <- ll_insample(res$language[k], t$z)
  if (length(full) == 1 && length(base) == 1)
    stopifnot(isTRUE(all.equal(res$dll[k], full - base, tolerance = 1e-6)))
}

saveRDS(res, file.path(out_dir, "perm_tests.rds"))
write.csv(select(res, -null), file.path(out_dir, "perm_tests.csv"),
          row.names = FALSE)

# Category bootstrap ---------------------------------------------------------
log_msg("Category bootstrap")
set.seed(seed + 1)
draws <- draw_categories(cats, B_boot)
boot <- run_boot(Ls, cats, draws, cores = cores)
saveRDS(list(boot = boot, observed = obs_74_23),
        file.path(out_dir, "boot_74v23.rds"))
write.csv(summarise_boot(boot, obs_74_23),
          file.path(out_dir, "boot_74v23.csv"), row.names = FALSE)

# Calibration ----------------------------------------------------------------
log_msg("Calibration")
set.seed(seed + 2)
cal <- calibrate(Ls, R_cal, B_cal, cores = cores)
saveRDS(cal, file.path(out_dir, "calibration.rds"))
write.csv(cal, file.path(out_dir, "calibration.csv"), row.names = FALSE)
log_msg("Calibration: MRQAP rejects ", round(100 * mean(cal$p_mrqap < .05), 1),
        "%, naive Wald rejects ", round(100 * mean(cal$p_naive < .05), 1),
        "% at alpha = .05")

writeLines(c(paste("B", B), paste("B_boot", B_boot), paste("R_cal", R_cal),
             paste("B_cal", B_cal), paste("seed", seed),
             paste("date", format(Sys.time())),
             capture.output(sessionInfo())),
           file.path(out_dir, "run_info.txt"))
log_msg("Done")
