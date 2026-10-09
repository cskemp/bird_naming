# Cross-validated model comparison for Analysis 1 (predicting name-sharing).
#
# Usage (from the repository root):
#   Rscript --vanilla cv_analyses/run_cv.R           # full run (R = 50 reps)
#   Rscript --vanilla cv_analyses/run_cv.R --quick   # smoke test (R = 2 reps)
#
# Results go to cv_analyses/output/ (or cv_analyses/output/quick/ with --quick).

suppressPackageStartupMessages(library(here))
source(here("cv_analyses", "cv_functions.R"))

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args

K <- 10
R <- if (quick) 2 else 50
seed <- 302
cores <- max(1, parallel::detectCores() - 1)
out_dir <- if (quick) here("cv_analyses", "output", "quick") else
  here("cv_analyses", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) message(format(Sys.time(), "%H:%M:%S "), ...)
log_msg("K = ", K, ", R = ", R, ", seed = ", seed, ", cores = ", cores)

pairs <- load_pairs(languages)
for (l in languages) log_msg(l, ": ", nrow(pairs[[l]]), " pairs, ",
                             sum(pairs[[l]]$SameName), " sharing a name")

# Fold assignments, drawn up front so the results don't depend on the number
# of cores
set.seed(seed)
folds <- lapply(pairs, function(d) {
  species <- sort(unique(c(d$Species1, d$Species2)))
  lapply(seq_len(R), function(r) make_folds(species, K))
})

# In-sample fits -------------------------------------------------------------
log_msg("In-sample fits")
insample <- in_sample(pairs, models)
saveRDS(insample, file.path(out_dir, "in_sample.rds"))
write.csv(insample, file.path(out_dir, "in_sample.csv"), row.names = FALSE)

# Species held out, within each language -------------------------------------
cv <- bind_rows(lapply(languages, function(l) {
  log_msg("Species-held-out CV: ", l)
  cv_language(pairs[[l]], folds[[l]], models, cores = cores) %>%
    mutate(language = l, .before = 1)
}))

# Species held out, pooled over languages ------------------------------------
log_msg("Species-held-out CV: pooled")
cv_all <- cv_pooled_species(pairs, folds, models, cores = cores) %>%
  mutate(language = "All", .before = 1)
cv <- bind_rows(cv, cv_all)

saveRDS(cv, file.path(out_dir, "cv_species.rds"))
write.csv(summarise_cv(cv), file.path(out_dir, "cv_species_summary.csv"),
          row.names = FALSE)
write.csv(summarise_contrasts(cv),
          file.path(out_dir, "cv_species_contrasts.csv"), row.names = FALSE)

# Leave one language out -----------------------------------------------------
log_msg("Leave-one-language-out")
lolo <- cv_lolo(pairs, models, cores = cores)
saveRDS(lolo, file.path(out_dir, "cv_lolo.rds"))
write.csv(lolo, file.path(out_dir, "cv_lolo.csv"), row.names = FALSE)
write.csv(bind_rows(summarise_lolo_contrasts(lolo, "ll_pop"),
                    summarise_lolo_contrasts(lolo, "ll_recal")),
          file.path(out_dir, "cv_lolo_contrasts.csv"), row.names = FALSE)

writeLines(c(paste("K", K), paste("R", R), paste("seed", seed),
             paste("date", format(Sys.time())),
             capture.output(sessionInfo())),
           file.path(out_dir, "run_info.txt"))
log_msg("Done")
