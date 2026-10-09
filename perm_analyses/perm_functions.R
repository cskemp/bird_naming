# Functions for the permutation tests of Analysis 1 (predicting name-sharing).
# See perm_analyses/README.md for the design. These functions have no side
# effects; run_perm.R calls them and saves the results.
#
# Data loading and the shared helpers (languages, load_pairs, zscore,
# loglik_pairs, predictors) come from cv_analyses/cv_functions.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
})
source(here("cv_analyses", "cv_functions.R"))

# Tests of "does X add to Z?". Z may be empty, in which case the test is plain
# QAP ("does X predict naming at all?").
#   a: overall effect of one measure
#   b: added value of one measure beyond others
#   c: non-nested comparisons, done as encompassing tests in both directions
perm_tests <- list(
  list(test = "74",          x = "Dist74",  z = character(0),         q = "a"),
  list(test = "23",          x = "Dist23",  z = character(0),         q = "a"),
  list(test = "Phy",         x = "DistPhy", z = character(0),         q = "a"),
  list(test = "Per",         x = "DistPer", z = character(0),         q = "a"),
  list(test = "Per | 74",     x = "DistPer", z = "Dist74",             q = "b"),
  list(test = "Phy | 74",     x = "DistPhy", z = "Dist74",             q = "b"),
  list(test = "Phy | 74+Per", x = "DistPhy", z = c("Dist74", "DistPer"), q = "b"),
  list(test = "Per | 74+Phy", x = "DistPer", z = c("Dist74", "DistPhy"), q = "b"),
  list(test = "74 | Phy+Per", x = "Dist74",  z = c("DistPhy", "DistPer"), q = "b"),
  list(test = "74 | 23",      x = "Dist74",  z = "Dist23",             q = "c"),
  list(test = "23 | 74",      x = "Dist23",  z = "Dist74",             q = "c"),
  list(test = "Per | Phy",    x = "DistPer", z = "DistPhy",            q = "c"),
  list(test = "Phy | Per",    x = "DistPhy", z = "DistPer",            q = "c")
)

# ---------------------------------------------------------------------------
# Data as species x species matrices

# Turn one language's pair data into symmetric species x species matrices, one
# per predictor (z-scored over all pairs) plus SameName. Pairs are indexed by
# the upper triangle (row < col); y is SameName in that order and never
# changes. Permuting species relabels the predictor matrices only.
make_language <- function(d) {
  d <- zscore(d)
  species <- sort(unique(c(d$Species1, d$Species2)))
  n <- length(species)
  i <- match(d$Species1, species)
  j <- match(d$Species2, species)
  mats <- lapply(c(predictors, "SameName"), function(p) {
    M <- matrix(0, n, n)
    M[cbind(i, j)] <- d[[p]]
    M[cbind(j, i)] <- d[[p]]
    M
  })
  names(mats) <- c(predictors, "SameName")
  ut <- which(upper.tri(diag(n)), arr.ind = TRUE)
  L <- list(species = species, n = n, row = ut[, 1], col = ut[, 2],
            mats = mats)
  L$y <- pair_vec(L, mats$SameName)
  L
}

# Upper-triangle pair vector of matrix M, optionally after relabelling the
# species by permutation perm (i.e. the pairs of M[perm, perm])
pair_vec <- function(L, M, perm = NULL) {
  if (is.null(perm)) M[cbind(L$row, L$col)] else
    M[cbind(perm[L$row], perm[L$col])]
}

# Fill a symmetric matrix from an upper-triangle pair vector
vec_to_mat <- function(L, v) {
  M <- matrix(0, L$n, L$n)
  M[cbind(L$row, L$col)] <- v
  M[cbind(L$col, L$row)] <- v
  M
}

# Stack the pair vectors of predictor columns `cols` over a list of languages
stack_cols <- function(Ls, cols) {
  if (length(cols) == 0) return(matrix(numeric(0), sum(sapply(Ls, function(L)
    length(L$y))), 0))
  do.call(rbind, lapply(Ls, function(L)
    sapply(cols, function(p) pair_vec(L, L$mats[[p]]))))
}

# Intercept columns: a single intercept within one language, a fixed intercept
# for each language when pooling (as in cv_pooled_species)
intercepts <- function(Ls) {
  n <- sapply(Ls, function(L) length(L$y))
  if (length(Ls) == 1) return(matrix(1, n, 1))
  lang <- factor(rep(names(Ls), n), levels = names(Ls))
  model.matrix(~ lang - 1)
}

# ---------------------------------------------------------------------------
# MRQAP with double-semi-partialling (Dekker, Krackhardt & Snijders 2007)

glm_bin <- function(X, y, start = NULL) {
  suppressWarnings(glm.fit(X, y, family = binomial(), start = start))
}

# Signed root deviance gain: positive when X has the expected (negative)
# slope, i.e. more distant birds are less likely to share a name
signed_root <- function(dev_z, dev_zx, beta) {
  unname(sign(-beta) * sqrt(max(0, dev_z - dev_zx)))
}

# Wald p-value for the last coefficient, treating pairs as independent (the
# test the paper's BIC comparisons implicitly rely on)
naive_p <- function(fit, X) {
  w <- fit$weights
  V <- solve(crossprod(X * sqrt(w)))
  k <- ncol(X)
  2 * pnorm(-abs(fit$coefficients[k] / sqrt(V[k, k])))
}

# Test whether predictor x adds to predictors z, over the languages in Ls
# (one language, or several for the pooled model).
# perms: list over permutations b of a list over languages of species
# permutations. Returns observed statistics and the null distribution.
mrqap <- function(Ls, x, z, perms, cores = 1) {
  y <- unlist(lapply(Ls, `[[`, "y"), use.names = FALSE)
  Z <- cbind(intercepts(Ls), stack_cols(Ls, z))
  xv <- drop(stack_cols(Ls, x))

  fit_z <- glm_bin(Z, y)
  X <- cbind(Z, xv)
  fit_zx <- glm_bin(X, y)
  beta <- fit_zx$coefficients[ncol(X)]
  s_obs <- signed_root(fit_z$deviance, fit_zx$deviance, beta)

  # Residual of x on Z (OLS over pairs), back into per-language matrices
  e <- lm.fit(Z, xv)$residuals
  n_pairs <- sapply(Ls, function(L) length(L$y))
  e_split <- split(e, rep(seq_along(Ls), n_pairs))
  E <- lapply(seq_along(Ls), function(i) vec_to_mat(Ls[[i]], e_split[[i]]))

  start <- c(fit_z$coefficients, 0)
  one_perm <- function(b) {
    ev <- unlist(lapply(seq_along(Ls), function(i)
      pair_vec(Ls[[i]], E[[i]], perms[[b]][[names(Ls)[i]]])),
      use.names = FALSE)
    f <- glm_bin(cbind(Z, ev), y, start = start)
    signed_root(fit_z$deviance, f$deviance, f$coefficients[ncol(Z) + 1])
  }
  null <- unlist(parallel::mclapply(seq_along(perms), one_perm,
                                    mc.cores = cores))
  B <- length(null)

  list(s_obs = s_obs,
       dll = (fit_z$deviance - fit_zx$deviance) / 2,
       beta = unname(beta),
       p_one = (1 + sum(null >= s_obs)) / (B + 1),
       p_two = (1 + sum(abs(null) >= abs(s_obs))) / (B + 1),
       p_naive = unname(naive_p(fit_zx, X)),
       n_pairs = length(y),
       null = null)
}

# Draw B species permutations for every language, up front, so that results
# don't depend on the number of cores. Every test reuses the same draws.
draw_perms <- function(Ls, B) {
  lapply(seq_len(B), function(b) lapply(Ls, function(L) sample.int(L$n)))
}

# Run every test in `tests` for each language and for the pooled model
run_mrqap <- function(Ls, tests, perms, cores = 1, log = message) {
  sets <- c(lapply(names(Ls), function(l) Ls[l]), list(Ls))
  names(sets) <- c(names(Ls), "All")
  out <- list()
  for (s in names(sets)) {
    log("MRQAP: ", s)
    for (t in tests) {
      r <- mrqap(sets[[s]], t$x, t$z, perms, cores)
      out[[length(out) + 1]] <- tibble(
        language = s, test = t$test, question = t$q, s_obs = r$s_obs,
        dll = r$dll, beta = r$beta, p_one = r$p_one, p_two = r$p_two,
        p_naive = r$p_naive, n_pairs = r$n_pairs,
        null_q025 = quantile(r$null, .025), null_q975 = quantile(r$null, .975),
        null = list(r$null))
    }
  }
  bind_rows(out)
}

# ---------------------------------------------------------------------------
# Category bootstrap for LL(74) - LL(23)

# Folk-category membership for one language. Each species is given a single
# "home" category for resampling: the smallest one it belongs to (as in
# Analysis 2). Overdifferentiated species keep their other memberships through
# the original SameName values (see boot_pairs).
home_categories <- function(bird_data, L, language) {
  m <- bird_data[bird_data$Language == language, c("CanonicalName", "IndexFG")]
  m$species <- sub("[.].*", "", m$CanonicalName)
  size <- table(m$IndexFG)
  m$size <- as.numeric(size[m$IndexFG])
  m <- m[order(m$species, m$size, m$IndexFG), ]
  m <- m[!duplicated(m$species), ]
  stopifnot(setequal(m$species, L$species))
  split(match(m$species, L$species), m$IndexFG)
}

# One bootstrap data set for one language. cats: list of species indices per
# home category; draw: indices of the sampled categories (with repeats).
# Nodes are species within sampled category copies. Two nodes in the same copy
# share a name. Nodes in different copies of the same category are not paired
# (they would be the same birds with different names). Other pairs keep their
# original SameName, which is 1 only for overdifferentiated species.
boot_pairs <- function(L, cats, draw, cols) {
  sp <- unlist(cats[draw], use.names = FALSE)
  copy <- rep(seq_along(draw), lengths(cats[draw]))
  cat_id <- rep(draw, lengths(cats[draw]))
  ut <- which(upper.tri(diag(length(sp))), arr.ind = TRUE)
  u <- ut[, 1]
  v <- ut[, 2]
  keep <- copy[u] == copy[v] | cat_id[u] != cat_id[v]
  u <- u[keep]
  v <- v[keep]
  y <- ifelse(copy[u] == copy[v], 1, L$mats$SameName[cbind(sp[u], sp[v])])
  X <- sapply(cols, function(p) L$mats[[p]][cbind(sp[u], sp[v])])
  list(y = y, X = matrix(X, ncol = length(cols), dimnames = list(NULL, cols)))
}

# LL(74) - LL(23) on a list of per-language data sets (from boot_pairs),
# with an intercept per language
dll_74_23 <- function(sets) {
  y <- unlist(lapply(sets, `[[`, "y"), use.names = FALSE)
  n <- sapply(sets, function(s) length(s$y))
  I <- if (length(sets) == 1) matrix(1, n, 1) else
    model.matrix(~ f - 1, data.frame(f = factor(rep(seq_along(sets), n))))
  X <- do.call(rbind, lapply(sets, `[[`, "X"))
  ll <- function(col) {
    D <- cbind(I, X[, col])
    f <- glm_bin(D, y)
    sum(loglik_pairs(drop(D %*% f$coefficients), y))
  }
  c(dll = ll("Dist74") - ll("Dist23"), n = length(y))
}

# Bootstrap within each language and pooled. draws: list over b of a list over
# languages of sampled category indices. The difference is rescaled to the
# original number of pairs so it is on the same scale as the observed one.
run_boot <- function(Ls, cats, draws, cores = 1) {
  cols <- c("Dist74", "Dist23")
  n_orig <- sapply(Ls, function(L) length(L$y))
  one <- function(b) {
    sets <- lapply(names(Ls), function(l)
      boot_pairs(Ls[[l]], cats[[l]], draws[[b]][[l]], cols))
    names(sets) <- names(Ls)
    per <- sapply(names(Ls), function(l) {
      r <- dll_74_23(sets[l])
      r[["dll"]] * n_orig[[l]] / r[["n"]]
    })
    r <- dll_74_23(sets)
    c(per, All = r[["dll"]] * sum(n_orig) / r[["n"]])
  }
  m <- do.call(rbind, parallel::mclapply(seq_along(draws), one,
                                         mc.cores = cores))
  as_tibble(m) %>% mutate(b = row_number(), .before = 1)
}

draw_categories <- function(cats, B) {
  lapply(seq_len(B), function(b) lapply(cats, function(cl)
    sample.int(length(cl), replace = TRUE)))
}

# Observed LL(74) - LL(23), from the same code path (every category drawn once)
observed_74_23 <- function(Ls, cats) {
  cols <- c("Dist74", "Dist23")
  sets <- lapply(names(Ls), function(l)
    boot_pairs(Ls[[l]], cats[[l]], seq_along(cats[[l]]), cols))
  names(sets) <- names(Ls)
  c(sapply(names(Ls), function(l) dll_74_23(sets[l])[["dll"]]),
    All = dll_74_23(sets)[["dll"]])
}

summarise_boot <- function(boot, observed) {
  tidyr::pivot_longer(boot, -b, names_to = "language", values_to = "dll") %>%
    group_by(language) %>%
    summarise(observed = observed[first(language)],
              boot_mean = mean(dll), boot_sd = sd(dll),
              ci_lo = quantile(dll, .025), ci_hi = quantile(dll, .975),
              pct_favour_74 = 100 * mean(dll > 0), .groups = "drop") %>%
    arrange(language == "All", language)
}

# ---------------------------------------------------------------------------
# Calibration: is MRQAP's false-positive rate right?

# For each language, test R "noise" predictors given 74. Each noise predictor
# is the Per matrix with its species labels shuffled, so it has realistic
# matrix structure but no true relation to naming. Each is tested by MRQAP
# (B permutations) and by the naive Wald test.
calibrate <- function(Ls, R, B, cores = 1) {
  bind_rows(lapply(names(Ls), function(l) {
    L <- Ls[[l]]
    noise_perm <- lapply(seq_len(R), function(r) sample.int(L$n))
    perms <- lapply(seq_len(R), function(r)
      lapply(seq_len(B), function(b) setNames(list(sample.int(L$n)), l)))
    one <- function(r) {
      L$mats$Noise <- L$mats$DistPer[noise_perm[[r]], noise_perm[[r]]]
      res <- mrqap(setNames(list(L), l), "Noise", "Dist74", perms[[r]])
      tibble(language = l, rep = r, p_mrqap = res$p_two,
             p_naive = res$p_naive)
    }
    bind_rows(parallel::mclapply(seq_len(R), one, mc.cores = cores))
  }))
}
