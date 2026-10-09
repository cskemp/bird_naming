# Functions for the cross-validated version of Analysis 1 (predicting
# name-sharing). See cv_analyses/README.md for the design. These functions
# have no side effects; run_cv.R calls them and saves the results.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(here)
  library(glmmTMB)
})

languages <- c("Anindilyakwa", "Innu", "Saami", "Tlingit",
               "Tobelo", "Tzeltal", "Zapotec")

# The nine models in Table 2 of the paper
models <- list(
  "null"       = character(0),
  "23"         = "Dist23",
  "74"         = "Dist74",
  "Phy"        = "DistPhy",
  "Per"        = "DistPer",
  "74+Phy"     = c("Dist74", "DistPhy"),
  "74+Per"     = c("Dist74", "DistPer"),
  "Phy+Per"    = c("DistPhy", "DistPer"),
  "74+Phy+Per" = c("Dist74", "DistPhy", "DistPer")
)

# Contrasts reported in the paper's discussion: LL(first) - LL(second)
contrasts <- list(
  "74 vs 23"              = c("74", "23"),
  "74+Per vs 74"          = c("74+Per", "74"),
  "74+Phy+Per vs 74+Per"  = c("74+Phy+Per", "74+Per"),
  "74+Phy vs 74"          = c("74+Phy", "74"),
  "Per vs Phy"            = c("Per", "Phy")
)

predictors <- c("Dist74", "Dist23", "DistPhy", "DistPer")

# ---------------------------------------------------------------------------
# Data

# Read the pair data and remove duplicates exactly as Analyses.Rmd does
# (section "Take out duplicates"): strip the .n suffix from overdifferentiated
# species, order each pair, keep SameName = TRUE if any copy of the pair is
# TRUE, and drop self-pairs.
load_pairs <- function(langs = languages) {
  out <- lapply(langs, function(language) {
    d <- read.csv(here("data", paste0("similarity_", language, ".csv")))
    d <- d %>%
      mutate(across(c(Species1, Species2), ~ gsub("\\..*", "", .x))) %>%
      mutate(s1 = pmin(Species1, Species2),
             s2 = pmax(Species1, Species2)) %>%
      select(-Species1, -Species2) %>%
      rename(Species1 = s1, Species2 = s2) %>%
      arrange(desc(SameName)) %>%
      distinct(Species1, Species2, .keep_all = TRUE) %>%
      filter(Species1 != Species2) %>%
      select(Species1, Species2, SameName, all_of(predictors))

    n_sp <- length(unique(c(d$Species1, d$Species2)))
    stopifnot(nrow(d) == choose(n_sp, 2), !anyNA(d))
    d$SameName <- as.numeric(d$SameName)
    d
  })
  names(out) <- langs
  out
}

# z-score predictor columns using the means and SDs of rows `ref`
zscore <- function(d, ref = rep(TRUE, nrow(d))) {
  for (p in predictors) {
    m <- mean(d[[p]][ref])
    s <- sd(d[[p]][ref])
    d[[p]] <- (d[[p]] - m) / s
  }
  d
}

# Randomly assign each species to one of K folds (as equal in size as possible)
make_folds <- function(species, K) {
  f <- sample(rep_len(seq_len(K), length(species)))
  names(f) <- species
  f
}

# ---------------------------------------------------------------------------
# Model fitting and scoring

design <- function(d, preds, intercept_cols = NULL) {
  X <- if (is.null(intercept_cols)) {
    matrix(1, nrow(d), 1, dimnames = list(NULL, "(Intercept)"))
  } else {
    intercept_cols
  }
  if (length(preds) > 0) X <- cbind(X, as.matrix(d[, preds, drop = FALSE]))
  X
}

fit_glm <- function(X, y) {
  fit <- suppressWarnings(glm.fit(X, y, family = binomial()))
  list(coef = fit$coefficients, converged = fit$converged)
}

# Log-probability of each observed outcome under linear predictor eta
# (computed on the log scale so it never underflows to -Inf)
loglik_pairs <- function(eta, y) {
  ifelse(y == 1, plogis(eta, log.p = TRUE), plogis(-eta, log.p = TRUE))
}

# Rank-based AUC (equivalent to the Mann-Whitney U statistic)
auc <- function(p, y) {
  r <- rank(p)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Summarise one rep's held-out scores. ll_sum/p_sum/count are per pair, summed
# over the folds in which that pair was held out. Each pair's score is the mean
# over those folds, so every pair counts once.
score_rep <- function(ll_sum, p_sum, count, y) {
  stopifnot(all(count >= 1))
  ll <- ll_sum / count
  p <- p_sum / count
  c(ll = sum(ll), auc = auc(p, y), brier = mean((p - y)^2))
}

# ---------------------------------------------------------------------------
# Cross-validation holding out species, within one language

# folds: a list of R fold assignments (named integer vectors, one per rep)
cv_language <- function(d, folds, models, cores = 1) {
  sp1 <- d$Species1
  sp2 <- d$Species2
  y <- d$SameName

  one_rep <- function(r) {
    f <- folds[[r]]
    f1 <- f[sp1]
    f2 <- f[sp2]
    K <- max(f)
    acc <- lapply(models, function(m) list(ll = numeric(length(y)),
                                           p = numeric(length(y))))
    count <- numeric(length(y))
    nonconv <- setNames(integer(length(models)), names(models))

    for (k in seq_len(K)) {
      test <- f1 == k | f2 == k
      train <- !test
      if (!any(test)) next
      # no training pair may involve a held-out species
      stopifnot(!any(f1[train] == k | f2[train] == k))
      count[test] <- count[test] + 1

      dz <- zscore(d, ref = train)
      for (m in names(models)) {
        X <- design(dz, models[[m]])
        fit <- fit_glm(X[train, , drop = FALSE], y[train])
        if (!fit$converged) nonconv[m] <- nonconv[m] + 1L
        eta <- drop(X[test, , drop = FALSE] %*% fit$coef)
        acc[[m]]$ll[test] <- acc[[m]]$ll[test] + loglik_pairs(eta, y[test])
        acc[[m]]$p[test] <- acc[[m]]$p[test] + plogis(eta)
      }
    }

    bind_rows(lapply(names(models), function(m) {
      s <- score_rep(acc[[m]]$ll, acc[[m]]$p, count, y)
      tibble(rep = r, model = m, ll = s[["ll"]], auc = s[["auc"]],
             brier = s[["brier"]], n_pairs = length(y),
             nonconverged = nonconv[[m]])
    }))
  }

  bind_rows(parallel::mclapply(seq_along(folds), one_rep, mc.cores = cores))
}

# ---------------------------------------------------------------------------
# Cross-validation holding out species, pooled across languages

# The pooled model has a separate intercept for each language (fixed effects)
# and common slopes. Predictors are z-scored within each language over all its
# pairs, as in the paper; this uses no outcome information. Folds are drawn
# within each language, and fold k holds out fold k of every language.
# folds_list: list over languages of lists of R fold assignments.
cv_pooled_species <- function(pairs, folds_list, models, cores = 1) {
  langs <- names(pairs)
  d <- bind_rows(lapply(langs, function(l) zscore(pairs[[l]]) %>%
                          mutate(Language = l)))
  d$Language <- factor(d$Language, levels = langs)
  y <- d$SameName
  L <- model.matrix(~ Language - 1, d)
  R <- length(folds_list[[1]])

  one_rep <- function(r) {
    f1 <- f2 <- integer(nrow(d))
    for (l in langs) {
      idx <- d$Language == l
      f <- folds_list[[l]][[r]]
      f1[idx] <- f[d$Species1[idx]]
      f2[idx] <- f[d$Species2[idx]]
    }
    K <- max(f1, f2)
    acc <- lapply(models, function(m) list(ll = numeric(length(y)),
                                           p = numeric(length(y))))
    count <- numeric(length(y))
    nonconv <- setNames(integer(length(models)), names(models))

    for (k in seq_len(K)) {
      test <- f1 == k | f2 == k
      train <- !test
      if (!any(test)) next
      count[test] <- count[test] + 1
      for (m in names(models)) {
        X <- design(d, models[[m]], intercept_cols = L)
        fit <- fit_glm(X[train, , drop = FALSE], y[train])
        if (!fit$converged) nonconv[m] <- nonconv[m] + 1L
        eta <- drop(X[test, , drop = FALSE] %*% fit$coef)
        acc[[m]]$ll[test] <- acc[[m]]$ll[test] + loglik_pairs(eta, y[test])
        acc[[m]]$p[test] <- acc[[m]]$p[test] + plogis(eta)
      }
    }

    bind_rows(lapply(names(models), function(m) {
      s <- score_rep(acc[[m]]$ll, acc[[m]]$p, count, y)
      tibble(rep = r, model = m, ll = s[["ll"]], auc = s[["auc"]],
             brier = s[["brier"]], n_pairs = length(y),
             nonconverged = nonconv[[m]])
    }))
  }

  bind_rows(parallel::mclapply(seq_len(R), one_rep, mc.cores = cores))
}

# ---------------------------------------------------------------------------
# Leave one language out

# Fit SameName ~ preds + (1 | Language) to six languages and score the
# seventh in two ways:
#   (a) ll_pop:   population-level intercept (random effect set to zero)
#   (b) ll_recal: slopes from the training languages, with only the intercept
#                 refitted on the held-out language. This removes differences
#                 in base rate (how many pairs share a name), which none of the
#                 measures is meant to predict.
cv_lolo <- function(pairs, models, cores = 1) {
  langs <- names(pairs)
  dz <- lapply(pairs, zscore)
  grid <- expand.grid(held_out = langs, model = names(models),
                      stringsAsFactors = FALSE)

  one_fit <- function(i) {
    h <- grid$held_out[i]
    m <- grid$model[i]
    preds <- models[[m]]
    train <- bind_rows(lapply(setdiff(langs, h), function(l)
      dz[[l]] %>% mutate(Language = l)))
    test <- dz[[h]]
    rhs <- if (length(preds) == 0) "1" else paste(preds, collapse = " + ")
    fit <- glmmTMB(as.formula(paste("SameName ~", rhs, "+ (1 | Language)")),
                   data = train, family = binomial())
    beta <- fixef(fit)$cond
    X <- design(test, preds)
    eta <- drop(X %*% beta[colnames(X)])

    slopes <- if (length(preds) == 0) rep(0, nrow(test)) else
      drop(as.matrix(test[, preds, drop = FALSE]) %*% beta[preds])
    recal <- suppressWarnings(glm(test$SameName ~ 1, offset = slopes,
                                  family = binomial()))
    eta_recal <- slopes + coef(recal)[[1]]

    tibble(held_out = h, model = m, n_pairs = nrow(test),
           ll_pop = sum(loglik_pairs(eta, test$SameName)),
           ll_recal = sum(loglik_pairs(eta_recal, test$SameName)),
           auc = if (length(preds) == 0) 0.5 else auc(eta, test$SameName),
           converged = fit$fit$convergence == 0)
  }

  bind_rows(parallel::mclapply(seq_len(nrow(grid)), one_fit,
                               mc.cores = cores))
}

# ---------------------------------------------------------------------------
# In-sample fits (for comparison with Table 2 of the paper)

in_sample <- function(pairs, models) {
  per_language <- bind_rows(lapply(names(pairs), function(l) {
    d <- zscore(pairs[[l]])
    bind_rows(lapply(names(models), function(m) {
      X <- design(d, models[[m]])
      fit <- fit_glm(X, d$SameName)
      ll <- sum(loglik_pairs(drop(X %*% fit$coef), d$SameName))
      k <- ncol(X)
      tibble(language = l, model = m, ll = ll,
             bic = -2 * ll + k * log(nrow(d)))
    }))
  }))

  # The paper's "All" column: glmmTMB with a random intercept for Language
  all <- bind_rows(lapply(pairs, zscore), .id = "Language")
  pooled <- bind_rows(lapply(names(models), function(m) {
    preds <- models[[m]]
    rhs <- if (length(preds) == 0) "1" else paste(preds, collapse = " + ")
    fit <- glmmTMB(as.formula(paste("SameName ~", rhs, "+ (1 | Language)")),
                   data = all, family = binomial())
    tibble(language = "All", model = m, ll = as.numeric(logLik(fit)),
           bic = BIC(fit))
  }))

  bind_rows(per_language, pooled) %>%
    group_by(language) %>%
    mutate(dbic = bic - min(bic)) %>%
    ungroup()
}

# ---------------------------------------------------------------------------
# Summaries

# Per language x model: mean and SD over reps of held-out LL, and the held-out
# deviance difference from the best model, dDev = 2 * (LL_best - LL_model).
# dDev is on the same scale as the paper's dBIC, without the parameter penalty.
summarise_cv <- function(cv) {
  cv %>%
    group_by(language, rep) %>%
    mutate(ddev = 2 * (max(ll) - ll)) %>%
    group_by(language, model) %>%
    summarise(ll_mean = mean(ll), ll_sd = sd(ll),
              ddev_mean = mean(ddev), ddev_sd = sd(ddev),
              pct_best = 100 * mean(ddev == 0),
              auc_mean = mean(auc), brier_mean = mean(brier),
              n_pairs = first(n_pairs),
              nonconverged = sum(nonconverged), .groups = "drop") %>%
    group_by(language) %>%
    mutate(ddev_of_means = 2 * (max(ll_mean) - ll_mean)) %>%
    ungroup()
}

# Per language x contrast: LL(first) - LL(second) in nats, summarised over
# reps. Positive values favour the first model.
summarise_contrasts <- function(cv) {
  bind_rows(lapply(names(contrasts), function(cn) {
    a <- contrasts[[cn]][1]
    b <- contrasts[[cn]][2]
    cv %>%
      filter(model %in% c(a, b)) %>%
      select(language, rep, model, ll) %>%
      pivot_wider(names_from = model, values_from = ll) %>%
      mutate(diff = .data[[a]] - .data[[b]]) %>%
      group_by(language) %>%
      summarise(contrast = cn, diff_mean = mean(diff), diff_sd = sd(diff),
                pct_favour_first = 100 * mean(diff > 0), .groups = "drop")
  }))
}

# LOLO contrasts: one value per held-out language, plus a sign test and a
# paired t-test across languages. Languages differ a lot in their number of
# pairs, so the t-test uses the difference per 1000 held-out pairs. With 7
# languages these tests have little power (see README).
summarise_lolo_contrasts <- function(lolo, score = c("ll_pop", "ll_recal")) {
  score <- match.arg(score)
  bind_rows(lapply(names(contrasts), function(cn) {
    a <- contrasts[[cn]][1]
    b <- contrasts[[cn]][2]
    w <- lolo %>%
      filter(model %in% c(a, b)) %>%
      select(held_out, model, n_pairs, all_of(score)) %>%
      pivot_wider(names_from = model, values_from = all_of(score)) %>%
      mutate(diff = .data[[a]] - .data[[b]],
             diff_per_1000 = 1000 * diff / n_pairs)
    n_pos <- sum(w$diff > 0)
    n <- sum(w$diff != 0)
    tt <- t.test(w$diff_per_1000)
    w %>%
      select(held_out, diff, diff_per_1000) %>%
      mutate(contrast = cn, score = score, n_favour_first = n_pos,
             n_languages = n,
             sign_test_p = binom.test(n_pos, n)$p.value,
             t_test_p = tt$p.value,
             mean_diff_per_1000 = mean(w$diff_per_1000),
             ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2])
  }))
}
