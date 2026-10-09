# Analysis plan: testing whether the four measures predict folk bird categories

This plan responds to issue **A1** in `REVIEW.md`: the pairwise models treat
non-independent observations as independent. It sets out how to test whether the
four measures predict which birds share a folk name, and how to compare them. The
measures are Clements 1974 taxonomy (`74`), Clements 2023 taxonomy (`23`),
BirdTree phylogeny (`Phy`) and AVONET perceptual distance (`Per`).

Status: §1 is implemented in `cv_analyses/`, and §2 in `perm_analyses/`, with
the changes to §2 described there. `analyses/Analyses.Rmd` is unchanged.

## The problem

Each language's folk taxonomy partitions its species into named categories.
Analysis 1 turns that partition into about 70k species pairs. Each pair is coded
`SameName` and regressed on the distances, and models are compared by BIC. The
pairs are not independent, for two reasons:

1. **Shared species.** A species appears in n − 1 pairs. A 19-member category
   contributes 171 "same" pairs, so a handful of large categories dominate the fit.
2. **Transitivity.** If A and B share a name, and B and C share a name, then A and
   C must share one too. The evidence comes from the categories, not from the
   pairs inside them.

BIC treats every pair as an independent observation, so its differences are far
too large. For example, 1974 beats 2023 by ΔBIC = 97 in Anindilyakwa, yet Hunn's
category-level count in the SI is only 17 categories to 12. Analysis 3 has the
same problem, because every member of a category gets the same companion count.
Analysis 2 (one row per species) is much less affected.

Random effects for species, e.g. a multiple-membership term
`(1 | mm(Species1, Species2))` in `brms`, deal with the first problem but not the
second. They are worth reporting as a sensitivity check, but inference should not
rest on them.

## Recommended approach

The primary analysis has three parts:

| Question | Method | Section |
| --- | --- | --- |
| Which measure, or combination of measures, predicts folk categories best? | Cross-validation that holds out species, plus leave-one-language-out | 1 |
| Does a measure add anything beyond the others? Is 1974 really better than 2023? | Permutation tests that keep the category structure | 2 |
| An intuitive check that needs no pair model | Recovering the folk categories by clustering | 3 |

### 1. Cross-validation that holds out species

**Aim.** Compare models (`74`, `23`, `Phy`, `Per` and their combinations) on how
well they predict pairs they were not fitted to, in place of in-sample BIC.

**Procedure, within each language.**

1. Randomly split the species into K = 10 folds. Overdifferentiated duplicates
   (`.n` suffixes) go in the same fold as their base species.
2. For each fold k:
   - Fit each candidate model to pairs where **both** species are outside fold k.
   - Compute the log-likelihood on pairs where **at least one** species is in
     fold k.
3. Sum the held-out log-likelihood across folds. Repeat with 20–50 different
   random fold splits, and report the mean and the spread across splits.
4. Report the difference in held-out log-likelihood between models (e.g. `74`
   vs `23`, `74` vs `74+Per`). AUC or Brier score can be reported alongside it.

**Across languages.** For the "All" model, also hold out one whole language at a
time, fitting on the other six. This tests directly whether a measure's
relationship to folk naming carries over to other cultures, which is the paper's
central claim. It also gives seven independent estimates of each model's
advantage, so a simple sign test or paired t-test across languages is valid. The
SI weighted-perceptual analysis (issue A4) should be redone this way, with any
weight optimisation done only on the training languages.

**Notes.**
- Hold out species, not pairs. If pairs are split at random, the training set
  almost always contains other pairs involving the same species, so the
  held-out score is inflated.
- Z-score the predictors using the training data only.
- Use plain `glm` (binomial) within a language. Keep the `(1 | Language)` random
  intercept only for the pooled model.

### 2. Permutation tests that keep the category structure

**Aim.** Get valid p-values for (a) whether a measure predicts naming at all,
(b) whether it adds anything beyond the other measures, and (c) whether 1974 fits
better than 2023.

**Which null for which question.** Shuffling species among the folk
categories (keeping the number of categories and their sizes fixed) preserves
all the dependence between pairs. But it only tests against "no folk structure
at all". That suits (a). It does not suit (b) or (c): a random taxonomy has no
structure, so every measure beats it easily, and the spread of LL(`74`) −
LL(`23`) under that null says nothing about whether the observed difference is
real. As implemented in `perm_analyses/`:

- (a) **Overall effect: QAP.** Relabel the species of X's distance matrix.
  This is equivalent to the species shuffle. The statistic is the signed root
  log-likelihood gain of `SameName ~ X` over the intercept-only model.
- (b) **Added value: MRQAP** with double-semi-partialling (Dekker, Krackhardt &
  Snijders, 2007). Regress X on the other measures Z, relabel the species of
  the residual, and refit `SameName ~ Z + residual`. This keeps the
  correlation between e.g. `74` and `Phy` intact, so it also tests the
  wrong-sign `Phy` coefficients in Innu and Saami (issue A2). The statistic is
  signed, so a wrong-sign effect does not count as support.
- (c) **1974 vs 2023: encompassing MRQAP plus a category bootstrap.** Test
  `74 | 23` and `23 | 74`. If only the first is significant, 1974 is better.
  For an effect size, resample the folk categories with replacement and give a
  percentile interval for LL(`74`) − LL(`23`).

Pooling over languages uses a model with a fixed intercept for each language
and common slopes, with species relabelled within each language. The existing
`bird_perm` code (`Analyses.Rmd` ~3385–3540) does the species shuffle, but with
a within-category variance statistic, so it was not reused.

**Notes.**
- With about 70k pairs, refitting thousands of times takes a while but is
  manageable with `glm.fit`. `SameName` stays fixed and only the species
  labels of the predictor (or residual) matrix change between permutations.
- Fix the random seed and save the null distributions to `output/`.

### 3. Recovering the folk categories by clustering

**Aim.** Ask the question in Hunn's own terms: if you cluster the birds using
measure X, how closely do the clusters match the folk categories? This needs no
pair-level model, and it is easy for readers to follow.

**Procedure, within each language.**

1. For `Per` and `Phy`, build an average-linkage hierarchical clustering (UPGMA)
   from the distance matrix. For `74` and `23`, use the taxonomic tree itself,
   cutting at species, genus, family or order level, or use UPGMA on the step
   distances with ties broken at random.
2. Cut each tree so it gives the same number of clusters as the language has folk
   categories. Also report results over a range of cut heights.
3. Score how well the clusters match the folk partition, using the adjusted Rand
   index (Hubert & Arabie, 1985) or variation of information (Meilă, 2007).
4. For combinations of measures, cluster on a weighted sum of the z-scored
   distances, with the weights fit by leave-one-language-out (section 1).
5. Use the shuffles from section 2 to get a null distribution for each score.

**Output.** A table or dot plot of ARI by measure and language, with null 95%
bands. This would make a natural SI figure, and possibly a main-text one.

## The other analyses

- **Analysis 2 (singletons).** Each row is a species, so the observations are
  close to independent and the current models are reasonable. As a cheap check,
  re-run the 1974-vs-2023 and added-value comparisons with the permutation tests from
  section 2 and with species-level cross-validation.
- **Analysis 3 (companion counts).** The outcome is shared by all members of a
  category. Either:
  - model one row per category, regressing category size on its members' mean
    (or minimum) isolation, or
  - keep the species-level model but take p-values from the shuffles in
    section 2.
- **Overdifferentiated species (issue A8).** Run each new analysis twice: once
  with each species kept in its smallest category (as now), and once with each
  species keeping all of its memberships.

## Sensitivity analyses to report alongside

1. A pair model with species random effects, e.g.
   `brms::brm(SameName ~ z74 + zPer + (1 | mm(Species1, Species2)), family = bernoulli())`,
   fit within each language. Report its coefficients and intervals, while noting
   that it does not account for transitivity.
2. The phylogenetic comparisons repeated with BirdTree **Sequenced Species**
   trees, or the Hackett backbone (issue A3). This checks that any phylogeny
   result doesn't come from species placed using taxonomic constraints.
3. Commonality or dominance analysis of the held-out log-likelihood. This
   splits the predictive gain into what is unique to `74`, unique to `Phy`, and
   shared between them, which makes the collinearity between them (issue A2)
   explicit instead of letting it show up as coefficient signs.

## Order of work

1. Section 2 for Analysis 1 (done: `perm_analyses/`).
2. Section 1: cross-validation holding out species and holding out languages
   (done: `cv_analyses/`). The redone weighted-perceptual analysis is still to
   do.
3. **Next:** bootstrap intervals on coefficients and contrasts (see "Next
   step" below).
4. Section 3: clustering-recovery figure.
5. Analysis 2 and 3 checks, and the sensitivity analyses.

## Next step: intervals from the category bootstrap

**Aim.** The permutation tests in `perm_analyses/` give valid p-values, but no
intervals. Report 95% intervals on the coefficients of the Analysis 1 models,
and on their log-likelihood contrasts, that reflect the dependence between
pairs. A simple Bayesian fit or the glm's Wald intervals won't do. Both rest on
the pair likelihood, which treats ~70k pairs as independent. In the
calibration check, naive tests at α = .05 found effects in 23% of noise
predictors.

**Procedure.** Extend the category bootstrap in `perm_analyses/perm_functions.R`
(`home_categories()`, `draw_categories()`, `boot_pairs()`), which currently
covers only LL(74) − LL(23):

1. Within each language, resample the folk categories with replacement,
   B = 2000, as now. For the pooled model, resample within each language and
   fit with a fixed intercept per language and common slopes.
2. For each resample, fit the nine models of Table 2 and record:
   - every slope coefficient;
   - LL differences for the contrasts in `cv_analyses/cv_functions.R`
     (`74 vs 23`, `74+Per vs 74`, `74+Phy+Per vs 74+Per`, `74+Phy vs 74`,
     `Per vs Phy`), rescaled to the original number of pairs.
3. Keep the predictors on their full-data z-scale (as in `make_language()`),
   so that coefficients are comparable across resamples.
4. Report the estimate, the bootstrap SD, the 95% percentile interval and the
   % of resamples with the expected sign. Also report the ratio of the
   bootstrap SD to the naive Wald SE, which measures how overconfident the
   pair model is.

**Checks.**
- Drawing every category once reproduces the full-data glm coefficients.
- **Coverage:** reuse the noise predictors from the calibration check (`Per`
  with its species labels shuffled, entered alongside `74`). The 95% bootstrap
  interval for the noise coefficient should include 0 in about 95% of cases.
  There are only 45–106 categories per language, so percentile intervals may
  be too narrow. If coverage is clearly below 95%, try BCa intervals, or
  report the observed coverage next to the intervals.

**Output.** A table of coefficients with intervals, by model and language,
next to the paper's Table 2. Also a figure of the contrasts with their
intervals, next to the permutation p-values. Add these to
`perm_analyses/Perm_Analysis1.Rmd` (or to a new `Boot_Analysis1.Rmd`).

**Caveats.** The bootstrap is approximate: duplicated categories bring
duplicated between-category pairs, and a few large categories dominate. Treat
the intervals as effect sizes alongside the permutation p-values, not as a
replacement for them.

## References

- Dekker, D., Krackhardt, D., & Snijders, T. A. B. (2007). Sensitivity of MRQAP
  tests to collinearity and autocorrelation conditions. *Psychometrika*, 72(4),
  563–581.
- Hubert, L., & Arabie, P. (1985). Comparing partitions. *Journal of
  Classification*, 2, 193–218.
- Meilă, M. (2007). Comparing clusterings — an information based distance.
  *Journal of Multivariate Analysis*, 98(5), 873–895.
- Hunn, E. (1976). Toward a perceptual model of folk biological classification.
  *American Ethnologist*, 3(3), 508–524.
