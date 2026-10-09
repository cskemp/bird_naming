# Analysis plan: testing whether the four measures predict folk bird categories

This plan responds to issue **A1** in `REVIEW.md`: the pairwise models treat
non-independent observations as independent. It sets out how to test whether the
four measures predict which birds share a folk name, and how to compare them. The
measures are Clements 1974 taxonomy (`74`), Clements 2023 taxonomy (`23`),
BirdTree phylogeny (`Phy`) and AVONET perceptual distance (`Per`).

Nothing here has been implemented yet. The plan is to add new notebooks, such as
`analyses/Robustness.Rmd`, and leave `analyses/Analyses.Rmd` unchanged.

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

**Null model: shuffle species among categories.** Within a language, randomly
reassign the species to the folk categories, keeping the number of categories and
their sizes fixed. Each shuffle produces a random taxonomy with the same structure
as the real one, so all dependence between pairs is preserved. The existing
`bird_perm` code (`Analyses.Rmd` ~3385–3540, currently `eval=FALSE`) already does
this shuffle, so it can be reused.

**Statistic.** For each shuffle, refit the pair models and record one of:

- (a) **Overall effect:** the log-likelihood improvement of `SameName ~ X` over
  the intercept-only model.
- (b) **Added value:** the log-likelihood improvement from adding X to a model
  that already includes the other measures, e.g. `74 + Per` vs `74`.
- (c) **1974 vs 2023:** LL(`74`) − LL(`23`). The null here is two-sided: under a
  random taxonomy, neither measure should be favoured.

The p-value is the proportion of shuffles, out of 1000–5000, whose statistic is
at least as large as the observed one. To pool over languages, either sum the
statistic across languages, shuffling within each language, or combine the
per-language p-values (e.g. with Fisher's method).

**Alternative for (b): MRQAP.** Shuffling species among categories tests against
"no structure at all". To test one predictor while controlling for others, it
is more standard to permute the *residuals* of that predictor, as in the
double-semi-partialling method of Dekker, Krackhardt & Snijders (2007). The
`sna` package (`netlogit`) implements this approach. It keeps the correlation
between `74` and `Phy` intact, so it also resolves the wrong-sign `Phy`
coefficients in Innu and Saami (issue A2): if `Phy` adds nothing real, its
contribution won't stand out against this null.

**Notes.**
- With about 70k pairs, refitting thousands of times takes a while but is
  manageable. Use `glm.fit` or `speedglm`, and cache the distance vectors so
  only `SameName` changes between shuffles.
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
  re-run the 1974-vs-2023 and added-value comparisons with the shuffles from
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

1. Section 2(c) and 2(b) for Analysis 1, reusing `bird_perm`. This is the
   quickest way to see whether the main claims survive.
2. Section 1: cross-validation holding out species and holding out languages,
   including the redone weighted-perceptual analysis.
3. Section 3: clustering-recovery figure.
4. Analysis 2 and 3 checks, and the sensitivity analyses.

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
