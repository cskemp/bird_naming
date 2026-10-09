# cv_analyses

Analysis 1 of the paper, redone with cross-validation: which of the four
measures (Clements 1974 taxonomy `74`, Clements 2023 taxonomy `23`, BirdTree
phylogeny `Phy` and AVONET perceptual distance `Per`) best predicts whether two
birds share a folk name?

## Why

Analysis 1 (`analyses/Analyses.Rmd`) fits logistic regressions to about 70k
species pairs and compares models by BIC. The pairs are not independent:

- each species appears in n − 1 pairs;
- if A and B share a name, and B and C do, then A and C must too.

BIC treats every pair as an independent observation, so its differences are
far too confident. See issue A1 in `REVIEW.md`, and §1 of
`README_analysis_plan.md`. This folder replaces in-sample BIC with
log-likelihood on held-out species and held-out languages. Nothing outside
this folder is changed.

## How to run

From the repository root, using system R (tested with R 4.4.3; needs
tidyverse, here, glmmTMB, knitr, kableExtra and rmarkdown):

```sh
Rscript --vanilla cv_analyses/run_cv.R --quick   # smoke test, 2 splits, ~1 min
Rscript --vanilla cv_analyses/run_cv.R           # full run, 50 splits, ~5 min
Rscript --vanilla -e 'rmarkdown::render("cv_analyses/CV_Analysis1.Rmd")'
```

`run_cv.R` does all the model fitting and caches the results. The notebook
only reads those results.

## Files

| File | Contents |
| --- | --- |
| `cv_functions.R` | Data loading, cross-validation, LOLO and summary functions |
| `run_cv.R` | Driver: runs everything and writes `output/` |
| `CV_Analysis1.Rmd` | Report: tables, figure and limitations |
| `output/in_sample.csv` | Full-data fits: LL, BIC and ΔBIC (reproduces Table 2) |
| `output/cv_species.rds` | Species-held-out results, one row per language × split × model |
| `output/cv_species_summary.csv` | The same results, summarised over splits |
| `output/cv_species_contrasts.csv` | Key model contrasts, summarised over splits |
| `output/cv_lolo.rds`, `.csv` | Leave-one-language-out LL per held-out language × model |
| `output/cv_lolo_contrasts.csv` | LOLO contrasts, with sign tests and t-tests across languages |
| `output/table_*.csv`, `table_ddev.tex`, `cv_contrasts.png` | Tables and figure written by the notebook |
| `output/run_info.txt`, `run_cv.log` | Settings, session info and log |
| `output/quick/` | Outputs from the `--quick` smoke test |

## Design

**Data.** The data are `data/similarity_<Language>.csv` for the 7 languages in
the paper. Duplicate pairs are removed exactly as in `Analyses.Rmd`: the `.n`
suffix of overdifferentiated species is stripped, and a pair counts as sharing
a name if any copy of it does. `in_sample.csv` reproduces the ΔBICs of Table 2
of the paper exactly, which confirms the data are the same.

**Models.** The nine models of Table 2: `null`, `23`, `74`, `Phy`, `Per`,
`74+Phy`, `74+Per`, `Phy+Per`, `74+Phy+Per`.

**Species held out (within each language).**
- Species are randomly split into K = 10 folds.
- For fold k, each model is fitted to the pairs with *neither* species in fold
  k, and scored on the pairs with *at least one* species in fold k. Splitting
  pairs at random would not work, because the training set would almost always
  contain other pairs involving the same species.
- A pair whose two species are in different folds is held out twice. Its score
  is the mean over those two folds, so every pair counts once and the summed
  held-out LL is on the same scale as the in-sample LL.
- This is repeated with 50 random splits. The fold assignments are drawn up
  front from seed 302, so the results don't depend on the number of cores.
- Reported: held-out LL, the deviance difference from the best model
  ΔDev = 2 × (LL_best − LL_model), which is on the same scale as ΔBIC, AUC, and
  the Brier score.

**Species held out (pooled, "All").** Folds are drawn within each language, and
fold k holds out fold k of every language. The model has a fixed intercept for
each language and common slopes. The paper's `(1 | Language)` random intercept
is replaced by fixed intercepts for speed; with every language in the training
set the two are nearly equivalent.

**z-scoring.** Within one language, z-scoring doesn't change an unpenalised
glm's predictions. Predictors are still z-scored using the training pairs, for
readable coefficients. In pooled models, each language is z-scored over all its
own pairs, as in the paper. This uses predictor values only, never outcomes, so
nothing leaks from the held-out data.

**Leave one language out (LOLO).** `SameName ~ predictors + (1 | Language)` is
fitted with glmmTMB to six languages and scored on the seventh, in two ways:
- **population intercept** (random effect set to zero). This is the strictest
  test, but it also depends on whether the base rate of shared names
  transfers;
- **recalibrated intercept**: the slopes are kept and only the intercept is
  refitted on the held-out language. This asks whether the *relationship*
  between distance and name-sharing transfers. The same single parameter is
  refitted for every model.

Differences between models are given per 1000 held-out pairs, and the
cross-language tests use these normalised values.

**Contrasts.** `74` vs `23`, `74+Per` vs `74`, `74+Phy+Per` vs `74+Per`,
`74+Phy` vs `74`, and `Per` vs `Phy`.

## Limitations

- **Leave-one-language-out has only 7 languages, so its inferential power is
  limited.**
  - A two-sided sign test can't reach p < .05 unless all 7 languages agree,
    and even then p ≈ .016. If 6 of 7 agree, p ≈ .125.
  - A paired t-test across languages has 6 df and low power.
  - The 7 languages are a convenience sample, not a random sample of the
    world's folk taxonomies. They come from different ethnographers and are
    clustered by region. So "transfers to a new language" means "to languages
    like these".

  LOLO results are therefore best read as effect sizes, language by language,
  rather than as a test.
- **The species-held-out CV isn't a significance test either.** The SD across
  splits shows how much the result depends on the random fold assignment, not
  sampling uncertainty.
- **What the CV does provide** is a fairer point estimate of predictive
  performance than in-sample BIC. Valid p-values for "does this measure add
  anything?" and "is 1974 better than 2023?" come from the permutation tests,
  MRQAP and category bootstrap in `perm_analyses/`.
