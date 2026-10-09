# perm_analyses

Permutation tests for Analysis 1 of the paper: does each of the four measures
(Clements 1974 taxonomy `74`, Clements 2023 taxonomy `23`, BirdTree phylogeny
`Phy` and AVONET perceptual distance `Per`) predict whether two birds share a
folk name, does it add to the others, and is 1974 better than 2023?

## Why

Analysis 1 (`analyses/Analyses.Rmd`) compares logistic regressions of about
70k dependent species pairs by BIC (issue A1 in `REVIEW.md`). `cv_analyses/`
gives fairer point estimates by cross-validation, but no p-values. This folder
gives p-values that respect the dependence between pairs. It implements §2 of
`README_analysis_plan.md`, revised as described below. Nothing outside this
folder is changed.

## How to run

From the repository root, using system R (tested with R 4.4.3; needs
tidyverse, here, glmmTMB, knitr, kableExtra and rmarkdown):

```sh
Rscript --vanilla perm_analyses/run_perm.R --quick   # smoke test, ~1 min
Rscript --vanilla perm_analyses/run_perm.R           # full run, ~25 min on 9 cores
Rscript --vanilla -e 'rmarkdown::render("perm_analyses/Perm_Analysis1.Rmd")'
```

`run_perm.R` does all the fitting and caches the results. The notebook only
reads those results (pass `params = list(quick = TRUE)` to render the smoke
test). Data loading and shared helpers come from `cv_analyses/cv_functions.R`.

## Files

| File | Contents |
| --- | --- |
| `perm_functions.R` | Matrix set-up, MRQAP, category bootstrap and calibration functions |
| `run_perm.R` | Driver: checks, then runs everything and writes `output/` |
| `Perm_Analysis1.Rmd` | Report: calibration, tables, figure, bootstrap and limitations |
| `output/perm_tests.rds`, `.csv` | One row per language × test: s, ΔLL, β, permutation and naive p-values (the `.rds` also holds the null distributions) |
| `output/boot_74v23.rds`, `.csv` | Category bootstrap of LL(74) − LL(23) |
| `output/calibration.rds`, `.csv` | False-positive check with noise predictors |
| `output/table_*.csv`, `table_perm.tex`, `perm_tests.png` | Tables and figure written by the notebook |
| `output/run_info.txt`, `run_perm.log` | Settings, session info and log |
| `output/quick/` | Outputs from the `--quick` smoke test |

## Design

**Why not shuffle species among categories for everything?** The original
plan (§2) used that one null for all three questions. It fits only question
(a). A random taxonomy has no structure, so both `74` and `23` beat it easily,
and the spread of LL(74) − LL(23) under that null says nothing about whether the
observed difference is real. The same goes for "does X add to Z?". So:

| Question | Test |
| --- | --- |
| (a) Does X predict naming at all? | QAP: relabel the species of X's distance matrix. This is the same as shuffling species among the folk categories with their sizes fixed. |
| (b) Does X add to Z? | MRQAP with double-semi-partialling (Dekker, Krackhardt & Snijders, 2007) |
| (c) 74 vs 23, Per vs Phy | Encompassing MRQAP in both directions (74 \| 23 and 23 \| 74), plus a category bootstrap of LL(74) − LL(23) |

**Data.** The data are `data/similarity_<Language>.csv`, deduplicated exactly as
in `Analyses.Rmd` (via `load_pairs()`), and turned into species × species
matrices. `SameName` is never changed, so all the dependence between pairs is
kept.

**MRQAP, for X given Z** (Z may be empty), within a language:
1. Fit `SameName ~ Z + X` and `~ Z` with `glm.fit`. The statistic is the signed
   root deviance gain s = sign(−β_X)·√(2ΔLL). s > 0 is the expected direction
   (more distant → less likely to share a name). Because the statistic is
   signed, a wrong-sign `Phy` effect (issue A2) is not counted as support.
2. Regress X on [1, Z] by OLS over the pairs, and keep the residual as a
   matrix.
3. Relabel the residual matrix's species B = 2000 times. Refit
   `SameName ~ Z + residual` each time and record s.
4. p_one = (1 + #{s_b ≥ s})/(B + 1); p_two uses |s|. The naive Wald p-value,
   which treats pairs as independent, is reported for comparison.

**Pooled ("All").** The model has a fixed intercept for each language and
common slopes. The residuals are taken on Z plus the language dummies, and
species are relabelled within each language. The same permutations (seed 302,
drawn up front) are used for every test.

**Tests.**
- (a) `74`, `23`, `Phy`, `Per`
- (b) `Per | 74`, `Phy | 74`, `Phy | 74+Per`, `Per | 74+Phy`, `74 | Phy+Per`
- (c) `74 | 23`, `23 | 74`, `Per | Phy`, `Phy | Per`

If X adds to Y but Y adds nothing to X, X is the better measure. If each adds
to the other, each carries something the other lacks, and the test doesn't
pick a winner.

**Category bootstrap of LL(74) − LL(23).**
- Within each language, the folk categories (`IndexFG` in
  `data/bird_data.csv`) are resampled with replacement, 2000 times.
- Birds in the same copy of a category share a name. Different copies of the
  same category are not paired. All other pairs keep their original
  `SameName`.
- Each overdifferentiated species is resampled with its smallest category
  (9 memberships in total).
- Drawing every category once reproduces the data exactly, and `run_perm.R`
  checks this.
- The differences are rescaled to the original number of pairs.

**Calibration.**
- For each language, 100 noise predictors are made by relabelling the species
  of the `Per` matrix. Each is tested given `74` by MRQAP (B = 199) and by the
  naive Wald test.
- MRQAP should reject about 5% at α = .05.

**Checks in `run_perm.R`.**
- The pair vectors line up with the pair data.
- A relabelled pair vector equals a direct lookup.
- The identity permutation reproduces the observed statistic.
- Observed ΔLLs for nested tests equal the LL differences in
  `cv_analyses/output/in_sample.csv`.
- The identity bootstrap reproduces LL(74) − LL(23).

## Limitations

- **The MRQAP null relabels species.** It keeps the dependence in the outcome
  exactly, but a relabelled residual doesn't keep the clade structure of the
  real one. For example, `74` and `23` differ only in particular clades. The
  calibration uses noise with realistic matrix structure, but not that kind of
  sparse, clade-shaped residual.
- **The tests and the bootstrap answer different questions.** MRQAP asks
  whether a measure's unique information matches the folk names better than
  chance. The bootstrap asks how stable the *size* of the 74-vs-23 advantage is
  when different folk categories are drawn. A measure can pass the first and
  still have a wide bootstrap interval, if its advantage rests on a few large
  categories.
- **The bootstrap is approximate.** Duplicated categories bring duplicated
  between-category pairs.
- **Seven languages.** Pooled tests treat the languages as fixed. For
  carrying over to new languages, see leave-one-language-out in
  `cv_analyses/`.

## Results (full run, B = 2000)

See `Perm_Analysis1.html` for the tables and figure, and its Summary section
for details.

- MRQAP's false-positive rate in the calibration check is 5.0%. The naive Wald
  test's is 23.4%.
- All four measures predict name-sharing in every language (p < .001).
- `Per` adds to `74` in 5 of 7 languages and pooled. It does not in Innu or
  Saami.
- `Phy` never adds to `74` in the expected direction. Its wrong-sign effect
  (A2) is significant in several languages and pooled.
- `74` adds to `23` everywhere. But `23` also adds to `74` in Tzeltal, Zapotec,
  Tobelo and pooled.
- The category-bootstrap 95% interval for LL(74) − LL(23) excludes 0 only in
  Tlingit. Pooled: 90% of resamples favour 74, interval [−53, 318].
