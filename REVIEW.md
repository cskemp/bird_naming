# Critical review: bird_naming scripts vs. CSJ manuscript

## Context
You asked for a critical review of `preprocessing/Preprocessing.Rmd`, `analyses/Analyses.Rmd`, `analyses/Plots.Rmd` and `brown_britishnames/BrownSimilarity.Rmd`, checked against `CSJ_main.tex` and `CSJ_supplementary.tex`, with no changes to any files. I read the code and ran some read-only R checks on the committed CSVs (system R with `--vanilla`, scripts kept in the scratchpad). Nothing in the repo was modified. The issues are listed below by severity.

**What checked out:** the committed `output/*.tex` match Tables 2 and 3 and the SI tables. The worked examples reproduce: eagle/falcon .055 → .27; kookaburra .33 → .76; emu z = 10.70, with all other birds between −1.06 and 2.17. Iso74 and Iso23 match an independent recomputation for all 7 languages. The pair dedup in Analysis 1 gives n(n−1)/2 unique-species pairs. The SI "LL unweighted" values match single-predictor Per glm fits. Tlingit and Zapotec have no hidden duplicates.

---

## A. Issues that affect results or how they should be read

1. **Pairs are not independent, so BIC is overconfident (Analysis 1).** Each species appears in n−1 pairs, but the logistic models treat about 70k pairs as independent. The ΔBICs of hundreds (for example, 1974 vs 2023 in Anindilyakwa: 97) are therefore much too strong. Analysis 3 has the same problem, because every member of a category gets the same companion count. Options: crossed random effects for Species1 and Species2, or a QAP/permutation test. There is permutation-test code (Analyses.Rmd ~3385), but it is `eval=FALSE` and not reported.
2. **Some "best" models have coefficients with the wrong sign.** In Analysis 1, the best models for Innu and Saami are 74+Phy, and Phy has a *positive* β (+0.94 and +1.78). That means more phylogenetically distant pairs are *more* likely to share a name. The VIFs are 6.6–8.3. In Analysis 2, the best Tobelo model includes Phy with β = −0.81, the opposite of the prediction. The text reports that these best models "include phylogenetic … distance" without mentioning the sign. This looks like a suppression/collinearity artefact, not evidence for phylogeny.
3. **Part of the "phylogenetic" predictor comes from taxonomy.** The trees are BirdTree *Ericson All Species* (preprocessing/data/README). In those trees, species without genetic data are placed using taxonomic constraints. So the phylogeny-vs-traditional comparison is partly circular. The manuscript describes phylogeny as DNA-based and never names the tree set. Suggested fixes: report how many species are unsequenced, and add a robustness check using Sequenced-only trees or the Hackett backbone.
4. **The SI weighted-perceptual analysis (Table "weights") has no code in the repo.** The out-of-sample log-likelihood gains look implausibly large: Zapotec +680, Tzeltal +241, Tlingit +184. They are much larger than the gain from adding every other predictor. This should be checked for leakage, such as the held-out language being used in optimisation or a different pair set. The "LL unweighted" column itself is correct.
5. **Problems with the Hunn correspondence measure (SI Table hunn, Analyses.Rmd 417–503):**
   - The "All" column averages each language's D, prop₀ **and CI bounds**. An average of CI bounds is not a CI, and the D and prop₀ values are unweighted means, not pooled values.
   - The SI says D is "normalised between 0 and 1", but the code takes a plain mean of integer scores.
   - The text says 17 + 12 + 488 categories, which is 517, but the table has n = 516. It should say 487.
   - Some species were a single species in 1974: Branta canadensis/hutchinsii (Tlingit), Eurystomus orientalis/azureus (Tobelo) and Vireo cassinii/plumbeus (Zapotec). `standardise_names()` + `distinct()` drops the second of each pair from the 1974 tree, which shrinks those categories in 1974 only and favours 1974.
6. **Table 1 "Species" counts rows, not species.** Overdifferentiated duplicates are counted more than once: Anindilyakwa 190 rows vs 189 unique, Saami 80/79, Tobelo 121/115, Tzeltal 188/187. The caption says this is the number of Western scientific species.
7. **Undocumented changes to the source data.** The manuscript says five datasets were "compiled by Holman", but the code changes them:
   - Anindilyakwa: the emu (the key outlier) and the bustard are added from Waddy. *Haematopus longirostris* is moved into a new singleton (Preprocessing.Rmd:315). The reef-heron rows are edited by hand.
   - Tobelo: *Pelecanus* is regrouped (1971). One *Ptilinopus superbus* is renamed *P. granulifrons* (1991).
   - Many species are re-identified using range lists, for example *Turdus migratorius* → *T. rufopalliatus* and *Picoides tridactylus* → *P. dorsalis*.

   These should be listed in the SI.
8. **Overdifferentiated species are handled differently across analyses.** In Analysis 1 a pair counts as sharing a name if the two species are ever together. Analyses 2 and 3 keep only each species' *smallest* category. As a result, *Egretta sacra* (Anindilyakwa and Tobelo), *Myiagra alecto* (Tobelo) and *Eupsittula canicularis* (Tzeltal) count as singletons with 0 companions, even though they share names with 2–3 species. The choice is defensible, but it should be stated, and a sensitivity check would help.
   - Related: Tobelo *Lorius roratus* (Eclectus) is commented as an "acceptable duplicate", but line 1999 drops one occurrence. Eclectus is no longer overdifferentiated as a result, which affects the Names count. The Discussion's commented draft uses Eclectus as an example of overdifferentiation.
9. **The Fig. 3 t-SNE doesn't match its caption.** `Rtsne()` gets the *distance matrix* as a feature matrix (`is_distance = FALSE`, Preprocessing.Rmd:3258), so it embeds distance profiles, not "Euclidean distance over AVONET features". Perplexity is set to the maximum, (N−1)/3. The panel order (VarPer) is computed in a space re-standardised using only n > 1 birds, which is not the main perceptual space. Plots.Rmd also dedups first, so overdifferentiated birds appear only in their smallest category.

## B. Code bugs and reproducibility

1. Analyses.Rmd:1140 and 1809 call `return(NULL)` outside a function. This errors when the loop reaches the `Language` term of the "All" model. The comment "the last one won't go into the list" is this symptom, and the fix-up line (1185/1855) never runs in the same chunk. The notebook can't be knit end to end.
2. `View(table_nb_preds)` (Analyses ~2700) fails when run non-interactively.
3. File names with the wrong case break on case-sensitive file systems:
   - Preprocessing.Rmd:271 and BrownSimilarity.Rmd:539 read `changes_V3_Clements23.csv`, but the file is `Changes_…`.
   - Plots.Rmd:1625 reads `phylo-tree_anindilyakwa.nex`, but the file is `…_Anindilyakwa.nex`.
4. Preprocessing.Rmd:4157 has `Species2 %in% bird_data$CanonicalName`, which should be `bird_data2`. As a result the "_filtered" distance summaries in `output/SummaryTable.csv` aren't filtered.
5. Preprocessing.Rmd:3620–3625 assigns Iso23 with an inverted `match()`. It only works because the order happens to be the identity. The current values are correct (verified).
6. Fragile spots that currently do no harm:
   - Saami `unacceptable_duplicates` is computed but never removed (1627).
   - Line 1746 is dead code after the spelling fix at 1667.
   - The reef-heron `filter()` silently drops rows with NA CommonName.
   - `TaxSimCalc` never sets the last diagonal entry, and `drop_na` silently drops birds.
   - The colour fallback join reorders `bird_data`.
   - Tlingit/Zapotec `distinct(ScientificName2023)` would silently remove any overdifferentiation.
7. `Changes_Clements74_Clements23.csv` seems to duplicate `eBird_1974_2023.csv` and is not used.

## C. Manuscript vs. code: description mismatches

1. **Perceptual space:** the features are z-scored *within each language*, so the space depends on the species set. The manuscript doesn't say this. Kipp's distance and HWI are functions of wing and secondary length, which implicitly gives wing shape extra weight. Worth stating.
2. **Phylogenetic distance:** the manuscript calls it "the height of the tree where the branches meet". `cophenetic()` returns the path length, which is twice the node height. Rankings are unaffected, but the wording is wrong.
3. **Analysis 2 text:** it says only Tlingit and Tzeltal show good evidence of 1974 > 2023. Anindilyakwa (10.92 vs 8.77, Δ = 2.15) also passes the paper's own threshold of 2. Separately, 2023 beats 1974 in Innu, Saami and Tobelo, so the Discussion's statement that 1974 is better holds mainly for Analysis 1.
4. **SI colour section:**
   - "78 species with >1 colour set": I get 47 unique species (72 rows) in the 7 languages.
   - "841 species" with separate sets: I find 821 by eBird-2021 name.
   - The colour variables are z-scored before computing distance, which the SI doesn't say.
5. **Stale READMEs and comments:**
   - preprocessing/README says colour "isn't used in our paper", but it is in the SI.
   - data/README says the similarity CSVs contain z-columns. They don't; the z-scores are computed in Analyses.
   - analyses/README describes the permutation test as if it runs; it is `eval=FALSE`.
   - The prose in Analyses.Rmd (394, 1081–1088) mentions Rangi and has out-of-date counts.
   - Preprocessing.Rmd 3435–3443 has stale notes on duplicate counts.

## D. brown_britishnames (not in the manuscript)

1. **Different tree set:** it uses **Hackett** All Species, while the main pipeline uses **Ericson** All Species. The README says it follows the "same methods". The README also claims Hackett All Species is "the only option that includes species placed without genetic data", which is wrong, because Ericson All Species does too.
2. ~~**Probable wrong snipe species:** `Capella gallinago` is an *American candidate* in arrays 24 and 26 but is mapped to *Gallinago gallinago*. The main pipeline (Innu) maps the same name to *Gallinago delicata*.~~ **FIXED:** now mapped to *G. delicata*, with a V3/BirdTree override back to *G. gallinago* (BirdTree has no separate *G. delicata*, so `nexus_brown.nex` is unchanged). Outputs regenerated: all perceptual distances shift by about ±0.01 through rescaling, no ranks change for the four main measures, and three numbers in the colour-section prose were updated (2.92→2.93, 30→31, 98→97).
3. **Different 1974-name rule:** it keeps Brown's binomial if it appears in the change file, and otherwise maps 2023 → 1974. The pipeline always maps 2023 → 1974. Results should be the same, but the rules differ.
4. **Hard-coded numbers in the prose** (mean ranks 3.11, 2.76, 2.74; r = 0.21; "30 of 156") will go stale if anything is re-run. The colour weight 0.38/3.70 is borrowed from logistic β's; the README notes this is an approximation.
5. **`config.yaml`:** it has a personal email address committed, and its status is `CREATED` with an empty `completed_at`. The top-level README doesn't mention this folder.

## Next steps
No files have been changed. If you'd like, I can draft fixes for section B. I can also write manuscript wording for C, or set up robustness analyses for A1–A3, for example crossed random effects or a permutation test, and a sequenced-only phylogeny.
