# brown_britishnames

Data and code for Cecil H. Brown's "British Names for American Birds" (*Journal of
Linguistic Anthropology* 2(1):30–50, 1992), plus perceptual, taxonomic and
phylogenetic similarity measures computed for it using the same methods as the rest
of this repository.

Everything in this folder is self-contained. `BrownSimilarity.Rmd` reads the
reference taxonomies under `preprocessing/data/` but writes only into this folder,
and nothing else in the repository depends on it.

## Source material

| File | What it is |
| --- | --- |
| `brown_britishnamesforamericanbirds.pdf` | The original article |
| `brown_britishnamesforamericanbirds.md` | Transcription of the article |
| `brown_table1.csv` | Brown's Table 1: the 87 British bird names applied to American birds, with his naming-pattern code (A/B/exceptional) |
| `brown_table3.csv` | Brown's Table 3: the similarity judgement experiment |

### brown_table3.csv

Brown showed 34 subjects 26 arrays of bird pictures. Each array holds one British
target bird and six candidate American birds, and subjects picked the American bird
that most resembled the target. One row per target–candidate pair (26 × 6 = 156).

- `array` — array number, 1–26
- `british_target`, `british_family`, `british_binomial` — the target bird
- `american_bird`, `american_family`, `american_binomial` — the candidate
- `n`, `percent` — how many of the 34 subjects chose that candidate
- `binomial_as_printed` — filled in only where the printed table has a typo or an
  abbreviation (e.g. `Accipter striatus`, `Cathatus guttatus`, `Xanthocephalus x.`);
  the `*_binomial` columns hold the corrected, expanded form

Scientific names throughout are as Brown gave them, i.e. Peterson (1980) /
AOU-checklist names of that era, so a good many are now retired
(`Dendroica striata`, `Quiscalus quiscalus`, `Ajaia ajaja`, …).

## BrownSimilarity.Rmd

Computes, for every pair of birds within each array, the similarity measures the
main analysis uses:

- **Perceptual similarity** — Euclidean distance in the 11-dimensional AVONET trait
  space (beak length culmen/nares, beak width, beak depth, tarsus length, wing
  length, Kipp's distance, secondary, hand-wing index, tail length, mass), log
  transformed except for the hand-wing index, then z-scored.
- **Taxonomic similarity** — number of steps up the taxonomic tree before two birds
  share a taxon, for **Clements 1974** and **Clements 2023**: 0 = same species,
  1 = same genus, 2 = same family, 3 = same order, 4 = neither.
- **Phylogenetic similarity** — mean cophenetic distance over the 250 BirdTree
  posterior samples in `nexus_brown.nex`. If that file is removed the notebook
  still runs and reports the other three.

It then ranks the six candidates in each array by each measure and draws two sets of
histograms, one panel per measure, of where the measures ranked the birds the
subjects picked. The dashed line at 3.5 marks the mean rank expected from a measure that
ordered the candidates at random.

- **Top choice only** — one observation per array, the rank of the candidate the
  most subjects chose. 27 observations, because array 14's modal choice is tied
  between Harris' sparrow and chipping sparrow and both are kept.
- **All choices** — every candidate contributes its rank once per subject who chose
  it, so a bird picked by 15 subjects contributes fifteen observations and one
  picked by a single subject contributes one. 883 observations (34 subjects × 26
  arrays, less the response missing from array 25).

Both summary tables carry 95% bootstrap confidence intervals on the mean rank, and a
further section compares every pair of measures on the same replicates (a paired
comparison; p values are uncorrected, and there are six comparisons per block).

**Two resampling schemes are reported side by side**, because they answer different
questions:

- **Over arrays** — draw 26 arrays with replacement. Asks how the measures would
  compare on a *different set of birds*, which is the claim of scientific interest.
- **Over subjects** — hold the arrays fixed and resample each array's 34 choices.
  Asks how they would compare if *different people* judged these same pictures.

The array intervals are several times wider, and the two schemes disagree. Under
array resampling only the gaps involving Clements 1974 are clear; perceptual,
phylogenetic and Clements 2023 cannot be separated from one another, so the
mean-rank ordering of those three should not be read as a ranking. Under subject
resampling perceptual does separate from Clements 2023 — but that only says these
particular arrays would behave the same way with other people, not that the ordering
would survive different birds.

An important caveat is printed alongside those plots. Perceptual distance nearly
always separates all six candidates, but the taxonomic measures average fewer than
two distinct values per array, and in some arrays all six candidates sit at the
same taxonomic distance from the target — Brown drew most candidates from the
target's own family, so taxonomy has little room to discriminate within an array.
An all-six tie gives every candidate a fractional rank of 3.5, which is the spike
in the two taxonomic histograms.

Both replicate `preprocessing/Preprocessing.Rmd`. `TaxSimCalc` is copied verbatim
from there (an .Rmd cannot be `source()`d, and the brief for this folder was to
leave the rest of the repository untouched). `MinVal` and the isolation measures
are not carried over — they are about a bird's nearest neighbour within a whole
language dataset, which has no counterpart in a 7-bird array.

Run it from RStudio, or:

```
Rscript -e 'rmarkdown::render(here::here("brown_britishnames/BrownSimilarity.Rmd"))'
```

It needs `tidyverse`, `here`, `rdist` and `ape`, all already in the project
`renv.lock` (`renv::restore()` if they are not installed).

### Where nexus_brown.nex came from

BirdTree phylogenies come from the tree pruner at <https://birdtree.org>, which
takes a species list and returns a nexus of posterior samples pruned to it. The
trees the main pipeline uses were downloaded that way, one job per language (the job
IDs are recorded in `Preprocessing.Rmd`). Brown's birds needed their own job: the
eight existing per-language trees between them cover only 35 of his 99 birds, and no
single one covers more than 11, so nothing already in the repo could be reused, and
cophenetic distances cannot be combined across separate downloads.

To reproduce it:

1. Knit the notebook once. It writes `brown_birdtree_species.csv`, the 99 birds with
   their BirdLife V3 names — all 99 checked against the BirdTree taxonomy, so the
   pruner finds every one.
2. Paste that file's `ScientificNameV3` column into the pruner and ask for **250
   trees** from **Hackett All Species**. That is the 9993-OTU set listed in
   `preprocessing/data/taxonomy_birdtree.csv`, and the only option that includes
   species placed without genetic data — Hackett Sequenced Species would silently
   drop any of Brown's birds that lack sequence data.
3. Save the download as `brown_britishnames/nexus_brown.nex` and re-knit.

If the file is absent the notebook drops `DistPhy`/`zDistPhy` from the pair files,
`RankPhy` from the table3 file, and the phylogenetic panel from each set of
histograms, rather than failing.

### Scaling decision

The AVONET traits are z-scored across **all 99 unique birds in the table at once**,
not array by array. Brown's 26 arrays are one stimulus set, so this is the analogue
of the pipeline's per-language scaling; scaling inside each 7-bird array would make
the standard deviations unstable and the distances incomparable between arrays.

For the same reason the `z*` columns in the output are z-scores over the 546
within-array pairs, not over all 4851 pairs of the 99 birds — cross-array pairs are
not part of Brown's design.

### Name resolution

Brown's binomials have to be mapped onto current names before they can be joined to
AVONET and the Clements taxonomies. Of the 99 unique birds:

- **84** are still valid Clements 2023 names and are kept as they are.
- **12** are retired 1974 names that `preprocessing/data/eBird_1974_2023.csv` maps
  forward (e.g. `Dendroica striata` → `Setophaga striata`,
  `Nuttallornis borealis` → `Contopus cooperi`).
- **3** resolve through neither file and are hard-coded in the notebook:

  | Brown | Clements 2023 | Why |
  | --- | --- | --- |
  | `Capella gallinago` | `Gallinago gallinago` | In the 1974 list as Common Snipe but with no 2023 partner, because that species was later split; the Old World bird Brown illustrates is *Gallinago gallinago* |
  | `Picoides pubescens` | `Dryobates pubescens` | In neither list; the Downy Woodpecker moved to *Dryobates* |
  | `Quiscalus quiscalus` | `Quiscalus quiscula` | In neither list; older spelling of the Common Grackle's name |

The 1974 name is Brown's own binomial wherever that appears in the 1974 checklist
(54 birds); otherwise the 2023 name is mapped back through the change file, falling
back to the 2023 name where nothing changed. `Family1974` comes from
`genus_match_1974.csv`, with two genera absent from it filled in by hand
(`Hylocichla` → Turdidae, `Limnodromus` → Scolopacidae), and `Order1974` from
`clements_1974_orders.csv`. The 2021 name, used only to key into AVONET, comes from
`Changes_Clements21_Clements23.csv`.

After resolution all 99 birds have complete AVONET measurements, so no bird is
dropped from any array.

Two of Brown's common names do not match the current ones, both harmlessly: his
"sparrow hawk" is *Falco sparverius*, now the American Kestrel (array 4
deliberately contrasts it with the British sparrowhawk, *Accipiter nisus*), and
"woodlark" is Clements' "Wood Lark".

## Generated files

Regenerate all three by running the notebook.

### brown_similarity_pairs.csv

546 rows — all 21 within-array pairs for each of the 26 arrays.

| Column | Meaning |
| --- | --- |
| `array` | Array number |
| `Species1`, `Species2` | Clements 2023 names, sorted alphabetically within the pair |
| `Species1Common`, `Species2Common` | Brown's common names, in the same order |
| `BrownName1`, `BrownName2` | Brown's original binomials, in the order the pair was constructed |
| `pair_type` | `target-candidate` (156 rows) or `candidate-candidate` (390 rows) |
| `DistPer`, `zDistPer` | Perceptual distance, raw and z-scored |
| `DistPhy`, `zDistPhy` | Phylogenetic distance — these two columns exist only when `nexus_brown.nex` is present |
| `Dist23`, `zDist23` | Taxonomic distance under Clements 2023 |
| `Dist74`, `zDist74` | Taxonomic distance under Clements 1974 |

### brown_table3_similarity.csv

156 rows — `brown_table3.csv` with the resolved 2023 names, the six distance
columns and four rank columns appended, so everything sits next to Brown's `n` and
`percent`. This is the file to use for comparing the measures against the subjects'
judgements.

The rank columns order the six candidates *within each array*, 1 being closest:
`RankPer`, `Rank23`, `Rank74` (plus `RankPhy` when the BirdTree nexus is present)
from the distance measures, and `RankHuman` from the vote counts. All are **fractional** (`ties.method = "average"`), which matters
because the taxonomic distances are small integers and tie heavily.

### brown_bird_taxonomy.csv

99 rows — the name resolution table: Brown's binomial and common name, how the 2023
name was arrived at (`NameSource`), the 1974/2021/2023 names with genus, family and
order, and the 11 AVONET traits. Lets the mapping be audited by hand without
re-running the notebook.

### brown_birdtree_species.csv

99 rows — the species list for the BirdTree tree pruner: `ScientificNameV3` with the
Clements 2023 name and Brown's common name alongside for checking.

## Two things to know about the data

Brown's Table 3 accounts for all 34 subjects in 25 of the 26 arrays. **Array 25**
(swift) sums to 33 responses / 97.1% in the published table itself — checked
against the PDF, so it is not a transcription error and nothing here corrects it.

The 1974 and 2023 taxonomic distances genuinely differ, which is why both are kept.
The spoonbill pair in array 24 is 1 step apart in 2023 (*Ajaia ajaja* became
*Platalea ajaja*, sharing a genus with *Platalea leucorodia*) but 2 apart in 1974,
and the woodpecker pair in array 22 is 2 apart in 2023 but 1 apart in 1974 (the
downy woodpecker was *Dendrocopos pubescens*, sharing a genus with *Dendrocopos
major*).
