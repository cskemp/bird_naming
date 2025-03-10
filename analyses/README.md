# Analyses

This folder includes two Rmd files - Analyses and Plots. The analyses file has all of the code I used for the analyses that appear in the paper, and the plots file has all of the code for the plots. And also a bit of extra stuff.

# Analyses.Rmd

This Analyses file contains all of the analyses from the paper. These include:

The Hunn Correspondence Measure: This directly compares a folk taxonomy to a Western scientific taxonomy using a method proposed by Eugene Hunn (1977). This compares the dissimilarity of each folk category with regard to the Western scientific taxonomy then returns a weighted dissimilarity measure for the entire taxonomy.

Logistic Regression Analyses predicting whether two birds share a name in a folk taxonomy, with perceptual similarity, phylogenetic similarity and proximity in 1974 and 2023 Clements taxonomy as predictors.

Logistic Regression Analyses predicting whether a bird will be assigned to its own singleton category or be placed into a category with other birds. Perceptual isolation, phylogenetic isolation and isolation in the Clements taxonomy are used as predictors. Throughout this Analyses.Rmd and Plots.Rmd I pretty often refer to singleton categories as monotypic categories; that's what we were originally calling them but we changed the terminology as we were writing the paper. 

Negative Binomial Regression Analyses predicting the number of companions a bird will have in its category (i.e., category size minus one). Same predictors as the previous analysis. This regression is also run without the singleton categories to see if the results hold up without their influence. There is also a bit that takes out the emu outlier from the Anindilyakwa data and re-runs the analysis.

Finally there is a permutation test which keeps the structure of the folk categories, shuffles the birds randomly into the categories and calculates a test score which it compares to that of the real data. It's not absolutely perfect because it keeps the duplicate birds in the data and it's possible for them to end up in a category with themselves. I may change that but at the moment that's what is going on.

The code also generates coefficient tables for each of these analyses and exports them into the output folder as csv files and LaTex-compatible text files. The model selection criterion we have chosen to use is Delta BIC, which is BIC adjusted so that the best model scores 0 (i.e., the best model's score is subtracted from every other model's BIC in the set.) Coefficients in the tables depend on the analysis performed. All tables for regression analyses have the best single-predictor model and the best overall model marked in the LaTex versions.

In addition, the regression analyses generate prediction tables, as the coefficients of logistic and negative binomial regressions are hard to interpret on their own. The models represented in the tables are the best models for each language, including just the ones with the most predictors if more than one is indicated. Each row of these tables show the probability of the outcome variable being TRUE (for logistic regressions) or the predicted count of the outcome variable (for negative binomial regressions) with changing values of the predictor variable, with other included variables held at 0.

## Datasets needed to run this code

Two datasets are loaded in from the data folder in this repository. These datasets were created in the preprocessing/Preprocessing.Rmd file in this repository, please refer to the corresponding readme file for details of how they were created. If you want to use your own data with my code, this is how it must be structured:

### bird_data

This dataframe comes from the data/bird_data.csv in this repository. It has the following columns used in the analyses:

-   Language - the names of the languages where the bird naming data is from.
-   Canonical Name - we have used Clements 2023 scientific names for our canonical names. This column should not contain duplicates (duplicate birds have a suffix)
-   IndexFG - indices or names for each folk category.
-   n - number of birds in the folk category
-   ScientificName2023 & ScientificName1974 - Names in the Clements 2023 and 1974 taxonomies. These are the taxonomies that we used in the Hunn Correspondence measure, but you can switch them for other ones that you want to compare. You'll just have to change the code a little bit to do so.
-   BeakLengthCulmen, BeakLengthNares, BeakWidth, BeakDepth, TarsusLength, WingLength, KippsDistance, Secondary, HandWingIndex, TailLength, Mass - These are all the measurements from AVONET (Tobias et al., 2022) that are used in the preprocessing file to calculate the perceptual distance and isolation measures. They only time they are used in the Analyses file is to run a shapiro test on them (these measures can be downloaded from <https://opentraits.org/datasets/avonet.html>).
-   IsoPer/zIsoPer - unscaled and scaled perceptual isolation scores for each bird (we calculated these by finding the Euclidean distance between each bird and its closest neighbour based on the scaled and log-transformed AVONET variables. See the preprocessing file for details.)
-   VarPer - variance for each folk category of birds based on IsoPer scores. This is only used in the Analyses file to compare the variances for each group across distance measures.
-   IsoPhy/zIsoPhy - unscaled and scaled phylogenetic isolation scores for each bird (we calculated these by taking cophenetic distance between each bird and its closest bird in a phylogenetic tree. See preprocessing file for details.)
-   VarPhy - variance for each folk category of birds based on IsoPhy scores. This is only used in the Analyses file to compare the variances for each group across distance measures.
-   Iso23/zIso23 & Iso74/zIso74 - unscaled and scaled isolation variables for each of the two Western scientific taxonomies (calculated by counting up the steps of the taxonomic tree until the bird and its closest bird share a taxon. See preprocessing file for details)

The following columns aren't used in the Analyses.Rmd code but are in the Plots.Rmd code (described below)

-   CommonName - self-explanatory
-   Genus2023/Family2023/Order2023 & Genus1974/Family1974/Order1974 - genuses, families and orders from Clements taxonomies
-   ScientificNameV3 - names from the BirdLife V3 taxonomy, with some small changes to align with the phylogenetic trees. Details at Jetz et al. (2012).
-   tsne1per/tsne2per - coordinates for dimension reduction. These coordinates were generated in preprocessing/Preprocessing.Rmd

The following columns are in our data but aren't used. (the colour data, if you're interested, comes from Santangeli (2023) - their data can be downloaded from <https://figshare.com/s/5602ec236f40db0eb162?file=42156828>. It's a really cool paper and you should give it a read.)

cBlack, cWhite, cYellow, cBlue, cRed, cGreen, cElaboration, cDiversity, IsoCol, zIsoCol, zcElaboration, zcDiversity

### bird_similarity:

This is a large list indexed by Language, imported from data/similarity\_[language].csv. Each of the list elements is a dataframe that is a pairwise list of all bird pair combinations in the folk taxonomy. These dataframes do include duplicate birds but that is dealt with in the code. Here are the columns used in these analyses:

-   Species1/Species2 - The scientific names of the two birds in the pair. In this data, there should be a line for every distinct pair of birds.
-   SameName - a logical column indicating whether a pair of birds has the same name in the folk taxonomy.
-   DistPer/zDistPer - unscaled and scaled perceptual distances between birds
-   DistPhy/zDistPhy - unscaled and scaled phylogenetic distances between birds
-   Dist23/zDist23 & Dist74/zDist74 - unscaled and scaled distances between birds in Clements taxonomy

(Our data also has columns for colour distances but we don't use them here)

# Plots.Rmd

Here you will find all of the code for generating the plots in our paper, plus probably a bunch of other ones that did not make it into the paper. All of the plots are outputted into the output folder.

Note on colours: Each language is assigned a colour and the panels associated with that language are consistently the same colour in every plot. The only plot where distinguishing the colours is necessary for interpreting the plot is for the plot under the heading "\\## Scatterplot of all languages together". The colours I have chosen are just for aesthetics and I haven't checked if they can be distinguished by individuals with colour vision differences or if they are printed in greyscale as I did not use the scatterplot in my paper. If you require this plot to be colourblind-friendly there is an option to do that (just for the facetted plot with all the languages, called mega_plot) by commenting in and out lines. The instructions are in the code, you just have to comb through.

## Data required

This file requires the same two dataframes described above, as well as the following data:

-   data/Predictions_nb.csv & data/Predictions_tnb.csv - these are predictions from single-predictor negative binomial and zero-truncated negative binomial models. They are generated and saved when you run Analyses.Rmd.

## Plots included

These are the plots included:

-   Correlations between similarity predictors - These are pairs plots for all of the similarity predictors for each language. The bottom triangle has scatterplots with pairs of birds as the datapoints. The diagonal has distributions for each level out the outcome variable. The upper triangle has correlation coefficients.

-   Density plots of all isolation predictors for all languages. These have free scales and the x scale is independent.

-   Density plots for all of the AVONET variables for each language, with both fixed and free scales.

-   Violin plots showing distributions of the distance predictors for pair of birds with the same name and different names. These plots also have data points for each pair of birds behind the violins.

-   Violin plots showing distributions of the isolation predictors for birds in singleton categories and birds who share their categories. Throughout this file singleton categories are often referred to as "monotypic categories".

-   scatter plots of real observations overlaid with predictions for companion count (the number of birds with which a single bird species shares its name). These plots have two regression lines: a dotted line with ribbon confidence intervals shows the predictions from the negative binomial model and 95% confidence intervals and a solid line with error bars shows the predictions from the zero-truncated negative binomial model with 95% confidence intervals. The plots are facetted by language and predictor variable. There are two different versions of this plot: one is tall and skinny, and the other one splits the plot into two and makes it wider.

-   scatter plots with negative binomial regression lines with all of the languages in one plot together. These plots use ggplot's built in argument for calculating the regression line rather than predictions from the models. These plots are pretty but not especially informative, I was just trying something out. There are a lot of weird hacks in there for things like putting the facets in a different order than they are plotted.

-   tsne plots - these are a two dimensional representation of the eleven dimensional perceptual space. There is one plot for each language and each facet of the plot highlights the members of one folk category with all other birds shown as pale grey dots. The final panel in each plot highlights all of the birds in singleton categories. These plots use the tsne coordinates generated in the preprocessing stage and included in bird_data.csv. There is code for an additional tsne plot that puts the names of the Anindilyakwa birds back into the plots.

-   dendrograms for phylogenetic trees. These use nexuses for each set of birds downloaded from <https://birdtree.org/subsets/>. This labels the tips of the tree branches with the 2023 Clements scientific name of each species, but if you wanted to you could change it to common name. There is also code in there to colour the branches of singleton categories, and to highlight any folk category with a most recent common ancestor whose children all belong to that same folk category, and to colour the tips so that birds in the same folk generic category are all coloured the same. And also options to mark nodes that correspond with 2023 Clements families and orders. This stuff is just fun for visualising how the different taxonomies interact with the phylogenetic data; I kind of just got on a roll and tried to see how much cool stuff I could do. 

-   The last few plots are just the mini-plots included in the figure with the birds in the paper - two taxonomic trees and a phylogenetic tree. There should also be a tsne plot but there's not.

References:

Hunn, E. (1977). Tzeltal folk zoology: The classification of discontinuities in nature. Academic Press.

Jetz, W., G. H. Thomas, J. B. Joy, K. Hartmann, and A. O. Mooers. 2012. The global diversity of birds in space and time. Nature 491:444-448. <https://doi.org/10.1038/nature11631>

Santangeli, A., Haukka, A., Morris, W., Arkkila, S., Delhey, K., Kempenaers, B., Valcu, M., Dale, J., Lehikoinen, A., & Mammola, S. (2023). What drives our aesthetic attraction to birds?. Npj Biodiversity, 2(1). <https://doi.org/10.1038/s44185-023-00026-2>

Tobias, J. A., Sheard, C., Pigot, A. L., Devenish, A. J. M., Yang, J., Sayol, F., Neate-Clegg, M. H. C., Alioravainen, N., Weeks, T. L., Barber, R. A., Walkden, P. A., Macgregor, H. E. A., Jones, S. E. I., Vincent, C., Phillips, A. G., Marples, N. M., Montaño-Centellas, F. A., Leandro-Silva, V., Claramunt, S., … Schleuning, M. (2022). AVONET: Morphological, ecological and geographical data for all birds. Ecology Letters, 25(3), 581–597. <https://doi.org/10.1111/ele.13898>
