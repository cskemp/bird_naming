# preprocessing

This has all the code I used to preprocess the datasets used in these analyses. The rawdata and data folders have all the data needed for this. They have their own readme files that explain what all the data is and where it came from.

Note: throughout Preprocessing.Rmd I used the term "Montagnais" to refer to the Innu language, as that is what Holman uses in his dataset. However I subsequently learned that Montagnais is an exonym and the preferred term for the language is Innu, so this term is used throughout the paper and all through this repository outside of this markdown file.

Also, I preprocessed all of the Rangi data and it was only after that when I read the original source for that data that we decided not to include it in the paper.

Here is a walkthrough of what is going on in the Preprocessing.Rmd:

## Functions

There are two custom functions used in this file:

### MinVal

This finds the minimum value of each row of a matrix, ignoring duplicates but allowing legitimate 0 values (like if there were two species of bird that are genuinely identical to one another in all the AVONET measurements).

### TaxSimCalc

Calculates taxonomic similarity. It essentially counts up the steps of a taxonomic tree until two organisms share a taxon, but it does it by comparing values in columns one at a time for each pair of birds and scoring them 0 if it's a match on Species, 1 for Genus, 2 for Family etc. It returns both a matrix and a pairwise list with all of the values.

## Adding all taxonomies to datasets

Each bird in the datasets needs to have names in four different taxonomies:
- 1974 Clements taxonomy (for 1974 traditional similarity)
- 2023 Clements taxonomy (for 2023 traditional similarity)
- BirdLife Version 3 taxonomy (for compatibility with BirdTree for phylogenetic similarity)
- 2021 Clements taxonomy (for compatibility with AVONET)
The first bit of preprocessing is just adding all of the names into the datasets. I think I commented it fairly well so you should be able to figure out what I am doing by going through. During this process there was a huge amount of googling of birds and avibase.com was really helpful here. 
After this we're left with a dataframe called bird_data with all of the birds from each language, their names in each of the taxonomies and a variable called "IndexFG" which identifies members of folk generic groups in all of the languages.

## Calculating all similarity measures 

### Perceptual similarity

I add the AVONET variables into the dataset and calculate the Euclidean distances between each pair of birds in each language dataset. I make a pairwise list for each language, and add these values, both scaled and unscaled. The perceptual isolation measure is the smallest of these distances for each bird, and I identify this and put it into bird_data.
I use tsne to plot the 11-dimensional perceptual space into two dimensions, and add the tsne coordinates to bird_data.
I also calculate variances for each of the folk generic categories. I do this in a kind of weird way, but I couldn't think of a less convoluted way at the time that I did it, and it works fine.


### Phylogenetic similarity

First the nexuses from BirdTree are processed. Each of the 250 samples is turned into a matrix with pairwise cophenetic distances between each pair of birds, and then the matrices are averaged. These are the phylogenetic similarity values which are put into a pairwise list. The smallest value for each bird is added back to bird_data as the phylogenetic isolation measure.
In this stage the phylogenetic data is also processed so it can be visualised later as dendrograms, and outputted into the data folder.
A bit later in this file I calculate variance for each folk generic category based on phylogenetic similarity. I do this in a way less convoluted way than how to did the perceptual variance.

### Taxonomic Similarity

I use the TaxSimCalc function to find the taxonomic similarity for 2023 and 1974 between each pair of birds and put the values into a pairwise list, and the minimum values into bird_data as taxonomic isolation.

### Colour similarity

This isn't used in our paper, but it was just something I tried out as an alternative perceptual similarity measure. It calculates the colour similarity of birds by calculating the euclidean distance based on comparing the percentage of the body covered by a each colour. These colour similarity values were put into a pairwise list, and the smallest values were added to bird_data as a Colour Isolation measure.

## bird_similarity file

I join all of the similarity pairwise lists into one big list for each language, and export them to the data folder.

## Taxonomic trees

I make phylo objects for all of the 1974 and 2023 taxonomic information so that it can be read as a phylogenetic tree or plotted as a dendrogram. These are exported as nexuses into the data folder.

## Exporting bird_data

bird_data is exported into the data folder.

## Summary

I make a summary dataframe that has mean, sd, min and max for all isolation and distance measures. I export this into the output folder.