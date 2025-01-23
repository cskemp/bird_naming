# preprocessing/data

This is all the data needed to run preprocessing.Rmd and create the data frames used in all of the analyses and plots in this project.

## anindilyakwa_waddy.csv

This is just the bird names from a book with all of the names of Anindilyakwa plants and animals (Waddy, 1988). I use it here to compare to Holman's data. As it aligns pretty well, I can assume that the rest of his data is also suitable for use in this project. It has the original scientific names that Waddy listed for the birds as well as 2023 Clements scientific names that I added when I preprocessed this for an earlier version of this project. Any spelling mistakes or deviations from Waddy's original text are my errors, as I typed it all out by hand.

## AVONET2_eBird.csv

This is the data from AVONET (Tobias at al. 2022), and contains morphological and behavioural measurements from all bird species in the world. It was downloaded from <https://figshare.com/s/b990722d72a26b5bfead>

## bird_data.csv

A dataset used in an earlier version of this project. I was using this to compare to the new data after it was all processed but it's not used in the current version.

## Changes_Clements21_Clements23.csv, Changes_Clements74_Clements23.csv & Changes_V3_Clements23.csv

These were all downloaded from <https://avibase.bsc-eoc.org/compare.jsp> and contain the differences between historical versions of bird taxonomies. They include instances where a bird concept remains the same but has a name change (either common name or scientific name) and any bird concept that appears in only one of the two compared taxonomies.

## clements_1974_orders.csv

This is a list of all families from the Clements 1974 taxonomy and their corresponding orders.

## clements_2023.csv

This is the full clements 2023 taxonomy, just used to align the families and orders. It's available from <https://www.birds.cornell.edu/clementschecklist/download/>

## colour_data.csv

This data is from <https://figshare.com/s/5602ec236f40db0eb162?file=42156828> and accompanies the paper by Santangeli et al (2023).

## eBird_2022_2023_NT.csv

This is used to align the Clements 2022 names with Clements 2023 names in the Anindilyakwa data. It's downloaded from <https://avibase.bsc-eoc.org/compare.jsp>

## genus_match_1974.csv

A list of all of the genera from the Clements 1974 taxonomy and their corresponding families

## holman_data.csv

The data from Holman (2002) compiled from the files in preprocessing/rawdata

## nexus\_[language].nex

Phylogenies downloaded from <https://birdtree.org/subsets/>. We selected the "Ericson All Species: a set of 10000 trees with 9993 OTUs each" option when we downloaded them. Each has a sample of 250 trees and contains the subset of birds that appear in the language datasets.

## species-list\_[language].csv

These are lists of the birds that should appear in the regions where the language data was collected. Downloaded from <https://avibase.bsc-eoc.org/checklist.jsp?lang=EN>. I selected the region corresponding to where the original language data was collected with as much precision as possible. I did not keep a record of the exact regions I chose when I downloaded these; that's my bad. I use them in the code to identify any bird that appears in the data that does not appear in the region, which could indicate that the bird used to be considered conspecific with another species of bird but subsequently had its name changed, but the original name stayed in the data. I didn't use these lists as gospel, I googled any birds that were flagged in this way.

## taxonomy_birdtree.csv

This is the taxonomy used by Jetz et al. (2012), that mostly follows the Birdlife V3 list but with a few change (details can be found at <https://birdtree.org/taxonomy/>). Downloaded from <https://birdtree.org/subsets/>. I use this to help align the names of birds so I can download the correct subset of phylogenies.

## taxonomy_V3.csv

Version 3 (2018) of the Handbook of Birds of the World and BirdLife Checklist of Birds of the World, downloaded from <https://datazone.birdlife.org/species/taxonomy>. This is the taxonomy that the BirdTree taxonomy is based off, I use this to align the names of the birds in the data with their equivalents in this system. I guess I made this decision before I knew that BirdTree uses an altered version of this that you can also download; I could have just used that instead.

## tlingit_data.csv & zapotec_data.csv

These datasets were created for an earlier version of this project and build on the work of Abbott and Kemp (2020). The code for Abbott and Kemp's original versions of these datasets is available at <https://github.com/joshabbott/birdnaming>

# References

Clements, J. F. (1974). Birds of the World: A Check List.New York, Two Continents Pub. Group.

Clements, J. F., Rasmussen,P. C.,Schulenberg, T. S.,Iliff, M. J., Fredericks, T. A., Gerbracht,J. A., Lepage,D.,Spencer,A.,Billerman,S. M.,Sullivan, B. L., and Wood, C. L. 2023. The eBird/Clements checklist of Birds of the World: v2023. Downloaded from <https://www.birds.cornell.edu/clementschecklist/download/>

Handbook of the Birds of the World and BirdLife International (2018). Handbook of the Birds of the World and BirdLife International digital checklist of the birds of the world. Version 3. Retrieved from [http://datazone.birdlife.org/userfiles/file/Species/Taxonomy/HBWBirdLife_Checklist_v3_Nov18.zip](http://datazone.birdlife.org/userfiles/file/Species/Taxonomy/HBWBirdLife_Checklist_v3_Nov18.zip){.uri}

Holman, E. W. (2002). The Relation Between Folk and Scientific Classification of Plants and Animals. Journal of Classification, 19(1), 131-159. <https://doi.org/10.1007/s00357-001-0036-8>

Jetz, W., G. H. Thomas, J. B. Joy, K. Hartmann, and A. O. Mooers. (2012). The global diversity of birds in space and time. Nature, 491, 444-448.

Santangeli, A., Haukka, A., Morris, W., Arkkila, S., Delhey, K., Kempenaers, B., Valcu, M., Dale, J., Lehikoinen, A., & Mammola, S. (2023). What drives our aesthetic attraction to birds?. Npj Biodiversity, 2(1). <https://doi.org/10.1038/s44185-023-00026-2>

Tobias, J. A., Sheard, C., Pigot, A. L., Devenish, A. J. M., Yang, J., Sayol, F., Neate-Clegg, M. H. C., Alioravainen, N., Weeks, T. L., Barber, R. A., Walkden, P. A., Macgregor, H. E. A., Jones, S. E. I., Vincent, C., Phillips, A. G., Marples, N. M., Montaño-Centellas, F. A., Leandro-Silva, V., Claramunt, S., … Schleuning, M. (2022). AVONET: Morphological, ecological and geographical data for all birds. Ecology Letters, 25(3), 581–597. <https://doi.org/10.1111/ele.13898>

Waddy, J. A. (1988). Classification of plants and animals from a Groote Eylandt Aboriginal point of view. Australian National University.

Abbott, J. T. & Kemp, C. (2020). Birds and Words: Exploring environmental influences on folk categorization. Proceedings of the 42nd Annual Meeting of the Cognitive Science Society.
