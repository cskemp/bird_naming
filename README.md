# bird_naming

This repository contains all of the code and data for:
Wilson, Ferdinand and Kemp, Perceptual Similarity and the Relationship Between Folk and Scientific Bird Classification.
The main notebook for reproducing our results is analyses/Analyses.Rmd

# Folder Structure

## analyses
This folder contains all of the code to run our analyses (Analyses.Rmd) and generate the plots (Plots.Rmd)

## data
This folder contains the preprocessed data used in the analyses.

## output
This folder contains all of the figures and tables generated by our analyses.

## preprocessing
This folder contains a notebook with all of the code used to prepare the raw data for the analyses. It also contains two data folders (data and rawdata)

# Installing R Libraries
From within the R console, run
>renv::install("gitcreds")

then

>renv::restore()

This will install all of the packages used by the code in this repository.
