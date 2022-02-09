############################################################
# Main script. 
# Comment in all scripts you would like to run, 
# and then they will run from this script. 
# To run all scripts relating to the TB application,
# the 'Data.RData' file needs to be loaded. 
############################################################


# Welcome to my specialisation project, a requirement for my study program of
# industrial mathematics at NTNU. 
# A considerable part of the code in this project (most of what is related to MCMC by using NIMBLE), 
# is produced by the authors of the paper "A Hierarchical Framework for Correcting Under-Reporting in Count Data", 
# and found in the supplementary material supporting that paper. 
# I have adapted the code to fit my application and model, but the credit should go to them. 
# This relates both to the simulation study and the application on Tuberculosis data. 



# Needed if running on Markov
Markov = TRUE
if(Markov == TRUE){
  .libPaths("~/Documents/Prosjektoppgave/R_packages")
} else {
  .libPaths("C:/Users/sara_/OneDrive/Documents/R/win-library/4.1")
}

#install.packages("devtools")
#install.packages("ggplot2")
#install.packages("nimble")
#install.packages("coda")
#install.packages("mgcv")
#install.packages("dplyr")
#install.packages("ngspatial")
#install.packages("sp")
#install.packages("spdep")
#install.packages("maps")
#install.packages("mapproj")
#install.packages("matrixStats")
#install.packages("RColorBrewer")
#install.packages("sn")
#install.packages("INLA", repos="https://inla.r-inla-download.org/R/testing")
#install.packages("inlabru")
#install.packages("installr")




library(devtools) 
library(ggplot2) # For reproducing plots seen in the paper.
library(nimble) # For MCMC computation using NIMBLE.
library(coda) # For manipulation of MCMC results.
library(mgcv)
library(dplyr)
library(ngspatial)
library(sp) # For plotting the micro-regions.
library(spdep) # For computing the neighbourhood and adjancency objects.
library(maps) # For adding a map to the plots of Brazil.
library(mapproj)
library(matrixStats)
library(RColorBrewer)
library(sn)
library(INLA)
library(inlabru)


# Model thinning. This model uses a lot of memory (20GB), so if this is not available to you, 
# then consider thinning the model. Suggestion: thin_multiplier <- 10, reduces the memory usage to about 8GB. 
# Consider reducing the number og chains and increase the thinning of the TB model to reduce the memory usage even more. 

thin_multiplier <- 1 



# Setting up the color palette for plots
col.pal <- "PuBuGn"
#display.brewer.pal(name = col.pal, n = 9)
col.pal.violet <- brewer.pal(9, col.pal)[4]
col.pal.blue <- brewer.pal(9, col.pal)[6]
col.pal.seagreen <- brewer.pal(9, col.pal)[7]
col.pal.green <- brewer.pal(9, col.pal)[9]
col.pal.red <- brewer.pal(9, "Spectral")[1]


# Load in some functions. 
# All credit to the authors of "A Hierarchical Framework for Correcting Under-Reporting in Count Data" 
# for these functions, needed to run NIMBLE. 
source("Functions.R")

print("Function.R has been executed")

# Seed for NIMBLE and inlabru, 794637 is used for results in 
# "A Hierarchical Framework for Correcting Under-Reporting in Count Data" 
seed <- 794637 

# Load in the data.
load("Data.RData") 
print("Loading Data.RData has been executed")


#Creating the synthetic data
source("Sim_data.R")
print("Sim_data.R has been executed")

# Running the sensistivity analysis and 
# plotting the coverage using MCMC and inlabru, 
# with correlation 1, 0.6 and 0.4 with the under-reporting covariate w_s
source("Sim_coverage.R")
print("Sim_coverage.R has been executed")

# Running the simulation experiment and plotting the results.
# Using both inlabru and MCMC to run the experiments. 
source("Sim_experiments.R")
print("Sim_experiments has been executed")





