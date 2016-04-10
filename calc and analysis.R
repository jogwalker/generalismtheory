# Generalism Score Calculation and Analysis
# 10 April 2016
###############################

# packages
library(dplyr)
library(tidyr)
library(ggplot2)


###############################
# Functions
###############################


###############################
# Read in clean data
###############################
load("~/dat/generalismtheory/pairs.traits.RData") # associations
load("~/dat/generalismtheory/host.traits.RData") # host traits
# load("~/dat/generalismtheory/host.taxa.RData") # host taxa
load("~/dat/generalismtheory/para.traits.RData") # parasite traits

###############################
# Run Main
###############################