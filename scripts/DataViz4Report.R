#####################################################################################################

#### Forest soils dataviz script                                                  ###################

#### mark.farrell@csiro.au        +61 8 8303 8664         31/05/2021 ################################

#####################################################################################################


#### Set working directory ####
setwd("/Users/markfarrell/OneDrive - CSIRO/Data/ForestSoils")


#### Packages ####
install.packages("ggtern")


library(tidyverse)
library(janitor)
library(PerformanceAnalytics)
library(corrplot)
library(RColorBrewer)
library(plotrix)
library(ggpmisc)
library(ggtern)
library(ggbluebadge)
library(magrittr)


#### data in ####
sum <- read_csv("data/processed/summary.csv")
all <- read_csv("data/processed/ChemAll_adm_OLremPLFA.csv")
