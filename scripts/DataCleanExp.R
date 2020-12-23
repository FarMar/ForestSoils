#####################################################################################################

#### Forest soils data cleaning and analysis script                               ###################

#### mark.farrell@csiro.au        +61 8 8303 8664         17/12/2020 ################################

#####################################################################################################


#### Set working directory ####
setwd("/Users/markfarrell/OneDrive - CSIRO/Data/ForestSoils")


#### Packages ####
install.packages("PerformanceAnalytics")
install.packages("corrplot")


library(tidyverse)
library(janitor)
library(PerformanceAnalytics)
library(corrplot)


#### Import data ####
dat <- read_csv("data/working/MasterFieldDataFC_NSW - Data.csv")
spec(dat)

mir <- read_csv("data/working/MasterFieldDataFC_NSW - MIR_raw.csv")
cols_condense(mir)
dim(mir)

#### Set factors ####
dat <- dat %>% 
  mutate_at(vars(CombID, UniqueID, PrelimID, Transect, Plot, Inun), 
            list(~ factor(.))
  )
str(dat)

mir <- mir %>% 
  mutate_at(vars(CombID, UniqueID, PrelimID, Transect, Plot, Inun), 
            list(~ factor(.))
  )


