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

#### Initial facet plot for proteolysis ####
facetlabs <- c("Transect 100",
               "Transect 101",
               "Transect 102",
               "Transect 103",
               "Transect 104",
               "Transect 105",
               "Transect 106",
               "Transect 107",
               "Transect 108",
               "Transect 109")
names(facetlabs) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")

ggplot(data = dat, aes(x = Date, y = Proteolysis, colour = Plot)) +
  geom_point(aes(size = Wet)) +
  geom_line() +
  scale_colour_manual(values = c("blue", "green", "orange", "brown")) +
  scale_size(range = c(0, 5)) +
    facet_wrap( ~ Transect, ncol = 2, scales='free', labeller = labeller(
    Transect = facetlabs
  )) +
  #scale_y_continuous(limits=c(0,16)) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.title.x=element_blank()) 
  
  
