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
library(lubridate)


#### data in ####
sum <- read_csv("data/processed/summary.csv")
all <- read_csv("data/processed/ChemAll_adm_OLremPLFA.csv")

sum %<>% mutate(across(c(CombID, UniqueID, PrelimID, Transect, Plot, Inun), as.factor)) 
str(sum)

all %<>% mutate(Date = dmy(Date)) %>% 
  group_by(Transect) %>% 
  mutate(PlotPos = dense_rank(desc(RTHeight))) %>%
  ungroup() %>% 
  mutate("Sampling Period" = case_when(
    Date >= as_date("2019-03-25") & Date <= as_date("2019-03-28") ~ "Autumn 2019",
    Date >= as_date("2019-07-29") & Date <= as_date("2019-07-31") ~ "Winter 2019",
    Date >= as_date("2019-11-04") & Date <= as_date("2019-11-06") ~ "At flooding",
    Date >= as_date("2020-02-03") & Date <= as_date("2020-02-05") ~ "3 months post flood",
    Date >= as_date("2020-10-13") & Date <= as_date("2020-10-15") ~ "11 months post flood"
  ) 
  ) %>% 
  relocate("Sampling Period", .after = Date) %>% 
  relocate(PlotPos, .after = Plot) %>% 
  mutate(across(c(CombID, UniqueID, PrelimID, Transect, Plot, Inun, PlotPos, "Sampling Period"), as.factor)) 
str(all)


#### Ternary plot ####
ggtern(data=sum, aes(Sand, Silt, Clay, color = Transect)) + 
  geom_point(size = 4) +
  theme_rgbw() +
  theme_hidetitles() +
  theme(text = element_text(size=20)) +
  theme(legend.key=element_blank())

