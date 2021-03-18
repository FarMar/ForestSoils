#####################################################################################################

#### Forest soils data cleaning and analysis script                               ###################

#### mark.farrell@csiro.au        +61 8 8303 8664         17/12/2020 ################################

#####################################################################################################


#### Set working directory ####
setwd("/Users/markfarrell/OneDrive - CSIRO/Data/ForestSoils")


#### Packages ####
install.packages("PerformanceAnalytics")
install.packages("corrplot")
install.packages("RColorBrewer")
install.packages("plotrix")

library(tidyverse)
library(janitor)
library(PerformanceAnalytics)
library(corrplot)
library(RColorBrewer)
library(plotrix)

#### Import data ####
dat <- read_csv("data/working/MasterFieldDataFC_NSW - Data.csv")
spec(dat)

mir <- read_csv("data/working/MasterFieldDataFC_NSW - MIR_raw.csv")
cols_condense(mir)
dim(mir)

#### Set factors ####
dat <- dat %>% 
  mutate(across(c(CombID, UniqueID, PrelimID, Transect, Plot, Inun), as.factor))
str(dat)

mir <- mir %>% 
  mutate(across(c(CombID, UniqueID, PrelimID, Transect, Plot, Inun), as.factor))


#### Final data cleaning ####
# DON
dat <- dat %>% 
  mutate(DON2 = DTN - NH4 - NO3,
         DON = case_when(
           DON2 >= 0 ~ DON2,
           TRUE ~ 0
         )) %>% 
  select(-DON2)

#MBC, MBN
dat <- dat %>% 
  mutate(MBC = ((FumDOC - DOC) / .45),
         MBN = ((FumDTN - DTN) / .5)) %>% 
  select(-c(FumDOC, FumDTN))

#Microbial yield
dat <- dat %>% 
  mutate(MicY = AAMin_b / (AAMin_a + AAMin_b))

#C:N
dat <- dat %>% 
  mutate(CN = TotOC / TotN)

#MicC:N
dat <- dat %>% 
  mutate(MicCN = MBC / MBN)

#Vulnerability
dat <- dat %>% 
  mutate(Vuln = POC / (HOC + ROC))

#BD g cm-3 
bdmoist <- dat %>% # Need to do this as cores were only taken at 1st sampling so only moisture from that sampling relevant
  slice(1:40) %>% 
  pull(Moisture)
bdmoist = rep(bdmoist, times = 5)
dat <- dat %>% 
  mutate(bdm = bdmoist)

dat <- dat %>%
  mutate(BD0_30 = (CoreMass * (1 - bdm)) / CoreVol)

## Air dry moisture correction for all variables derived from dry soil
adm <- c("CEC", "AvailP", "ExCa", "ExMg", "ExNa",  "ExK", "Al", "As", "B", "Ca", "Cd", "Co", "Cr", "Cu", "Fe", "K", "Mg", "Mn", "Mo", "Na", "Ni", "P", "Pb", "S", 
         "Sb", "Se", "Zn", "TotOC", "TotN", "MTOC", "POC", "HOC", "ROC", "WHC", "Proteolysis")

dat_adm <- dat %>% 
  mutate(across(all_of(adm), ~ . * MoisFAD))

write_csv(dat_adm, "data/processed/ChemAll_adm.csv")


## Perform outlier removal in Excel using conditional formatting and sort by variable, then re-import and re-apply factors

OL_cor <- read_csv("data/processed/ChemAll_adm_OLrem.csv")

OL_cor <- OL_cor %>% 
  mutate(across(c(CombID, UniqueID, PrelimID, Transect, Plot, Inun), as.factor))
str(OL_cor)


## Stocks conversion - to 30 cm, per sq m - mostly either mg or g
StockVars <- c("CEC", "AvailP", "ExCa", "ExMg", "ExNa",  "ExK", "Al", "As", "B", "Ca", "Cd", "Co", "Cr", "Cu", "Fe", "K", "Mg", "Mn", "Mo", "Na", "Ni", "P", "Pb", "S", 
               "Sb", "Se", "Zn", "TotOC", "TotN", "MTOC", "POC", "HOC", "ROC", "WHC", "Proteolysis", "DOC", "DTN", "NO3", "NH4", "FAA", "DON", "MBC", "MBN")

stocks <- OL_cor %>% 
  mutate(across(all_of(StockVars), ~ . * BD0_30 * ((30 * 100 * 100) / 1000)))

write_csv(stocks, "data/processed/ChemAll_adm_stocks.csv")

#### Summarise over sampling (mean, SEM) ####
# Aim is to have a final df with the dynamic variables presented as March-19, mean, SEM alongside those only measured in March-19 e.g., CEC

tempvars <- c("NDVI",	"VH",	"VV",	"Wet", "Moisture", "pHw",	"pHc",	"EC", "AvailP", "TotOC",	"TotN",	"d13C",	"d15N",	"MTOC",	
              "POC",	"HOC",	"ROC",	"DOC",	"DTN",	"NO3",	"NH4",	"FAA", "Proteolysis",	"AAMin_k1",	"AAMin_k2",	"AAMin_a",	
              "AAMin_b",	"DON",	"MBC",	"MBN",	"MicY",	"CN",	"MicCN",	"Vuln")

summary <- OL_cor %>% 
  group_by(PrelimID) %>% 
  summarise(across(all_of(tempvars), 
                   list(mean = ~ mean(.x, na.rm = TRUE), sem = ~ std.error(.x, na.rm = TRUE)))) # {plotrix} is required for `std.error`

t1 <- OL_cor %>% 
  slice(1:40)

t1_summary <- left_join(t1, summary, by = "PrelimID")

t1_summary %>% write_csv("data/processed/summary.csv")

