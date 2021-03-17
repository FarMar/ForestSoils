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


library(tidyverse)
library(janitor)
library(PerformanceAnalytics)
library(corrplot)
library(RColorBrewer)


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











## Stocks conversion - to 30 cm, per sq m - mostly either mg or g
StockVars <- c("CEC", "AvailP", "ExCa", "ExMg", "ExNa",  "ExK", "Al", "As", "B", "Ca", "Cd", "Co", "Cr", "Cu", "Fe", "K", "Mg", "Mn", "Mo", "Na", "Ni", "P", "Pb", "S", 
               "Sb", "Se", "Zn", "TotOC", "TotN", "MTOC", "POC", "HOC", "ROC", "WHC", "Proteolysis", "DOC", "DTN", "NO3", "NH4", "FAA", "DON", "MBC", "MBN")

stocks <- dat_adm %>% 
  mutate(across(StockVars, ~ . * BD0_30 * ((30 * 100 * 100) / 1000)))

write_csv(stocks, "data/processed/ChemAll_adm_stocks.csv")

#### Summarise over sampling (mean, SEM) ####
# Aim is to have a final df with the dynamic variables presented as March-19, mean, SEM alongside those only measured in March-19 e.g., CEC





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
  geom_point(aes(size = Moisture)) +
  geom_line() +
  scale_colour_manual(values = brewer.pal(n = 4, name = "BrBG")) + # (values = c("blue", "dark green", "brown", "dark brown")) +
  scale_size(range = c(0, 6)) +
    facet_wrap( ~ Transect, ncol = 2, scales='free', labeller = labeller(
    Transect = facetlabs
  )) +
  #scale_y_continuous(limits=c(0,16)) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.title.x=element_blank()) 
  
  
#### Looping ####
# limit columns to just the factors required for the plot and the response variables from all five time points
trim <- dat %>% select(
  UniqueID,
  Date,
  Transect,
  Plot,
  NDVI,
  Wet,
  Moisture,
  pHw,
  pHc,
  EC,
  AvailP,
  DOC,
  DTN,
  FumDOC,
  FumDTN,
  NO3,
  NH4,
  FAA,
  Proteolysis,
  AAMin_k1,
  AAMin_k2,
  AAMin_a,
  AAMin_b
)

#Names for response and explanatory vars
#https://aosmith.rbind.io/2018/08/20/automating-exploratory-plots/
response = names(trim)[5:23]
expl = names(trim)[1:7]

response = set_names(response)
response

expl = set_names(expl)
expl

exp.fun = function(x, y, z1, z2, z3) {
  ggplot(data = trim, aes(x = .data[[x]], y = .data[[y]], colour = .data[[z2]])) +
    geom_point(aes(size = .data[[z3]])) +
    geom_line() +
    scale_colour_manual(values = brewer.pal(n = 4, name = "BrBG")) + 
    scale_size(range = c(0, 6)) +
    facet_wrap( ~ .data[[z1]], ncol = 2, scales='free') + # labeller won't work with the .data for some reason
    #scale_y_continuous(limits=c(0,16)) +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.title.x=element_blank()) 
 }

exp.fun("Date", "Proteolysis", "Transect", "Plot", "Moisture")

exp_plots = map(response, ~exp.fun("Date", .x, "Transect", "Plot", "Moisture") )

pdf("outputs/all_scatterplots.pdf")
exp_plots
dev.off()

