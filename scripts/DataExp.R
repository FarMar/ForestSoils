#####################################################################################################

#### Forest soils data exploration                                                ###################

#### mark.farrell@csiro.au        +61 8 8303 8664         18/03/2021 ################################

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
library(lubridate)
library(ggpmisc)
library(vegan)
library(ape)
library(RVAideMemoire)
library(BiodiversityR)
library(ggbluebadge)

#### Import data ####
OL_cor <- read_csv("data/processed/ChemAll_adm_OLrem.csv")
OL_cor <- OL_cor %>% 
  mutate(across(c(CombID, UniqueID, PrelimID, Transect, Plot, Inun), as.factor)) %>% 
  mutate(Date = dmy(Date))
str(OL_cor)


t1_summary <- read_csv("data/processed/summary.csv")
t1_summary <- t1_summary %>% 
  mutate(across(c(CombID, UniqueID, PrelimID, Transect, Plot, Inun), as.factor))
str(t1_summary)
t1_summary <- t1_summary %>% 
  relocate(where(is.character))

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

ggplot(data = OL_cor, aes(x = Date, y = Proteolysis, colour = Plot)) +
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
trim <- OL_cor %>% select(UniqueID, Date, Transect, Plot, NDVI, VH,	VV,	Wet, Moisture, pHw,	pHc,	EC, AvailP, 
                       DOC,	DTN,	NO3,	NH4,	FAA, Proteolysis,	AAMin_k1,	AAMin_k2,	AAMin_a, AAMin_b,	DON,	MBC,	
                       MBN,	MicY,	MicCN)


#Names for response and explanatory vars
#https://aosmith.rbind.io/2018/08/20/automating-exploratory-plots/
response = names(trim)[5:28]
expl = names(trim)[1:9]

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


#### summary plots ####

my.formula <- y ~ x

ggplot(data = t1_summary) +
  geom_point(aes(x = RTHeight, y = Proteolysis_mean, colour = Transect), size = 4) +
  stat_smooth(aes(x = RTHeight, y = Proteolysis_mean, colour = Transect), method = lm, se = FALSE, formula = my.formula, linetype = "longdash") +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  stat_smooth(aes(x = RTHeight, y = Proteolysis_mean), method = lm, formula = my.formula, colour = "black", size = 2) +
  stat_poly_eq(formula = my.formula, 
               aes(x = RTHeight, y = Proteolysis_mean, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  theme_classic() +
  theme(strip.background = element_blank())


response2 = names(t1_summary)[15:164]
expl2 = names(t1_summary)[1:14]

response2 = set_names(response2)
response2

expl2 = set_names(expl2)
expl2

exp.fun2 = function(x, y, z1) {
  ggplot(data = t1_summary) +
    geom_point(aes(x = .data[[x]], y = .data[[y]], colour = .data[[z1]]), size = 4) +
    stat_smooth(aes(x = .data[[x]], y = .data[[y]], colour = .data[[z1]]), method = lm, formula = my.formula, se = FALSE, linetype = "longdash") +
    scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
    stat_smooth(aes(x = .data[[x]], y = .data[[y]]), method = lm, formula = my.formula, colour = "black", size = 2) +
    stat_poly_eq(formula = my.formula, 
                 aes(x = .data[[x]], y = .data[[y]], label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE) +
    theme_classic() +
    theme(strip.background = element_blank())
}

exp.fun2("RTHeight", "Proteolysis_mean", "Transect")

exp_plots2 = map(response2, ~exp.fun2("RTHeight", .x, "Transect") )


pdf("outputs/all_summaries.pdf")
exp_plots2
dev.off()


#### Create date grouping column ####
OL_cor <- OL_cor %>%
  mutate("Sampling Period" = case_when(
    Date >= as_date("2019-03-25") & Date <= as_date("2019-03-28") ~ "Autumn 2019",
    Date >= as_date("2019-07-29") & Date <= as_date("2019-07-31") ~ "Winter 2019",
    Date >= as_date("2019-11-04") & Date <= as_date("2019-11-06") ~ "At flooding",
    Date >= as_date("2020-02-03") & Date <= as_date("2020-02-05") ~ "3 months post flood",
    Date >= as_date("2020-10-13") & Date <= as_date("2020-10-15") ~ "11 months post flood"
  ) 
         ) %>% 
  relocate("Sampling Period", .after = Date) 
OL_cor$`Sampling Period` <- as.factor(OL_cor$`Sampling Period`)

str(OL_cor)
levels(OL_cor$`Sampling Period`)

OL_cor <- OL_cor %>% 
  mutate(`Sampling Period` = fct_relevel(`Sampling Period`, #remember the back-ticks (would probably have solved factor palaver too)
                                         "Autumn 2019",
                                         "Winter 2019",
                                         "At flooding",
                                         "3 months post flood",
                                         "11 months post flood"
                                         ))



#### Data selection ####
## Temporal
temporal <- OL_cor %>% 
  select(UniqueID, Date, `Sampling Period`, Transect, Plot, Easting, Northing, Height, RHeight, RTHeight, Inun,
         NDVI, VH,	VV,	Wet, Moisture, pHw,	pHc,	EC, AvailP, 
         DOC,	DTN,	NO3,	NH4,	FAA, Proteolysis,	AAMin_k1,	DON,	MBC,	
         MBN,	MicY,	MicCN)

## Biogeochem
bgc_mean <- t1_summary %>% 
  select(UniqueID, Transect, Plot, Easting, Northing, Height, RHeight, RTHeight, Inun,
         Sand, Silt, Clay, CEC, ExCa, ExMg, ExNa, ExK, Al, B, Ca, Co, Cr, Cu, Fe, K, Mg, Mn, Na, Ni, P, Pb, S, Zn,
         WHC, BD0_30, NDVI_mean, Wet_mean, Moisture_mean, pHc_mean, EC_mean, AvailP_mean, TotOC_mean, TotN_mean,
         d13C_mean, d15N_mean, POC_mean, HOC_mean, ROC_mean, DOC_mean, DTN_mean, NO3_mean, NH4_mean, FAA_mean, Proteolysis_mean,
         AAMin_k1_mean, DON_mean, MBC_mean, MBN_mean, MicY_mean)

## MIR
# MIR import

mir <- read_csv("data/working/MasterFieldDataFC_NSW - MIR_raw.csv")
cols_condense(mir)
dim(mir)

mir <- mir %>%
  mutate("Sampling Period" = case_when(
    Date >= as_date("2019-03-25") & Date <= as_date("2019-03-28") ~ "Autumn 2019",
    Date >= as_date("2019-07-29") & Date <= as_date("2019-07-31") ~ "Winter 2019",
    Date >= as_date("2019-11-04") & Date <= as_date("2019-11-06") ~ "At flooding",
    Date >= as_date("2020-02-03") & Date <= as_date("2020-02-05") ~ "3 months post flood",
    Date >= as_date("2020-10-13") & Date <= as_date("2020-10-15") ~ "11 months post flood"
  ) 
  ) %>% 
  relocate("Sampling Period", .after = Date) 
mir$`Sampling Period` <- as.factor(mir$`Sampling Period`)

str(mir)
levels(mir$`Sampling Period`)

mir <- mir %>% 
  mutate(`Sampling Period` = fct_relevel(`Sampling Period`, #remember the back-ticks (would probably have solved factor palaver too)
                                         "Autumn 2019",
                                         "Winter 2019",
                                         "At flooding",
                                         "3 months post flood",
                                         "11 months post flood"
  ))

# initial check plot
spec <- mir %>% 
  select(2, 26:1996)

waves <- seq(7999.27979, 401.121063, by = -3.8569)
colnames(spec[,2:1972]) <- waves

matplot(x = waves, 
        y = t(spec[2:1972]), 
        ylim = c(0, 3.5),
        type = "l",
        lty = 1,
        main = "Raw spectra",
        xlab = "Wavenumber (cm-1)",
        ylab = "Absorbance",
        col = rep(palette(), each = 3)
)

# Interpolation
mirinterp <- spec

mirinterp1 <- new("hyperSpec", # makes the hyperspec object
                  spc = mirinterp[, grep('[[:digit:]]', colnames(mirinterp))], 
                  wavelength = as.numeric(colnames(mirinterp)[grep	('[[:digit:]]', colnames(mirinterp))]),
                  label = list(.wavelength = "Wavenumber",
                               spc = "Intensity"))

mirinterp3 <- hyperSpec::spc.loess(mirinterp1, c(seq(6000, 600, -4)))
# plot(mirinterp3, "spc", wl.reverse = T, col = rep(palette(), each = 3))
output <- mirinterp3[[]]

waves_l <- seq(6000, 600, by = -4)
colnames(output) <- waves_l

ID <- as.data.frame(mir$UniqueID)
final <- cbind(ID, output) #This is now the re-sampled df. Still needs baselining.

matplot(x = waves_l, y = t(final[,2:1352]), ylim=c(0,3), type = "l", lty = 1, 
        main = "Absorbance - 600 to 6000 & reample with resolution of 4", xlab = "Wavelength (nm)",
        ylab = "Absorbance", col = rep(palette(), each = 3))


# baseline offset
spoffs2 <- function (spectra) 
{
  if (missing(spectra)) {
    stop("No spectral data provided")
  }
  if (spectra[1, 1] < spectra[1, dim(spectra)[2]]) {
    spectra <- t(apply(spectra, 1, rev))
  }
  s <- matrix(nrow = dim(spectra)[1], ncol = dim(spectra)[2])
  for (i in 1:dim(spectra)[1]) {
    s[i, ] <- spectra[i, ] - min(spectra[i, ])
  }
  output <- rbind(spectra[1, ], s)
  output <- output[-1,]
}

spec_a_bc_d <- spoffs2(final[,2:1352])
dim(spec_a_bc_d)
head(spec_a_bc_d)

waves_ss <- seq(600, 6000, by=4)

matplot(x = waves_ss, y = t(spec_a_bc_d), ylim=c(0,2), xlim=rev(c(600, 6000)), type = "l", lty = 1,
        main = "Absorbance - baseline corrected", xlab = expression("Wavenumber" ~ (cm^{-1})),
        ylab = "Absorbance", col = rep(palette(), each = 3))


finalb <- cbind(ID, spec_a_bc_d) %>% #This is now the baselined and re-sampled df.
  rename(UniqueID = "mir$UniqueID")

# combine data
mir_meta <- temporal %>%
  select(UniqueID, Date, `Sampling Period`, Transect, Plot, Easting, Northing, Height, RHeight, RTHeight, Inun, Moisture)

mir_proc <- left_join(mir_meta, finalb, by = "UniqueID")


#### Multivariate Exploration and Analysis ####
## MIR
# Prep
tmir <- mir_proc %>% 
  mutate(across(c(13:1363), ~((.+10)^(1/4))))

z.fn <- function(x) {
  (x-mean(x))/sd(x)
}

stmir <- tmir %>% 
  mutate(across(c(13:1363), ~z.fn(.)))

fmir <- stmir %>% 
  select(1:12)
dmir <- stmir %>% 
  select(13:1363)

distmir <- vegdist(dmir, method = "manhattan", na.rm = TRUE)
pmir <- pcoa(distmir)
pmir$values$Relative_eig[1:10]
barplot(pmir$values$Relative_eig[1:10])

mir_points <- bind_cols(fmir, (as.data.frame(pmir$vectors)))

# Plot
ggplot(mir_points) + 
  geom_point(aes(x=Axis.1, y=Axis.2, colour = Transect), size = 4) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(
    x = "PCoA Axis 1; 81.0%",
    y = "PCoA Axis 2; 7.9%")

# Permanova
set.seed(1983)
perm_mir <- adonis2(distmir~Transect*`Sampling Period`, data = stmir, permutations = 9999, method = "manhattan")
perm_mir #strong impact of transect, weak of sampling time
permpt_mir <- pairwise.perm.manova(distmir, stmir$Transect, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpt_mir
permpd_mir <- pairwise.perm.manova(distmir, stmir$`Sampling Period`, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpd_mir #sniff of significance for last sampling vs 1st three samplings

perm_mirh <- adonis2(distmir~Transect*RTHeight, data = stmir, permutations = 9999, method = "manhattan")
perm_mirh #strong height interaction

# CAP by transect

