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
library(magrittr)

#### Import data ####
OL_cor <- read_csv("data/processed/ChemAll_adm_OLrem.csv")
OL_cor <- OL_cor %>% 
  group_by(Transect) %>% 
  mutate(PlotPos = dense_rank(desc(RTHeight))) %>%
  ungroup() %>% 
  relocate(PlotPos, .after = Plot) %>% 
  mutate(across(c(CombID, UniqueID, PrelimID, Transect, Plot, Inun, PlotPos), as.factor)) %>% 
  mutate(Date = dmy(Date)) 
str(OL_cor)


t1_summary <- read_csv("data/processed/summary.csv")
t1_summary <- t1_summary %>% 
  group_by(Transect) %>% 
  mutate(PlotPos = dense_rank(desc(RTHeight))) %>%
  ungroup() %>% 
  relocate(PlotPos, .after = Plot) %>% 
  mutate(across(c(CombID, UniqueID, PrelimID, Transect, Plot, Inun, PlotPos), as.factor))
str(t1_summary)
t1_summary <- t1_summary %>% 
  relocate(where(is.character))


plfa <- read_csv("data/working/MasterFieldDataFC_NSW - PLFAs.csv")
plfa <- plfa %>%
  mutate(Date = dmy(Date)) %>% 
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

plfa <- plfa %>% 
  mutate(`Sampling Period` = fct_relevel(`Sampling Period`, #remember the back-ticks (would probably have solved factor palaver too)
                                         "Autumn 2019",
                                         "Winter 2019",
                                         "At flooding",
                                         "3 months post flood",
                                         "11 months post flood"
  ))
  
str(plfa)

OLP_cor <- read_csv("data/processed/ChemAll_adm_OLremPLFA.csv")
OLP_cor <- OLP_cor %>%
  mutate(Date = dmy(Date)) %>% 
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
str(OLP_cor)



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

ggplot(data = OL_cor, aes(x = Date, y = Proteolysis, colour = PlotPos)) +
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
trim <- OL_cor %>% select(UniqueID, Date, Transect, Plot, PlotPos, NDVI, VH,	VV,	Wet, Moisture, pHw,	pHc,	EC, AvailP, 
                       DOC,	DTN,	NO3,	NH4,	FAA, Proteolysis,	AAMin_k1,	AAMin_k2,	AAMin_a, AAMin_b,	DON,	MBC,	
                       MBN,	MicY,	MicCN)

str(trim)
#Names for response and explanatory vars
#https://aosmith.rbind.io/2018/08/20/automating-exploratory-plots/
response = names(trim)[6:29]
expl = names(trim)[1:10]

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
  select(UniqueID, Date, `Sampling Period`, Transect, Plot, PlotPos, Easting, Northing, Height, RHeight, RTHeight, Inun,
         NDVI, VH,	VV,	Wet, Moisture,	pHc,	EC, AvailP, 
         DOC,	DTN,	NO3,	NH4,	FAA, Proteolysis,	AAMin_k1,	DON,	MBC,	
         MBN,	MicY,	MicCN)

## Biogeochem
bgc_mean <- t1_summary %>% 
  select(UniqueID, Transect, Plot, PlotPos, Easting, Northing, Height, RHeight, RTHeight, Inun,
         Clay, CEC, WHC, BD0_30, NDVI_mean, Wet_mean, Moisture_mean, pHc_mean, EC_mean, AvailP_mean, CN_mean, Vuln_mean,
         d13C_mean, d15N_mean, DOC_mean, NO3_mean, NH4_mean, FAA_mean, Proteolysis_mean,
         AAMin_k1_mean, DON_mean, MBC_mean, MBN_mean, MicY_mean)

##Temporal + PLFA
temporalP <- OLP_cor %>% 
  select(UniqueID, Date, `Sampling Period`, Transect, Plot, PlotPos, Easting, Northing, Height, RHeight, RTHeight, Inun,
         NDVI, VH,	VV,	Wet, Moisture,	pHc,	EC, AvailP, 
         DOC,	DTN,	NO3,	NH4,	FAA, Proteolysis,	AAMin_k1,	DON,	MBC,	
         MBN,	MicY,	MicCN, TotalPLFA, F_B, Gp_Gn, Act_Gp)

#### MIR ####
# MIR import

mir <- read_csv("data/working/MasterFieldDataFC_NSW - MIR_raw.csv")
cols_condense(mir)
dim(mir)

mir <- mir %>%
  group_by(Transect) %>% 
  mutate(PlotPos = dense_rank(desc(RTHeight))) %>%
  ungroup() %>% 
  relocate(PlotPos, .after = Plot) %>% 
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
  select(2, 27:1997)

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
  select(UniqueID, Date, `Sampling Period`, Transect, Plot, PlotPos, Easting, Northing, Height, RHeight, RTHeight, Inun, Moisture)

mir_proc <- left_join(mir_meta, finalb, by = "UniqueID")


## Multivariate Exploration and Analysis
## MIR
# Prep
tmir <- mir_proc %>% 
  mutate(across(c(14:1364), ~((.+10)^(1/4))))

z.fn <- function(x) {
  (x-mean(x))/sd(x)
}

stmir <- tmir %>% 
  mutate(across(c(14:1364), ~z.fn(.)))

fmir <- stmir %>% 
  select(1:13)
dmir <- stmir %>% 
  select(14:1363)

distmir <- vegdist(dmir, method = "manhattan", na.rm = TRUE)
pmir <- pcoa(distmir)
pmir$values$Relative_eig[1:10]
barplot(pmir$values$Relative_eig[1:10])

mir_points <- bind_cols(fmir, (as.data.frame(pmir$vectors)))

# Plot
ggplot(mir_points) + 
  geom_point(aes(x=Axis.1, y=Axis.2, colour = Transect, shape = PlotPos), size = 4) +
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
stmir <- as.data.frame(stmir) 
cap_mirt <- CAPdiscrim(distmir~Transect, data = stmir, axes = 10, m = 0, mmax = 10, add = FALSE, permutations = 999)
round(cap_mirt$F/sum(cap_mirt$F), digits=3)
barplot(cap_mirt$F/sum(cap_mirt$F))

cap_mirt_points <- bind_cols((as.data.frame(cap_mirt$x)), fmir) 
glimpse(cap_mirt_points)

ggplot(cap_mirt_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = Transect, shape = PlotPos), size = 4) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(
    x = "CAP Axis 1; 41.2%",
    y = "CAP Axis 2; 35.3%")

# CAP + spider
mir_cent <- aggregate(cbind(LD1, LD2) ~ Transect, data = cap_mirt_points, FUN = mean)

mir_segs <- merge(cap_mirt_points, setNames(cent, c('Transect', 'oLD1', 'oLD2')), by = 'Transect', sort = FALSE)

ggplot(cap_mirt_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = Transect, shape = PlotPos), size = 3, alpha = .7) +
  geom_segment(data = mir_segs, mapping = aes(x = LD1, y = LD2, xend = oLD1, yend = oLD2, colour = Transect), alpha = .7, size = .25) +
  geom_point(data = mir_cent, mapping = aes(x = LD1, y = LD2, colour = Transect), size = 5) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(
    x = "CAP Axis 1; 41.2%",
    y = "CAP Axis 2; 35.3%")



#### BGC ####
#pre-prep - PCA of total emlements to reduce dimenstions
tot_elms <- t1_summary %>% 
  select(47:66) %>% 
  select(!c(As, B, Cd, Mo, Sb, Se))
chart.Correlation(tot_elms)

ttot_elms <- tot_elms %>% 
  mutate(P = log1p(P),
         Na = log1p(Na),
         Mg = log1p(Mg),
         K = log1p(K),
         Co = log1p(Co),
         Ca = log1p(Ca))
chart.Correlation(ttot_elms)

pca_elms <- princomp(ttot_elms, cor = TRUE, scores = TRUE)
biplot(pca_elms, choices = c(1,2))
summary(pca_elms) #PC1 = 59.2%, PC2 = 11.7%
scores_elms <- as.data.frame(pca_elms[["scores"]]) %>% 
  select(1:2)

#prep
bgc_mean <- cbind(bgc_mean, scores_elms)
bgc_cor <- select(bgc_mean, 11:36)
chart.Correlation(bgc_cor, histogram=TRUE, pch=19)


tbgc_mean <- bgc_mean %>% 
  mutate(MBN_mean = log1p(MBN_mean),
         NH4_mean = log1p(NH4_mean),
         AvailP_mean = log1p(AvailP_mean),
         EC_mean = log1p(EC_mean),
         pHc_mean = log1p(pHc_mean),
         BD0_30 = log1p(BD0_30))

stbgc_mean <- tbgc_mean %>% 
  mutate(across(c(11:36), ~z.fn(.)))

fbgc <- stbgc_mean %>% 
  select(1:10)
dbgc <- stbgc_mean %>% 
  select(11:36)

# PCoA
distbgc <- vegdist(dbgc, method = "euclidean", na.rm = TRUE)
pbgc <- pcoa(distbgc)
pbgc$values$Relative_eig[1:10]
barplot(pbgc$values$Relative_eig[1:10])

bgc_points <- bind_cols(fbgc, (as.data.frame(pbgc$vectors)))

compute.arrows = function (given_pcoa, orig_df) {
  orig_df = orig_df #can be changed to select columns of interest only
  n <- nrow(orig_df)
  points.stand <- scale(given_pcoa$vectors)
  S <- cov(orig_df, points.stand) #compute covariance of variables with all axes
  pos_eigen = given_pcoa$values$Eigenvalues[seq(ncol(S))] #select only +ve eigenvalues
  U <- S %*% diag((pos_eigen/(n - 1))^(-0.5)) #Standardise value of covariance
  colnames(U) <- colnames(given_pcoa$vectors) #Get column names
  given_pcoa$U <- U #Add values of covariates inside object
  return(given_pcoa)
}
pbgc = compute.arrows(pbgc, dbgc)

pbgc_arrows_df <- as.data.frame(pbgc$U*10) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")

# Plot
ggplot(bgc_points) + 
  geom_point(aes(x=Axis.1, y=Axis.2, colour = Transect, shape = PlotPos), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = pbgc_arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = pbgc_arrows_df, aes(x=Axis.1, y=Axis.2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "PCoA Axis 1; 25.6%",
    y = "PCoA Axis 2; 16.2%")

# Permanova
set.seed(1983)
perm_bgc <- adonis2(distbgc~Transect+PlotPos, data = stbgc_mean, permutations = 9999, method = "euclidean")
perm_bgc #strong impact of transect and plot
permpt_bgc <- pairwise.perm.manova(distbgc, stbgc_mean$Transect, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpt_bgc #.098 is lowest possible  - several pairwise comps have this
permpp_bgc <- pairwise.perm.manova(distbgc, stbgc_mean$PlotPos, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpp_bgc #4 is sig diff from 1&2. 3 borderline diff from 1&2. 1 borderline diff from 2


# CAP by transect
stbgc_mean <- as.data.frame(stbgc_mean) 
cap_bgct <- CAPdiscrim(distbgc~Transect, data = stbgc_mean, axes = 10, m = 0, mmax = 10, add = FALSE, permutations = 999)
cap_bgct <- add.spec.scores(cap_bgct, dbgc, method = "cor.scores", multi = 1, Rscale = F, scaling = "1")
round(cap_bgct$F/sum(cap_bgct$F), digits=3)
barplot(cap_bgct$F/sum(cap_bgct$F))

cap_bgct_points <- bind_cols((as.data.frame(cap_bgct$x)), fbgc) 
glimpse(cap_bgct_points)

cap_bgct_arrows <- as.data.frame(cap_bgct$cproj*5) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")

ggplot(cap_bgct_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = Transect, shape = PlotPos), size = 4) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_bgct_arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_bgct_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 56.7%",
    y = "CAP Axis 2; 23.0%")

# CAP by transect + spider
bgc_centt <- aggregate(cbind(LD1, LD2) ~ Transect, data = cap_bgct_points, FUN = mean)

bgc_segst <- merge(cap_bgct_points, setNames(bgc_centt, c('Transect', 'oLD1', 'oLD2')), by = 'Transect', sort = FALSE)

ggplot(cap_bgct_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = Transect, shape = PlotPos), size = 3, alpha = .6) +
  geom_segment(data = bgc_segst, mapping = aes(x = LD1, y = LD2, xend = oLD1, yend = oLD2, colour = Transect), alpha = .7, size = .25) +
  geom_point(data = bgc_centt, mapping = aes(x = LD1, y = LD2, colour = Transect), size = 5) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_bgct_arrows,
               x = 0, y = 0, alpha = 0.3,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_bgct_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 56.7%",
    y = "CAP Axis 2; 23.0%")


# CAP by plotpos
stbgc_mean <- as.data.frame(stbgc_mean) 
cap_bgcp <- CAPdiscrim(distbgc~PlotPos, data = stbgc_mean, axes = 10, m = 3, mmax = 10, add = FALSE, permutations = 999)
cap_bgcp <- add.spec.scores(cap_bgcp, dbgc, method = "cor.scores", multi = 1, Rscale = F, scaling = "1")
round(cap_bgcp$F/sum(cap_bgcp$F), digits=3)
barplot(cap_bgcp$F/sum(cap_bgcp$F))

cap_bgcp_points <- bind_cols((as.data.frame(cap_bgcp$x)), fbgc) 
glimpse(cap_bgcp_points)

cap_bgcp_arrows <- as.data.frame(cap_bgcp$cproj*3) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")

ggplot(cap_bgcp_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = PlotPos), size = 4) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_bgcp_arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_bgcp_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 76.3%",
    y = "CAP Axis 2; 23.7%")

# CAP by plot + spider
bgc_centp <- aggregate(cbind(LD1, LD2) ~ Plot, data = cap_bgcp_points, FUN = mean)

bgc_segsp <- merge(cap_bgcp_points, setNames(bgc_centp, c('Plot', 'oLD1', 'oLD2')), by = 'Plot', sort = FALSE)

ggplot(cap_bgcp_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = Plot), size = 3, alpha = .6) +
  geom_segment(data = bgc_segsp, mapping = aes(x = LD1, y = LD2, xend = oLD1, yend = oLD2, colour = Plot), alpha = .9, size = .3) +
  geom_point(data = bgc_centp, mapping = aes(x = LD1, y = LD2, colour = Plot), size = 5) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_bgcp_arrows,
               x = 0, y = 0, alpha = 0.3,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_bgcp_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 56.7%",
    y = "CAP Axis 2; 23.0%")

#### Temporal ####
# Data for this are in `temporal`
glimpse(temporal)
temporal %<>% relocate(Inun, .after = PlotPos)
temporal <- temporal %>% 
  mutate(Inun = fct_relevel(`Inun`,
                            "y",
                            "m",
                            "n"))

# Quick correlation plot for evaluation
chart.Correlation(temporal[, 8:32], histogram = TRUE, pch = 19)

# Drop and transform
ttemporal <- temporal %>% 
  select(-c(VH, VV, DTN)) %>% 
  mutate(across(c(Moisture, pHc, EC, AvailP, NO3, NH4, FAA, Proteolysis, DON, MBC, MBN, MicCN), ~log1p(.)))
  
chart.Correlation(ttemporal[, 8:29], histogram = TRUE, pch = 19)

#prep
sttemporal <- ttemporal %>% 
  drop_na() %>% 
  mutate(across(c(13:29), ~z.fn(.)))

ftemp <- sttemporal %>% 
  select(1:12)
dtemp <- sttemporal %>% 
  select(13:29)

#PCoA
disttemp <- vegdist(dtemp, method = "euclidean", na.rm = TRUE)
ptemp <- pcoa(disttemp)
ptemp$values$Relative_eig[1:10]
barplot(ptemp$values$Relative_eig[1:10])

temp_points <- bind_cols(ftemp, (as.data.frame(ptemp$vectors)))

compute.arrows = function (given_pcoa, orig_df) {
  orig_df = orig_df #can be changed to select columns of interest only
  n <- nrow(orig_df)
  points.stand <- scale(given_pcoa$vectors)
  S <- cov(orig_df, points.stand) #compute covariance of variables with all axes
  pos_eigen = given_pcoa$values$Eigenvalues[seq(ncol(S))] #select only +ve eigenvalues
  U <- S %*% diag((pos_eigen/(n - 1))^(-0.5)) #Standardise value of covariance
  colnames(U) <- colnames(given_pcoa$vectors) #Get column names
  given_pcoa$U <- U #Add values of covariates inside object
  return(given_pcoa)
}
ptemp = compute.arrows(ptemp, dtemp)

ptemp_arrows_df <- as.data.frame(ptemp$U*10) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")

# Plot
ggplot(temp_points) + #Some separation by date, transect# seems noisy
  geom_point(aes(x=Axis.1, y=Axis.2, colour = Transect, shape = `Sampling Period`), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = ptemp_arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = ptemp_arrows_df, aes(x=Axis.1, y=Axis.2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "PCoA Axis 1; 19.7%",
    y = "PCoA Axis 2; 17.0%")

ggplot(temp_points) + #A bit more informative, definite axis1 trend of transect. Date clustering a bit more obvious
  geom_point(aes(x=Axis.1, y=Axis.2, colour = PlotPos, shape = `Sampling Period`), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = ptemp_arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = ptemp_arrows_df, aes(x=Axis.1, y=Axis.2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "PCoA Axis 1; 19.7%",
    y = "PCoA Axis 2; 17.0%")

ggplot(temp_points) + #Seems to clearly show separation
  geom_point(aes(x=Axis.1, y=Axis.2, colour = PlotPos, shape = Inun), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  scale_shape_manual(values = c(15, 18, 0)) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = ptemp_arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = ptemp_arrows_df, aes(x=Axis.1, y=Axis.2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "PCoA Axis 1; 19.7%",
    y = "PCoA Axis 2; 17.0%")


# Permanova
set.seed(1983)
perm_temptp <- adonis2(disttemp~Transect*`Sampling Period`, data = sttemporal, permutations = 9999, method = "euclidean")
perm_temptp #strong impact of transect and sampling period, no interaction
perm_temppp <- adonis2(disttemp~PlotPos*`Sampling Period`, data = sttemporal, permutations = 9999, method = "euclidean")
perm_temppp #strong impact of plot position and sampling period, no interaction
perm_temptpp <- adonis2(disttemp~Transect+PlotPos+`Sampling Period`, data = sttemporal, permutations = 9999, method = "euclidean")
perm_temptpp #strong impact of transect, plot position and sampling period in additive model
permpt_temp <- pairwise.perm.manova(disttemp, sttemporal$Transect, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpt_temp #All differ except 0&8, 3&9, 5&7
permpp_temp <- pairwise.perm.manova(disttemp, sttemporal$PlotPos, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpp_temp #All differ except 2&3
permpp_temp <- pairwise.perm.manova(disttemp, sttemporal$`Sampling Period`, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpp_temp #All differ

# CAP by transect
sttemporal <- as.data.frame(sttemporal) 
cap_tempt <- CAPdiscrim(disttemp~Transect, data = sttemporal, axes = 10, m = 0, mmax = 10, add = FALSE, permutations = 9)
cap_tempt <- add.spec.scores(cap_tempt, dtemp, method = "cor.scores", multi = 1, Rscale = F, scaling = "1")
round(cap_tempt$F/sum(cap_tempt$F), digits=3)
barplot(cap_tempt$F/sum(cap_tempt$F))

cap_tempt_points <- bind_cols((as.data.frame(cap_tempt$x)), ftemp) 
glimpse(cap_tempt_points)

cap_tempt_arrows <- as.data.frame(cap_tempt$cproj*5) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")

ggplot(cap_tempt_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = Transect, shape = PlotPos), size = 4) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_tempt_arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_tempt_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 62.1%",
    y = "CAP Axis 2; 18.5%")

# CAP by transect + spider
temp_centt <- aggregate(cbind(LD1, LD2) ~ Transect, data = cap_tempt_points, FUN = mean)

temp_segst <- merge(cap_tempt_points, setNames(temp_centt, c('Transect', 'oLD1', 'oLD2')), by = 'Transect', sort = FALSE)

ggplot(cap_tempt_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = Transect, shape = PlotPos), size = 3, alpha = .6) +
  geom_segment(data = temp_segst, mapping = aes(x = LD1, y = LD2, xend = oLD1, yend = oLD2, colour = Transect), alpha = .7, size = .25) +
  geom_point(data = temp_centt, mapping = aes(x = LD1, y = LD2, colour = Transect), size = 5) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_tempt_arrows,
               x = 0, y = 0, alpha = 0.3,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_tempt_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 62.1%",
    y = "CAP Axis 2; 18.5%")


# CAP by plotpos
cap_tempp <- CAPdiscrim(disttemp~PlotPos, data = sttemporal, axes = 10, m = 0, mmax = 10, add = FALSE, permutations = 9)
cap_tempp <- add.spec.scores(cap_tempp, dtemp, method = "cor.scores", multi = 1, Rscale = F, scaling = "1")
round(cap_tempp$F/sum(cap_tempp$F), digits=3)
barplot(cap_tempp$F/sum(cap_tempp$F))

cap_tempp_points <- bind_cols((as.data.frame(cap_tempp$x)), ftemp) 
glimpse(cap_tempp_points)

cap_tempp_arrows <- as.data.frame(cap_tempp$cproj*5) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")

ggplot(cap_tempp_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = PlotPos), size = 4) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_tempp_arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_tempp_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 79.5%",
    y = "CAP Axis 2; 20.0%")

# CAP by plot + spider
temp_centp <- aggregate(cbind(LD1, LD2) ~ PlotPos, data = cap_tempp_points, FUN = mean)

temp_segsp <- merge(cap_tempp_points, setNames(temp_centp, c('PlotPos', 'oLD1', 'oLD2')), by = 'PlotPos', sort = FALSE)

ggplot(cap_tempp_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = PlotPos), size = 3, alpha = .6) +
  geom_segment(data = temp_segsp, mapping = aes(x = LD1, y = LD2, xend = oLD1, yend = oLD2, colour = PlotPos), alpha = .9, size = .3) +
  geom_point(data = temp_centp, mapping = aes(x = LD1, y = LD2, colour = PlotPos), size = 5) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_tempp_arrows,
               x = 0, y = 0, alpha = 0.3,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_tempp_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 79.5%",
    y = "CAP Axis 2; 20.0%")

# CAP by SamplingPeriod
cap_tempps <- CAPdiscrim(disttemp~`Sampling Period`, data = sttemporal, axes = 10, m = 0, mmax = 10, add = FALSE, permutations = 9)
cap_tempps <- add.spec.scores(cap_tempps, dtemp, method = "cor.scores", multi = 1, Rscale = F, scaling = "1")
round(cap_tempps$F/sum(cap_tempps$F), digits=3)
barplot(cap_tempps$F/sum(cap_tempps$F))

cap_tempps_points <- bind_cols((as.data.frame(cap_tempps$x)), ftemp) 
glimpse(cap_tempps_points)

cap_tempps_arrows <- as.data.frame(cap_tempps$cproj*5) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")

ggplot(cap_tempps_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = `Sampling Period`), size = 4) +
  scale_colour_manual(values = brewer.pal(n = 6, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_tempps_arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_tempps_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 66.8%",
    y = "CAP Axis 2; 21.3%")

# CAP by SamplingPeriod + spider
temp_centps <- aggregate(cbind(LD1, LD2) ~ `Sampling Period`, data = cap_tempps_points, FUN = mean)

temp_segsps <- merge(cap_tempps_points, setNames(temp_centps, c('Sampling Period', 'oLD1', 'oLD2')), by = 'Sampling Period', sort = FALSE)

ggplot(cap_tempps_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = `Sampling Period`), size = 3, alpha = .6) +
  geom_segment(data = temp_segsps, mapping = aes(x = LD1, y = LD2, xend = oLD1, yend = oLD2, colour = `Sampling Period`), alpha = .9, size = .3) +
  geom_point(data = temp_centps, mapping = aes(x = LD1, y = LD2, colour = `Sampling Period`), size = 5) +
  scale_colour_manual(values = brewer.pal(n = 6, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_tempps_arrows,
               x = 0, y = 0, alpha = 0.3,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_tempps_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 66.8%",
    y = "CAP Axis 2; 21.3%")



#### PLFAs ####
# Data for this are in `plfa`
glimpse(plfa) #remember these have been standardised by total plfa already
plfa %<>% relocate(Inun, .after = PlotPos)
plfa <- plfa %>% 
  mutate(Inun = fct_relevel(`Inun`,
                            "y",
                            "m",
                            "n"))

# Quick correlation plot for evaluation
chart.Correlation(plfa[, 15:41], histogram = TRUE, pch = 19)

stplfa1 <- filter(plfa, `Sampling Period` == "Autumn 2019")
stplfa2 <- filter(plfa, `Sampling Period` == "Winter 2019")
stplfa3 <- filter(plfa, `Sampling Period` == "At flooding")
stplfa4 <- filter(plfa, `Sampling Period` == "3 months post flood")
stplfa5 <- filter(plfa, `Sampling Period` == "11 months post flood")

#prep
stplfa <- plfa 
fplfa <- stplfa %>% 
  select(1:14)
dplfa <- stplfa %>% 
  select(15:41)

fplfa1 <- stplfa1 %>% 
  select(1:14)
dplfa1 <- stplfa1 %>% 
  select(15:41)

fplfa2 <- stplfa2 %>% 
  select(1:14)
dplfa2 <- stplfa2 %>% 
  select(15:41)

fplfa3 <- stplfa3 %>% 
  select(1:14)
dplfa3 <- stplfa3 %>% 
  select(15:41)

fplfa4 <- stplfa4 %>% 
  select(1:14)
dplfa4 <- stplfa4 %>% 
  select(15:41)

fplfa5 <- stplfa5 %>% 
  select(1:14)
dplfa5 <- stplfa5 %>% 
  select(15:41)


#PCoA
distplfa <- vegdist(dplfa, method = "euclidean", na.rm = TRUE)
distplfa1 <- vegdist(dplfa1, method = "euclidean", na.rm = TRUE)
distplfa2 <- vegdist(dplfa2, method = "euclidean", na.rm = TRUE)
distplfa3 <- vegdist(dplfa3, method = "euclidean", na.rm = TRUE)
distplfa4 <- vegdist(dplfa4, method = "euclidean", na.rm = TRUE)
distplfa5 <- vegdist(dplfa5, method = "euclidean", na.rm = TRUE)

pplfa <- pcoa(distplfa)
pplfa1 <- pcoa(distplfa1)
pplfa2 <- pcoa(distplfa2)
pplfa3 <- pcoa(distplfa3)
pplfa4 <- pcoa(distplfa4)
pplfa5 <- pcoa(distplfa5)

pplfa$values$Relative_eig[1:10]
pplfa1$values$Relative_eig[1:10]
pplfa2$values$Relative_eig[1:10]
pplfa3$values$Relative_eig[1:10]
pplfa4$values$Relative_eig[1:10]
pplfa5$values$Relative_eig[1:10]

barplot(pplfa$values$Relative_eig[1:10])

plfa_points <- bind_cols(fplfa, (as.data.frame(pplfa$vectors)))
plfa_points1 <- bind_cols(fplfa1, (as.data.frame(pplfa1$vectors)))
plfa_points2 <- bind_cols(fplfa2, (as.data.frame(pplfa2$vectors)))
plfa_points3 <- bind_cols(fplfa3, (as.data.frame(pplfa3$vectors)))
plfa_points4 <- bind_cols(fplfa4, (as.data.frame(pplfa4$vectors)))
plfa_points5 <- bind_cols(fplfa5, (as.data.frame(pplfa5$vectors)))

compute.arrows = function (given_pcoa, orig_df) {
  orig_df = orig_df #can be changed to select columns of interest only
  n <- nrow(orig_df)
  points.stand <- scale(given_pcoa$vectors)
  S <- cov(orig_df, points.stand) #compute covariance of variables with all axes
  pos_eigen = given_pcoa$values$Eigenvalues[seq(ncol(S))] #select only +ve eigenvalues
  U <- S %*% diag((pos_eigen/(n - 1))^(-0.5)) #Standardise value of covariance
  colnames(U) <- colnames(given_pcoa$vectors) #Get column names
  given_pcoa$U <- U #Add values of covariates inside object
  return(given_pcoa)
}
pplfa = compute.arrows(pplfa, dplfa)
pplfa1 = compute.arrows(pplfa1, dplfa1)
pplfa2 = compute.arrows(pplfa2, dplfa2)
pplfa3 = compute.arrows(pplfa3, dplfa3)
pplfa4 = compute.arrows(pplfa4, dplfa4)
pplfa5 = compute.arrows(pplfa5, dplfa5)

pplfa_arrows_df <- as.data.frame(pplfa$U*.2) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")
pplfa1_arrows_df <- as.data.frame(pplfa1$U*.05) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")
pplfa2_arrows_df <- as.data.frame(pplfa2$U*.1) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")
pplfa3_arrows_df <- as.data.frame(pplfa3$U*.2) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")
pplfa4_arrows_df <- as.data.frame(pplfa4$U*.2) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")
pplfa5_arrows_df <- as.data.frame(pplfa5$U*.2) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")

# Plot
ggplot(plfa_points) + 
  geom_point(aes(x=Axis.1, y=Axis.2, colour = Transect, shape = PlotPos), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = pplfa_arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = pplfa_arrows_df, aes(x=Axis.1, y=Axis.2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "PCoA Axis 1; 29.2%",
    y = "PCoA Axis 2; 17.7%")

ggplot(plfa_points1) + 
  geom_point(aes(x=Axis.1, y=Axis.2, colour = Transect, shape = PlotPos), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = pplfa1_arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = pplfa1_arrows_df, aes(x=Axis.1, y=Axis.2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "PCoA Axis 1; 41.6%",
    y = "PCoA Axis 2; 25.0%")

ggplot(plfa_points2) + 
  geom_point(aes(x=Axis.1, y=Axis.2, colour = Transect, shape = PlotPos), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = pplfa2_arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = pplfa2_arrows_df, aes(x=Axis.1, y=Axis.2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "PCoA Axis 1; 33.0%",
    y = "PCoA Axis 2; 20.5%")

ggplot(plfa_points3) + 
  geom_point(aes(x=Axis.1, y=Axis.2, colour = Transect, shape = PlotPos), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = pplfa3_arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = pplfa3_arrows_df, aes(x=Axis.1, y=Axis.2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "PCoA Axis 1; 36.1%",
    y = "PCoA Axis 2; 21.2%")

ggplot(plfa_points4) + 
  geom_point(aes(x=Axis.1, y=Axis.2, colour = Transect, shape = PlotPos), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = pplfa4_arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = pplfa4_arrows_df, aes(x=Axis.1, y=Axis.2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "PCoA Axis 1; 40.3%",
    y = "PCoA Axis 2; 21.9%")

ggplot(plfa_points5) + 
  geom_point(aes(x=Axis.1, y=Axis.2, colour = Transect, shape = PlotPos), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = pplfa5_arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = pplfa5_arrows_df, aes(x=Axis.1, y=Axis.2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "PCoA Axis 1; 35.9%",
    y = "PCoA Axis 2; 19.2%")

# Not exactly much clear here

ggplot(plfa_points) + 
  geom_point(aes(x=Axis.1, y=Axis.2, colour = PlotPos, shape = `Sampling Period`), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = pplfa_arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = pplfa_arrows_df, aes(x=Axis.1, y=Axis.2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "PCoA Axis 1; 29.2%",
    y = "PCoA Axis 2; 17.7%")

ggplot(plfa_points) + 
  geom_point(aes(x=Axis.1, y=Axis.2, colour = PlotPos, shape = Inun), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  scale_shape_manual(values = c(15, 18, 0)) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = pplfa_arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = pplfa_arrows_df, aes(x=Axis.1, y=Axis.2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "PCoA Axis 1; 29.2%",
    y = "PCoA Axis 2; 17.7%")


# Permanova
set.seed(1983)
perm_plfatp <- adonis2(distplfa~Transect*`Sampling Period`, data = stplfa, permutations = 9999, method = "euclidean")
perm_plfatp #strong impact of transect and sampling period, STRONG interaction
perm_plfapp <- adonis2(distplfa~PlotPos*`Sampling Period`, data = stplfa, permutations = 9999, method = "euclidean")
perm_plfapp #strong impact of plot position and sampling period, no interaction
perm_plfatpp <- adonis2(distplfa~Transect*`Sampling Period`+PlotPos, data = stplfa, permutations = 9999, method = "euclidean")
perm_plfatpp #strong impact of transect, plot position, date and transect*date

permpt_plfa <- pairwise.perm.manova(distplfa, stplfa$Transect, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpt_plfa #Driffer: 0&3,9; 1&9; 2&3,8,9; 3&4,6,7,9; 4&8,9; 6&9  
permpp_plfa <- pairwise.perm.manova(distplfa, stplfa$PlotPos, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpp_plfa #Differ: 1&2,3,4; 2&4
permpp_plfa <- pairwise.perm.manova(distplfa, stplfa$`Sampling Period`, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpp_plfa #1&2,3,4,5; 2&3,4,5

# Split by sampling period required
perm_plfa1tp <- adonis2(distplfa1~Transect+PlotPos, data = stplfa1, permutations = 9999, method = "euclidean")
perm_plfa1tp 
permpt_plfa1 <- pairwise.perm.manova(distplfa1, stplfa1$Transect, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpt_plfa1 #NS
permpp_plfa1 <- pairwise.perm.manova(distplfa1, stplfa1$PlotPos, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpp_plfa1 #1 differs from 2 & 4 assuming p <0.1

perm_plfa2tp <- adonis2(distplfa2~Transect+PlotPos, data = stplfa2, permutations = 9999, method = "euclidean")
perm_plfa2tp 
permpt_plfa2 <- pairwise.perm.manova(distplfa2, stplfa2$Transect, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpt_plfa2 #0 differs from all but site 1
permpp_plfa2 <- pairwise.perm.manova(distplfa2, stplfa2$PlotPos, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpp_plfa2 #NS

perm_plfa3tp <- adonis2(distplfa3~Transect+PlotPos, data = stplfa3, permutations = 9999, method = "euclidean")
perm_plfa3tp 
permpt_plfa3 <- pairwise.perm.manova(distplfa3, stplfa3$Transect, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpt_plfa3 #NS
permpp_plfa3 <- pairwise.perm.manova(distplfa3, stplfa3$PlotPos, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpp_plfa3 #1&2,3,4; 2&4; 3&4

perm_plfa4tp <- adonis2(distplfa4~Transect+PlotPos, data = stplfa4, permutations = 9999, method = "euclidean")
perm_plfa4tp 
permpt_plfa4 <- pairwise.perm.manova(distplfa4, stplfa4$Transect, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpt_plfa4 #0&3,8; 1&5,7,8; 2&3,8; 3&4,5,6,7,8,9; 4&8
permpp_plfa4 <- pairwise.perm.manova(distplfa4, stplfa4$PlotPos, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpp_plfa4 #1&2,3,4

perm_plfa5tp <- adonis2(distplfa5~Transect+PlotPos, data = stplfa5, permutations = 9999, method = "euclidean")
perm_plfa5tp 
permpt_plfa5 <- pairwise.perm.manova(distplfa5, stplfa5$Transect, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpt_plfa5 #??


# CAP by transect - Does not make a huge amount of sense going off the above. 
# Instead will incorporate total PLFA, F:B, G+:G-, G+:Actino into temporal data and re-run those analyses tomorrow

#### Temporal + PLFA ####
# Data for this are in `temporalP`
glimpse(temporalP)
temporalP %<>% relocate(Inun, .after = PlotPos)
temporalP <- temporalP %>% 
  mutate(Inun = fct_relevel(`Inun`,
                            "y",
                            "m",
                            "n"))

# Quick correlation plot for evaluation
chart.Correlation(temporalP[, 8:36], histogram = TRUE, pch = 19)

# Drop and transform
ttemporalP <- temporalP %>% 
  select(-c(VH, VV, DTN)) %>% 
  mutate(across(c(Moisture, pHc, EC, AvailP, NO3, NH4, FAA, Proteolysis, DON, MBC, MBN, MicCN, TotalPLFA, F_B), ~log1p(.)))

chart.Correlation(ttemporalP[, 8:33], histogram = TRUE, pch = 19)

#prep
sttemporalP <- ttemporalP %>% 
  drop_na() %>% 
  mutate(across(c(13:33), ~z.fn(.)))

ftempP <- sttemporalP %>% 
  select(1:12)
dtempP <- sttemporalP %>% 
  select(13:33)

#PCoA
disttempP <- vegdist(dtempP, method = "euclidean", na.rm = TRUE)
ptempP <- pcoa(disttempP)
ptempP$values$Relative_eig[1:10]
barplot(ptempP$values$Relative_eig[1:10])

tempP_points <- bind_cols(ftempP, (as.data.frame(ptempP$vectors)))

compute.arrows = function (given_pcoa, orig_df) {
  orig_df = orig_df #can be changed to select columns of interest only
  n <- nrow(orig_df)
  points.stand <- scale(given_pcoa$vectors)
  S <- cov(orig_df, points.stand) #compute covariance of variables with all axes
  pos_eigen = given_pcoa$values$Eigenvalues[seq(ncol(S))] #select only +ve eigenvalues
  U <- S %*% diag((pos_eigen/(n - 1))^(-0.5)) #Standardise value of covariance
  colnames(U) <- colnames(given_pcoa$vectors) #Get column names
  given_pcoa$U <- U #Add values of covariates inside object
  return(given_pcoa)
}
ptempP = compute.arrows(ptempP, dtempP)

ptempP_arrows_df <- as.data.frame(ptempP$U*10) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")

# Plot
ggplot(tempP_points) + #Some separation by date, transect# seems noisy
  geom_point(aes(x=Axis.1, y=Axis.2, colour = Transect, shape = `Sampling Period`), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = ptempP_arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = ptempP_arrows_df, aes(x=Axis.1, y=Axis.2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "PCoA Axis 1; 18.6%",
    y = "PCoA Axis 2; 15.7%")

ggplot(tempP_points) + #A bit more informative, definite axis1 trend of transect. Date clustering a bit more obvious
  geom_point(aes(x=Axis.1, y=Axis.2, colour = PlotPos, shape = `Sampling Period`), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = ptempP_arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = ptempP_arrows_df, aes(x=Axis.1, y=Axis.2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "PCoA Axis 1; 18.6%",
    y = "PCoA Axis 2; 15.7%")

ggplot(tempP_points) + #Seems to clearly show separation
  geom_point(aes(x=Axis.1, y=Axis.2, colour = PlotPos, shape = Inun), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  scale_shape_manual(values = c(15, 18, 0)) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = ptempP_arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = ptempP_arrows_df, aes(x=Axis.1, y=Axis.2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "PCoA Axis 1; 18.6%",
    y = "PCoA Axis 2; 15.7%")


# Permanova
set.seed(1983)
perm_tempPtp <- adonis2(disttempP~Transect*`Sampling Period`, data = sttemporalP, permutations = 9999, method = "euclidean")
perm_tempPtp #strong impact of transect and sampling period, no interaction
perm_tempPpp <- adonis2(disttempP~PlotPos*`Sampling Period`, data = sttemporalP, permutations = 9999, method = "euclidean")
perm_tempPpp #strong impact of plot position and sampling period, no interaction
perm_tempPtpp <- adonis2(disttempP~Transect+PlotPos+`Sampling Period`, data = sttemporalP, permutations = 9999, method = "euclidean")
perm_tempPtpp #strong impact of transect, plot position and sampling period in additive model
permpt_tempP <- pairwise.perm.manova(disttempP, sttemporalP$Transect, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpt_tempP #All differ except 0&8, 1&8, 3&9, 5&7
permpp_tempP <- pairwise.perm.manova(disttempP, sttemporalP$PlotPos, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpp_tempP #All differ except 2&3
permps_tempP <- pairwise.perm.manova(disttempP, sttemporalP$`Sampling Period`, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permps_tempP #All differ

# CAP by transect
sttemporalP <- as.data.frame(sttemporalP) 
cap_temptP <- CAPdiscrim(disttempP~Transect, data = sttemporalP, axes = 10, m = 0, mmax = 10, add = FALSE, permutations = 9)
cap_temptP <- add.spec.scores(cap_temptP, dtempP, method = "cor.scores", multi = 1, Rscale = F, scaling = "1")
round(cap_temptP$F/sum(cap_temptP$F), digits=3)
barplot(cap_temptP$F/sum(cap_temptP$F))

cap_temptP_points <- bind_cols((as.data.frame(cap_temptP$x)), ftempP) 
glimpse(cap_temptP_points)

cap_temptP_arrows <- as.data.frame(cap_temptP$cproj*5) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")

ggplot(cap_temptP_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = Transect, shape = PlotPos), size = 4) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_temptP_arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_temptP_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 57.0%",
    y = "CAP Axis 2; 16.7%")

# CAP by transect + spider
tempP_centt <- aggregate(cbind(LD1, LD2) ~ Transect, data = cap_temptP_points, FUN = mean)

tempP_segst <- merge(cap_temptP_points, setNames(tempP_centt, c('Transect', 'oLD1', 'oLD2')), by = 'Transect', sort = FALSE)

ggplot(cap_temptP_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = Transect, shape = PlotPos), size = 3, alpha = .6) +
  geom_segment(data = tempP_segst, mapping = aes(x = LD1, y = LD2, xend = oLD1, yend = oLD2, colour = Transect), alpha = .7, size = .25) +
  geom_point(data = tempP_centt, mapping = aes(x = LD1, y = LD2, colour = Transect), size = 5) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_temptP_arrows,
               x = 0, y = 0, alpha = 0.3,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_temptP_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 57.0%",
    y = "CAP Axis 2; 16.7%")


# CAP by plotpos
cap_temppP <- CAPdiscrim(disttempP~PlotPos, data = sttemporalP, axes = 10, m = 0, mmax = 10, add = FALSE, permutations = 9)
cap_temppP <- add.spec.scores(cap_temppP, dtempP, method = "cor.scores", multi = 1, Rscale = F, scaling = "1")
round(cap_temppP$F/sum(cap_temppP$F), digits=3)
barplot(cap_temppP$F/sum(cap_temppP$F))

cap_temppP_points <- bind_cols((as.data.frame(cap_temppP$x)), ftempP) 
glimpse(cap_temppP_points)

cap_temppP_arrows <- as.data.frame(cap_temppP$cproj*5) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")

ggplot(cap_temppP_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = PlotPos), size = 4) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_temppP_arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_temppP_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 80.2%",
    y = "CAP Axis 2; 18.7%")

# CAP by plot + spider
tempP_centp <- aggregate(cbind(LD1, LD2) ~ PlotPos, data = cap_temppP_points, FUN = mean)

tempP_segsp <- merge(cap_temppP_points, setNames(tempP_centp, c('PlotPos', 'oLD1', 'oLD2')), by = 'PlotPos', sort = FALSE)

ggplot(cap_temppP_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = PlotPos), size = 3, alpha = .6) +
  geom_segment(data = tempP_segsp, mapping = aes(x = LD1, y = LD2, xend = oLD1, yend = oLD2, colour = PlotPos), alpha = .9, size = .3) +
  geom_point(data = tempP_centp, mapping = aes(x = LD1, y = LD2, colour = PlotPos), size = 5) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_temppP_arrows,
               x = 0, y = 0, alpha = 0.3,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_temppP_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 80.2%",
    y = "CAP Axis 2; 18.7%")

# CAP by SamplingPeriod
cap_temppsP <- CAPdiscrim(disttempP~`Sampling Period`, data = sttemporalP, axes = 10, m = 0, mmax = 10, add = FALSE, permutations = 9)
cap_temppsP <- add.spec.scores(cap_temppsP, dtempP, method = "cor.scores", multi = 1, Rscale = F, scaling = "1")
round(cap_temppsP$F/sum(cap_temppsP$F), digits=3)
barplot(cap_temppsP$F/sum(cap_temppsP$F))

cap_temppsP_points <- bind_cols((as.data.frame(cap_temppsP$x)), ftempP) 
glimpse(cap_temppsP_points)

cap_temppsP_arrows <- as.data.frame(cap_temppsP$cproj*5) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")
cap_temppsP_arrows

ggplot(cap_temppsP_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = `Sampling Period`), size = 4) +
  scale_colour_manual(values = brewer.pal(n = 6, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_temppsP_arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_temppsP_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 65.2%",
    y = "CAP Axis 2; 22.6%")

# CAP by SamplingPeriod + spider
tempP_centps <- aggregate(cbind(LD1, LD2) ~ `Sampling Period`, data = cap_temppsP_points, FUN = mean)

tempP_segsps <- merge(cap_temppsP_points, setNames(tempP_centps, c('Sampling Period', 'oLD1', 'oLD2')), by = 'Sampling Period', sort = FALSE)

ggplot(cap_temppsP_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = `Sampling Period`, shape = PlotPos), size = 3, alpha = .6) +
  geom_segment(data = tempP_segsps, mapping = aes(x = LD1, y = LD2, xend = oLD1, yend = oLD2, colour = `Sampling Period`), alpha = .9, size = .3) +
  geom_point(data = tempP_centps, mapping = aes(x = LD1, y = LD2, colour = `Sampling Period`), size = 5) +
  scale_colour_manual(values = brewer.pal(n = 6, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_temppsP_arrows,
               x = 0, y = 0, alpha = 0.3,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(2, "mm"))) +
  ggrepel::geom_text_repel(data = cap_temppsP_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 4
  ) +
  labs(
    x = "CAP Axis 1; 65.2%",
    y = "CAP Axis 2; 22.6%")



