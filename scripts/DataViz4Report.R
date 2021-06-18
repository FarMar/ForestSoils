#####################################################################################################

#### Forest soils dataviz script                                                  ###################

#### mark.farrell@csiro.au        +61 8 8303 8664         31/05/2021 ################################

#####################################################################################################


#### Set working directory ####
setwd("/Users/markfarrell/OneDrive - CSIRO/Data/ForestSoils")


#### Packages ####
install.packages("ggtern")
install.packages("ggdist")


library(tidyverse)
library(janitor)
library(PerformanceAnalytics)
library(corrplot)
library(RColorBrewer)
library(plotrix)
library(ggpmisc)
#library(ggtern)
library(ggbluebadge)
library(ggdist)
library(magrittr)
library(lubridate)
library(vegan)
library(ape)
library(RVAideMemoire)
library(BiodiversityR)

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
# This is best run standalone as {ggtern} masks a lot of ggplot 
ggtern(data=sum, aes(Sand,Clay,Silt, color = Transect)) + 
  geom_point(size = 4) +
  theme_rgbw() +
  theme_hidetitles() +
  theme(text = element_text(size=20)) +
  theme(legend.key=element_blank())


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
mir_meta <- all %>%
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
saveRDS(cap_mirt, file = "outputs/MIRCAP.rds")
readRDS("outputs/MIRCAP.rds")
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
mir_segs <- merge(cap_mirt_points, setNames(mir_cent, c('Transect', 'oLD1', 'oLD2')), by = 'Transect', sort = FALSE)

ggplot(cap_mirt_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = Transect, shape = PlotPos), size = 3, alpha = .7) +
  geom_segment(data = mir_segs, mapping = aes(x = LD1, y = LD2, xend = oLD1, yend = oLD2, colour = Transect), alpha = .5, size = .25) +
  geom_point(data = mir_cent, mapping = aes(x = LD1, y = LD2, colour = Transect), size = 5, alpha = 1.0) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(
    x = "CAP Axis 1; 41.2%",
    y = "CAP Axis 2; 35.3%")


#### Metals PCA ####
metals <- sum %>% 
  select(c(1:11, 45:65)) %>% 
  select(-c(As, Cd, Mo, Sb, Se)) %>% 
  group_by(Transect) %>%
  mutate(PlotPos = dense_rank(desc(RTHeight))) %>%
  ungroup() %>% 
  relocate(PlotPos, .after = Plot) %>% 
  mutate(across(c(CombID, UniqueID, PrelimID, Transect, Plot, PlotPos), as.factor))

metals %<>% 
  mutate(P = log1p(P),
         Na = log1p(Na),
         Mg = log1p(Mg),
         K = log1p(K),
         Co = log1p(Co),
         Ca = log1p(Ca))
chart.Correlation(metals[13:28])

pca_metals <- princomp(metals[13:28], cor = TRUE, scores = TRUE)
biplot(pca_metals, choices = c(1,2))
summary(pca_metals) #PC1 = 58.3%, PC2 = 13.9%
scores_metals <- as.data.frame(pca_metals[["scores"]]) %>% 
  select(1:2) 

metals_plot <- bind_cols(metals, scores_metals)

metals_cent <- aggregate(cbind(Comp.1, Comp.2) ~ Transect, data = metals_plot, FUN = mean)
metals_segs <- merge(metals_plot, setNames(metals_cent, c('Transect', 'PC1', 'PC2')), by = 'Transect', sort = FALSE)

ggplot(metals_plot) + 
  geom_point(aes(x=Comp.1, y=Comp.2, colour = Transect, shape = PlotPos), size = 3, alpha = .7) +
  geom_segment(data = metals_segs, mapping = aes(x = Comp.1, y = Comp.2, xend = PC1, yend = PC2, colour = Transect), alpha = .5, size = .25) +
  geom_point(data = metals_cent, mapping = aes(x = Comp.1, y = Comp.2, colour = Transect), size = 5, alpha = 1.0) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(
    x = "PCA Axis 1; 58.3%",
    y = "PCA Axis 2; 13.9%")


#### Chem CAP ####



#### BW ####


