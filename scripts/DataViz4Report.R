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
library(patchwork)

#### Colours ####
# No margin
par(mar=c(0,0,1,0))

# Classic palette Spectral, with 11 colors
coul <- brewer.pal(11, "Spectral") 
# Add more colors to this palette :
coul17 <- colorRampPalette(coul)(17)
# Plot it
pie(rep(1, length(coul17)), col = coul17 , main="") 


# Classic palette Spectral, with 11 colors
coul <- brewer.pal(11, "Spectral") 
# Add more colors to this palette :
coul11 <- colorRampPalette(coul)(11)
# Plot it
pie(rep(1, length(coul11)), col = coul11 , main="") 

# Classic palette Spectral, with 11 colors
coul <- brewer.pal(11, "Spectral") 
# Add more colors to this palette :
coul8 <- colorRampPalette(coul)(8)
# Plot it
pie(rep(1, length(coul8)), col = coul8 , main="") 

# Output the palettes for reference
x<-list(coul8, coul11, coul17)
y<-tibble(column1= map_chr(x, str_flatten, " "))
write_csv(y, "colours.csv")


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


#### BW ####

# Landscape data plots

RTHeight <- ggplot(sum) +
  stat_halfeye(aes(y = RTHeight),
    adjust = .5,
    width = .6,
    .width = 0,
    justification = -.3,
    point_colour = NA,
    fill = "#9E0142") +
    geom_point(aes(x = 0, y = RTHeight, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
             ) +
  geom_boxplot(aes(y = RTHeight),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = "Relative height in toposequence (m)",
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

TWI <- ggplot(sum) +
  stat_halfeye(aes(y = TWI),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#D53E4F") +
    geom_point(aes(x = 0, y = TWI, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = TWI),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = "Topographic wetness index",
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


TPI <- ggplot(sum) +
  stat_halfeye(aes(y = TPI),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#F46D43") +
    geom_point(aes(x = 0, y = TPI, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = TPI),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = "Topographic position index",
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Slope <- ggplot(sum) +
  stat_halfeye(aes(y = Slope),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#FDAE61") +
  geom_point(aes(x = 0, y = Slope, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = Slope),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = "Slope",
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

planCurv <- ggplot(sum) +
  stat_halfeye(aes(y = planCurv),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#FEE08B") +
  geom_point(aes(x = 0, y = planCurv, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = planCurv),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = "Plan curvature",
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

proCurv <- ggplot(sum) +
  stat_halfeye(aes(y = proCurv),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#FFFFBF") +
  geom_point(aes(x = 0, y = proCurv, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = proCurv),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = "Profile curvature",
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

NDVI <- ggplot(all) +
  stat_halfeye(aes(y = NDVI),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#E6F598") +
  geom_point(aes(x = 0, y = NDVI, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = NDVI),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = "Normalised difference vegetation index (NDVI)",
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Wet <- ggplot(all) +
  stat_halfeye(aes(y = Wet),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#ABDDA4") +
  geom_point(aes(x = 0, y = Wet, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = Wet),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = "Soil moisture by synthetic aperture radar (Sentinel)",
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Moisture <- ggplot(all) +
  stat_halfeye(aes(y = Moisture),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#66C2A5") +
  geom_point(aes(x = 0, y = Moisture, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = Moisture),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Soil moisture (g"~g^-1~" dry weight)"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

WHC <- ggplot(sum) +
  stat_halfeye(aes(y = WHC),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#3288BD") +
  geom_point(aes(x = 0, y = WHC, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = WHC),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Water holding capacity (g"~g^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

BD0_30 <- ggplot(sum) +
  stat_halfeye(aes(y = BD0_30),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#5E4FA2") +
  geom_point(aes(x = 0, y = BD0_30, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = BD0_30),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Bulk density (g"~cm^-3~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


RTHeight + TWI + TPI + Slope + planCurv + proCurv + NDVI + 
  Wet + Moisture + WHC + BD0_30 + guide_area() +
  plot_layout(ncol = 6, guides = 'collect') +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag.position = c(1, 1),
        plot.tag = element_text(size = 16, hjust = 4, vjust = 2))

#y = expression ("Bulk density g"~cm^-3)

# Chem data
pHc <- ggplot(all) +
  stat_halfeye(aes(y = pHc),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#9E0142") +
  geom_point(aes(x = 0, y = pHc, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = pHc),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression (~pH[CaCl[2]]),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


EC <- ggplot(all) +
  stat_halfeye(aes(y = EC),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#C0247A") +
  geom_point(aes(x = 0, y = EC, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = EC),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Electrical conductivity (dS "~m^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

CEC <- ggplot(sum) +
  stat_halfeye(aes(y = CEC),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#DC494C") +
  geom_point(aes(x = 0, y = CEC, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = CEC),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Cation exchange capacity ("~cmol^+~" "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


PC1 <- ggplot(metals_plot) +
  stat_halfeye(aes(y = Comp.1),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#F06744") +
  geom_point(aes(x = 0, y = Comp.1, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = Comp.1),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Total elements principal component 1, 58.3% of variance"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


PC2 <- ggplot(metals_plot) +
  stat_halfeye(aes(y = Comp.2),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#F88D51") +
  geom_point(aes(x = 0, y = Comp.2, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = Comp.2),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Total elements principal component 2, 13.9% of variance"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


P <- ggplot(sum) +
  stat_halfeye(aes(y = P),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#FDB466") +
  geom_point(aes(x = 0, y = P, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = P),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Total phosphorus (mg "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


K <- ggplot(sum) +
  stat_halfeye(aes(y = K),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#FDD380") +
  geom_point(aes(x = 0, y = K, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = K),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Total potassium (mg "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


S <- ggplot(sum) +
  stat_halfeye(aes(y = S),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#FEEB9E") +
  geom_point(aes(x = 0, y = S, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = S),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Total sulphur (mg "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


TotOC <- sum %>% drop_na(TotOC_mean) %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = TotOC_mean),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#FFFFBF") +
  geom_point(aes(x = 0, y = TotOC_mean, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = TotOC_mean),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Total organic carbon (g "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


TotN <- sum %>% drop_na(TotN_mean) %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = TotN_mean),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#EFF8A6") +
  geom_point(aes(x = 0, y = TotN_mean, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = TotN_mean),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Total nitrogen (g "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



CN <- sum %>% drop_na(CN_mean) %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = CN_mean),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#D7EF9B") +
  geom_point(aes(x = 0, y = CN_mean, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = CN_mean),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("C:N ratio"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



d13C <- sum %>% drop_na(d13C_mean) %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = d13C_mean),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#B2E0A2") +
  geom_point(aes(x = 0, y = d13C_mean, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = d13C_mean),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression (paste(delta^{13}, "C (\u2030)")),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


d15N <- sum %>% drop_na(d15N_mean) %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = d15N_mean),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#88CFA4") +
  geom_point(aes(x = 0, y = d15N_mean, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = d15N_mean),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression (paste(delta^{15}, "N (\u2030)")),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


POC <- sum %>% drop_na(POC_mean) %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = POC_mean),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#5FBAA8") +
  geom_point(aes(x = 0, y = POC_mean, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = POC_mean),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Particulate organic carbon (g "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


HOC <- sum %>% drop_na(HOC_mean) %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = HOC_mean),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#3F96B7") +
  geom_point(aes(x = 0, y = HOC_mean, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = HOC_mean),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Humus organic carbon (g "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


ROC <- sum %>% drop_na(ROC_mean) %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = ROC_mean),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#4272B2") +
  geom_point(aes(x = 0, y = ROC_mean, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = ROC_mean),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Resistant organic carbon (g "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


Vuln <- sum %>% drop_na(Vuln_mean) %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = Vuln_mean),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#5E4FA2") +
  geom_point(aes(x = 0, y = Vuln_mean, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = Vuln_mean),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Organic carbon vulnerability"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


pHc + EC + CEC + PC1 + PC2 + P + K + S + TotOC +
  TotN + CN + d13C + d15N + POC + HOC + ROC + Vuln + guide_area() +
  plot_layout(ncol = 6, guides = 'collect') +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag.position = c(1, 1),
        plot.tag = element_text(size = 16, hjust = 2, vjust = 2))


### Dynamic 
NO3 <- all %>% drop_na(NO3) %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = NO3),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#9E0142") +
  geom_point(aes(x = 0, y = NO3, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = NO3),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Extractable "~NO[3]^{"-"}~"-N (mg "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


NH4 <- all %>% drop_na(NH4) %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = NH4),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#E25249") +
  geom_point(aes(x = 0, y = NH4, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = NH4),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Extractable "~NH[4]^{"+"}~"-N (mg "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


FAA <- all %>% drop_na(FAA) %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = FAA),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#FBA45C") +
  geom_point(aes(x = 0, y = FAA, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = FAA),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Extractable free amino acid-N (mg "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


DON <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = DON),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#FEE899") +
  geom_point(aes(x = 0, y = DON, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = DON),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Dissolved organic N (mg "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


DOC <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = DOC),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#EDF7A3") +
  geom_point(aes(x = 0, y = DOC, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = DOC),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Dissolved organic C (mg "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


MBC <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = MBC),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#A1D9A4") +
  geom_point(aes(x = 0, y = MBC, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = MBC),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Microbial biomass C (mg "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


MBN <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = MBN),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#48A0B2") +
  geom_point(aes(x = 0, y = MBN, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = MBN),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Microbial biomass N (mg "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


AvailP <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = AvailP),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#5E4FA2") +
  geom_point(aes(x = 0, y = AvailP, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = AvailP),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Olsen-extractable P (mg "~kg^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


NO3 + NH4 + FAA + DON + DOC + MBC + MBN + AvailP +
  plot_layout(ncol = 4, guides = 'collect') +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 16, hjust = -12, vjust = 2))


### Microbial
Proteolysis <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = Proteolysis),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#9E0142") +
  geom_point(aes(x = 0, y = Proteolysis, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = Proteolysis),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Proteolysis rate (mg AA-N"~kg^-1~h^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


AAMin_k1 <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = AAMin_k1),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#D53E4F") +
  geom_point(aes(x = 0, y = AAMin_k1, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = AAMin_k1),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Rate of initial AA mineralisation ("~h^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


MicY <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = MicY),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#F46D43") +
  geom_point(aes(x = 0, y = MicY, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = MicY),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Microbial yield"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


TotalPLFA <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = TotalPLFA),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#FDAE61") +
  geom_point(aes(x = 0, y = TotalPLFA, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = TotalPLFA),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Total PLFA (nmol "~g^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


Bac <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = Bac),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#FEE08B") +
  geom_point(aes(x = 0, y = Bac, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = Bac),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Bacterial PLFA (nmol "~g^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


Fun <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = Fun),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#FFFFBF") +
  geom_point(aes(x = 0, y = Fun, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = Fun),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Fungal PLFA (nmol "~g^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


Gpos <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = Gpos),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#E6F598") +
  geom_point(aes(x = 0, y = Gpos, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = Gpos),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("G+ bacterial PLFA (nmol "~g^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


Gneg <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = Gneg),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#ABDDA4") +
  geom_point(aes(x = 0, y = Gneg, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = Gneg),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("G- bacterial PLFA (nmol "~g^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


Act <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = Act),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#66C2A5") +
  geom_point(aes(x = 0, y = Act, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = Act),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Actinomycete PLFA (nmol "~g^-1~")"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


F_B <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = F_B),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#3288BD") +
  geom_point(aes(x = 0, y = F_B, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = F_B),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Fungal:Bacterial ratio"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


Gp_Gn <- all %>% drop_na() %>% # Neat little hack to drop NA samples
  ggplot() + # Also need to drop the df call here
  stat_halfeye(aes(y = Gp_Gn),
               adjust = .5,
               width = .6,
               .width = 0,
               justification = -.3,
               point_colour = NA,
               fill = "#5E4FA2") +
  geom_point(aes(x = 0, y = Gp_Gn, colour = Transect),
             shape = 21,
             stroke = 1,
             size = 3,
             position = position_jitter(
               seed = 1,
               width = 0.1
             )
  ) +
  geom_boxplot(aes(y = Gp_Gn),
               alpha = 0,
               width = .25,
               outlier.shape = NA
  ) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  labs(y = expression ("Gram+:Gram- ratio"),
       colour = "Toposequence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


Proteolysis + AAMin_k1 + MicY + TotalPLFA + Bac + Fun + 
  Gpos + Gneg + Act + F_B + Gp_Gn + guide_area() +
  plot_layout(ncol = 6, guides = 'collect') +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag.position = c(1, 1),
        plot.tag = element_text(size = 16, hjust = 4, vjust = 2))

#### xy plots ####

# Add plot position
sum %<>% group_by(Transect) %>% 
  mutate(PlotPos = dense_rank(desc(RTHeight))) %>%
  ungroup() %>%
  relocate(PlotPos, .after = Plot) %>% 
  mutate(across(c(CombID, UniqueID, PrelimID, Transect, Plot, Inun, PlotPos), as.factor)) 
str(sum)

# isotopes
#CN
cn_c <- ggplot(sum) + 
  geom_point(aes(x=CN_mean, y=d13C_mean, colour = PlotPos), size = 3) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(
    x = "C:N ratio",
    y = expression (paste(delta^{13}, "C (\u2030)")),
    colour = "Plot position")

cn_n <- ggplot(sum) + 
  geom_point(aes(x=CN_mean, y=d15N_mean, colour = PlotPos), size = 3) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(
    x = "C:N ratio",
    y = expression (paste(delta^{15}, "N (\u2030)")),
    colour = "Plot position")
#vuln
vuln_c <- ggplot(sum) + 
  geom_point(aes(x=Vuln_mean, y=d13C_mean, colour = PlotPos), size = 3) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(
    x = "SOC vulnerability",
    y = expression (paste(delta^{13}, "C (\u2030)")),
    colour = "Plot position")

vuln_n <- ggplot(sum) + 
  geom_point(aes(x=Vuln_mean, y=d15N_mean, colour = PlotPos), size = 3) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(
    x = "SOC vulnerability",
    y = expression (paste(delta^{15}, "N (\u2030)")),
    colour = "Plot position")
#iso only
iso <- ggplot(sum) + 
  geom_point(aes(x=d13C_mean, y=d15N_mean, colour = PlotPos), size = 3) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(
    x = expression (paste(delta^{13}, "C (\u2030)")),
    y = expression (paste(delta^{15}, "N (\u2030)")),
    colour = "Plot position")

cn_c + cn_n + iso + vuln_c + vuln_n + guide_area() +
  plot_layout(ncol = 3, guides = 'collect') +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 16, hjust = -5, vjust = 1))

#### local scale ####
#### biogeochem ####

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

bgc_mean <- t1_summary %>% 
  select(UniqueID, Transect, Plot, PlotPos, Easting, Northing, Height, RHeight, RTHeight, Inun,
         Clay, CEC, WHC, BD0_30, NDVI_mean, Wet_mean, Moisture_mean, pHc_mean, EC_mean, AvailP_mean, CN_mean, Vuln_mean,
         d13C_mean, d15N_mean, DOC_mean, NO3_mean, NH4_mean, FAA_mean, Proteolysis_mean,
         AAMin_k1_mean, DON_mean, MBC_mean, MBN_mean, MicY_mean)

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
  geom_point(aes(x=Axis.1, y=Axis.2, colour = PlotPos), size = 6) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
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
saveRDS(cap_bgct, file = "data/processed/CAP_bgct.rds")


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

cap_bgct_fig <- ggplot(cap_bgct_points) + 
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
    y = "CAP Axis 2; 23.0%",
    colour = "Toposequence",
    shape = "Plot position")


# CAP by plotpos
stbgc_mean <- as.data.frame(stbgc_mean) 
cap_bgcp <- CAPdiscrim(distbgc~PlotPos, data = stbgc_mean, axes = 10, m = 3, mmax = 10, add = FALSE, permutations = 999)
cap_bgcp <- add.spec.scores(cap_bgcp, dbgc, method = "cor.scores", multi = 1, Rscale = F, scaling = "1")
saveRDS(cap_bgcp, file = "data/processed/CAP_bgcp.rds")
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
bgc_centp <- aggregate(cbind(LD1, LD2) ~ PlotPos, data = cap_bgcp_points, FUN = mean)

bgc_segsp <- merge(cap_bgcp_points, setNames(bgc_centp, c('PlotPos', 'oLD1', 'oLD2')), by = 'PlotPos', sort = FALSE)

cap_bgcpfig <- ggplot(cap_bgcp_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = PlotPos), size = 3, alpha = .6) +
  geom_segment(data = bgc_segsp, mapping = aes(x = LD1, y = LD2, xend = oLD1, yend = oLD2, colour = PlotPos), alpha = .9, size = .3) +
  geom_point(data = bgc_centp, mapping = aes(x = LD1, y = LD2, colour = PlotPos), size = 5) +
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
    x = "CAP Axis 1; 76.3%",
    y = "CAP Axis 2; 23.7%",
    colour = "Plot position")

cap_bgct_fig + cap_bgcpfig +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 16, hjust = -5, vjust = 1))


#### temporal ####
OL_cor <- read_csv("data/processed/ChemAll_adm_OLrem.csv")
OL_cor <- OL_cor %>% 
  group_by(Transect) %>% 
  mutate(PlotPos = dense_rank(desc(RTHeight))) %>%
  ungroup() %>% 
  relocate(PlotPos, .after = Plot) %>% 
  mutate(across(c(CombID, UniqueID, PrelimID, Transect, Plot, Inun, PlotPos), as.factor)) %>% 
  mutate(Date = dmy(Date)) 
str(OL_cor)

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

temporalP <- OLP_cor %>% 
  select(UniqueID, Date, `Sampling Period`, Transect, Plot, PlotPos, Easting, Northing, Height, RHeight, RTHeight, Inun,
         NDVI, VH,	VV,	Wet, Moisture,	pHc,	EC, AvailP, 
         DOC,	DTN,	NO3,	NH4,	FAA, Proteolysis,	AAMin_k1,	DON,	MBC,	
         MBN,	MicY,	MicCN, TotalPLFA, F_B, Gp_Gn, Act_Gp)

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
cap_temptP <- CAPdiscrim(disttempP~Transect, data = sttemporalP, axes = 10, m = 0, mmax = 10, add = FALSE, permutations = 99)
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
cap_temppsP <- CAPdiscrim(disttempP~`Sampling Period`, data = sttemporalP, axes = 10, m = 0, mmax = 10, add = FALSE, permutations = 999)
cap_temppsP <- add.spec.scores(cap_temppsP, dtempP, method = "cor.scores", multi = 1, Rscale = F, scaling = "1")
saveRDS(cap_temppsP, file = "outputs/cap_temppsP.rds")
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
  geom_point(aes(x=LD1, y=LD2, colour = `Sampling Period`, shape = PlotPos), size = 2.5, alpha = .4) +
  geom_segment(data = tempP_segsps, mapping = aes(x = LD1, y = LD2, xend = oLD1, yend = oLD2, colour = `Sampling Period`), alpha = .9, size = .3) +
  geom_point(data = tempP_centps, mapping = aes(x = LD1, y = LD2, colour = `Sampling Period`), size = 8) +
  scale_colour_manual(values = brewer.pal(n = 5, name = "Set1")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  geom_segment(data = cap_temppsP_arrows,
               x = 0, y = 0, alpha = 0.6,
               mapping = aes(xend = LD1, yend = LD2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = cap_temppsP_arrows, aes(x=LD1, y=LD2, label = variable), 
                           # colour = "#72177a", 
                           size = 5
  ) +
  labs(
    x = "CAP Axis 1; 65.2%",
    y = "CAP Axis 2; 22.6%", 
    shape = "Plot position")

#### temporal trends ####
#This needs to be a multi-panel figure(s) y = var, x = date, colour = plot position, thick lines and points = mean, hairlines = toposequences
# 1) TICK - make a df with only vars of interest
# 2) TICK - Make summary df with means by landscape position
# 3) TICK - Plot individuals with feint lines, colours by landscape position
# 4) TICK - Overlay points and thicker lines, colours by landscape position

seasonal <- temporalP %>% 
  select(-c(VH, VV, pHc, EC, DTN, MBC)) %>% 
  unite("Tr_PP", Transect:PlotPos, remove = FALSE)

seasonal_vars <- c("Date", "Moisture", "FAA", "NO3", "DON", "NH4", "AvailP", "DOC", "NDVI", "Wet", "Proteolysis", "AAMin_k1", "Gp_Gn", "F_B", "TotalPLFA", "MBN", "MicCN", "Act_Gp", "MicY")

seasonal_sum <- seasonal %>% 
  group_by(`Sampling Period`, PlotPos) %>%
  summarise(across(all_of(seasonal_vars),
                   list(mean = ~ mean(.x, na.rm = TRUE)))) %>% 
  ungroup() 
  
prot <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = Proteolysis, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = Proteolysis_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = Proteolysis_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("Proteolysis rate"),
    colour = "Plot position")

moist <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = Moisture, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = Moisture_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = Moisture_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("MC (g "~g^-1~")"),
    colour = "Plot position")

faa <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = FAA, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = FAA_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = FAA_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("FAA-N (mg N "~kg^-1~")"),
    colour = "Plot position")

no3 <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = NO3, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = NO3_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = NO3_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression (~NO[3]^{"-"}~"-N (mg "~kg^-1~")"),
    colour = "Plot position")

nh4 <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = NH4, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = NH4_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = NH4_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression (~NH[4]^{"+"}~"-N (mg "~kg^-1~")"),
    colour = "Plot position")

don <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = DON, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = DON_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = DON_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("DON (mg "~kg^-1~")"),
    colour = "Plot position")

doc <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = DOC, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = DOC_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = DOC_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("DOC (mg "~kg^-1~")"),
    colour = "Plot position")

availp <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = AvailP, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = AvailP_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = AvailP_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("Available P (mg "~kg^-1~")"),
    colour = "Plot position")

aak1 <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = AAMin_k1, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = AAMin_k1_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = AAMin_k1_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("AA min ("~h^-1~")"),
    colour = "Plot position")

cue <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = MicY, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = MicY_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = MicY_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("Amino acid CUE"),
    colour = "Plot position")

gpgn <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = Gp_Gn, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = Gp_Gn_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = Gp_Gn_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("G+ : G- ratio"),
    colour = "Plot position")

actgp <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = Act_Gp, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = Act_Gp_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = Act_Gp_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("Actinomycete : G+ ratio"),
    colour = "Plot position")

fb <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = F_B, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = F_B_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = F_B_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("Fungal : Bacterial ratio"),
    colour = "Plot position")

mbn <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = MBN, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = MBN_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = MBN_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("MBN (mg "~kg^-1~")"),
    colour = "Plot position")

miccn <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = MicCN, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = MicCN_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = MicCN_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("Microbial biomass C:N ratio"),
    colour = "Plot position")

totp <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = TotalPLFA, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = TotalPLFA_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = TotalPLFA_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("Total PLFA (nmol "~g^-1~")"),
    colour = "Plot position")

ndvi <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = NDVI, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = NDVI_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = NDVI_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("NDVI"),
    colour = "Plot position")

wet <- ggplot() + 
  geom_line(data = seasonal, aes(group = Tr_PP, x = Date, y = Wet, colour = PlotPos), size = 0.05) +
  geom_line(data = seasonal_sum, aes(x = Date_mean, y = Wet_mean, colour = PlotPos), size = 1) +
  geom_point(data = seasonal_sum, aes(x = Date_mean, y = Wet_mean, colour = PlotPos), size = 2) +
  scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y") +
  labs(
    x = "",
    y = expression ("Wetness index"),
    colour = "Plot position")

no3 + nh4 + faa + don + doc + availp + prot + aak1 + cue + moist +
  plot_annotation(tag_levels = 'a') + 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 16, hjust = 0, vjust = 1)) +
  plot_layout(ncol = 2, guides = 'collect') & theme(legend.position = 'bottom')

