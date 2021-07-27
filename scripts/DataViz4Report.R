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
