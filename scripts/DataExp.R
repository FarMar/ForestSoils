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
         Clay, CEC, WHC, BD0_30, NDVI_mean, Wet_mean, Moisture_mean, pHc_mean, EC_mean, AvailP_mean, CN_mean, Vuln_mean,
         d13C_mean, d15N_mean, DOC_mean, NO3_mean, NH4_mean, FAA_mean, Proteolysis_mean,
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
stmir <- as.data.frame(stmir) 
cap_mirt <- CAPdiscrim(distmir~Transect, data = stmir, axes = 10, m = 0, mmax = 10, add = FALSE, permutations = 999)
round(cap_mirt$F/sum(cap_mirt$F), digits=3)
barplot(cap_mirt$F/sum(cap_mirt$F))

cap_mirt_points <- bind_cols((as.data.frame(cap_mirt$x)), fmir) 
glimpse(cap_mirt_points)

ggplot(cap_mirt_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = Transect), size = 4) +
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
  geom_point(aes(x=LD1, y=LD2, colour = Transect, shape = `Sampling Period`), size = 3, alpha = .7) +
  geom_segment(data = mir_segs, mapping = aes(x = LD1, y = LD2, xend = oLD1, yend = oLD2, colour = Transect), alpha = .7, size = .25) +
  geom_point(data = mir_cent, mapping = aes(x = LD1, y = LD2, colour = Transect), size = 5) +
  scale_colour_manual(values = brewer.pal(n = 10, name = "Spectral")) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(
    x = "CAP Axis 1; 41.2%",
    y = "CAP Axis 2; 35.3%")



##BGC
#pre-prep - PCA of total emlements to reduce dimenstions
tot_elms <- t1_summary %>% 
  select(46:65) %>% 
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
bgc_cor <- select(bgc_mean, 10:35)
chart.Correlation(bgc_cor, histogram=TRUE, pch=19)


tbgc_mean <- bgc_mean %>% 
  mutate(MBN_mean = log1p(MBN_mean),
         NH4_mean = log1p(NH4_mean),
         AvailP_mean = log1p(AvailP_mean),
         EC_mean = log1p(EC_mean),
         pHc_mean = log1p(pHc_mean),
         BD0_30 = log1p(BD0_30))

stbgc_mean <- tbgc_mean %>% 
  mutate(across(c(10:35), ~z.fn(.)))

fbgc <- stbgc_mean %>% 
  select(1:9)
dbgc <- stbgc_mean %>% 
  select(10:35)

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
  geom_point(aes(x=Axis.1, y=Axis.2, colour = Transect, shape = Plot), size = 6) +
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
    x = "PCoA Axis 1; 59.1%",
    y = "PCoA Axis 2; 11.7%")

# Permanova
set.seed(1983)
perm_bgc <- adonis2(distbgc~Transect+Plot, data = stbgc_mean, permutations = 9999, method = "euclidean")
perm_bgc #strong impact of transect and plot
permpt_bgc <- pairwise.perm.manova(distbgc, stbgc_mean$Transect, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpt_bgc #.089 is lowest possible  - several pairwise comps have this
permpp_bgc <- pairwise.perm.manova(distbgc, stbgc_mean$Plot, nperm = 9999, progress = TRUE, p.method = "fdr", F = TRUE, R2 = TRUE)
permpp_bgc #sniff of significance for last sampling vs 1st three samplings


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
  geom_point(aes(x=LD1, y=LD2, colour = Transect), size = 4) +
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
  geom_point(aes(x=LD1, y=LD2, colour = Transect), size = 3, alpha = .6) +
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


# CAP by plot
stbgc_mean <- as.data.frame(stbgc_mean) 
cap_bgcp <- CAPdiscrim(distbgc~Plot, data = stbgc_mean, axes = 10, m = 6, mmax = 10, add = FALSE, permutations = 999)
cap_bgcp <- add.spec.scores(cap_bgcp, dbgc, method = "cor.scores", multi = 1, Rscale = F, scaling = "1")
round(cap_bgcp$F/sum(cap_bgcp$F), digits=3)
barplot(cap_bgcp$F/sum(cap_bgcp$F))

cap_bgcp_points <- bind_cols((as.data.frame(cap_bgcp$x)), fbgc) 
glimpse(cap_bgcp_points)

cap_bgcp_arrows <- as.data.frame(cap_bgcp$cproj*5) %>% #Pulls object from list, scales arbitrarily and makes a new df
  rownames_to_column("variable")

ggplot(cap_bgcp_points) + 
  geom_point(aes(x=LD1, y=LD2, colour = Plot), size = 4) +
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
    x = "CAP Axis 1; 85.4%",
    y = "CAP Axis 2; 14.0%")

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

## Temporal
