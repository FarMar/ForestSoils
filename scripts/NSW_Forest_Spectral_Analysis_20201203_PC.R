##############################################################################################################################################################################
########### Development Infrared Spectroscopic Model for NSW Forest soils contents ############################################################################################
###############################################################################################################################################################################

##################################################################################################################################
## Authors: 
## Senani Karunaratne | senanni.karunaratne@csiro.au
## Mark Farrell | mark.farrell@csiro.au 
## Start date: 20201202
## End date: 20201203
## TTD
# Check the baseline correction function
# additional spectral pre-processing functions
## ###################################################################################################################################


### Required libraries
library(plyr)
library(dplyr)
library(pls)
library(clhs)
library(prospectr)
library(resemble)
#install.packages("remotes")
remotes::install_github("philipp-baumann/simplerspec")
library(simplerspec)
library(baseline)
library(hydroGOF)

#######################################################################################################################################
### Preparation of data sets

### Set working directory
setwd('Z:\\reference\\Team\\Senani\\NSW_Forest')
getwd()
dir()


### Read analytical data
data <- read.csv('MasterFieldDataFC_NSW - Data.csv')
head(data)
dim(data)
summary(data$Clay)

### Read spectral data
spec <- read.csv('MasterFieldDataFC_NSW - MIR_raw.csv')
head(spec)
dim(spec)
names(spec)


### Select only spec and ID
spec_s <- spec %>%
   dplyr::select(2, 17:1987)
head(spec_s)
dim(spec_s)

#remove the 'X' from headers (MIR wave numbers)
colnames(spec_s) <- sub("X","",colnames(spec_s))  


### Quick plot to see the spectra
waves <- seq(7999.28, 401.1211, by= -3.8569)
colnames(spec_s[,2:1972]) <- waves

matplot(x = waves, y = t(spec_s[,2:1972]), ylim=c(0,3.5), type = "l", lty = 1,
        main = "Raw spectra", xlab = "Wavelength (nm)",
        ylab = "Absorbance", col = rep(palette(), each = 3))

### Convert absorbance to reflectance 
##spec_a <- log(1/spec_s[,2:1972])

spec_a <- spec_s[,2:1972]


matplot(x = waves, y = t(spec_a), ylim=c(0,3.5), type = "l", lty = 1,
        main = "Absorbance", xlab = "Wavelength (nm)",
        ylab = "Absorbance", col = rep(palette(), each = 3))

### Trim data between 6000 and 600
spec_a_d_6000_600 <- spec_a[,519:1919]

waves_s <- seq(6001.388, 601.6816, by= -3.8569)
colnames(spec_a_d_6000_600) <- waves_s

matplot(x = waves_s, y = t(spec_a_d_6000_600), ylim=c(0,3), type = "l", lty = 1,
        main = "Absorbance - 600 to 6000", xlab = "Wavelength (nm)",
        ylab = "Absorbance", col = rep(palette(), each = 3))

### Baseline correction
### Senani: Still need to check this function 
spec_a_d_6000_600 <- as.matrix(spec_a_d_6000_600)
spec_a_bc <- baseline(spec_a_d_6000_600)
str(spec_a_bc)
spec_a_bc_d <- getCorrected(spec_a_bc)
spec_a_bc_d <- setNames(spec_a_bc_d,names(spec_a_d_6000_600))

spec_a_bc_d <- as.data.frame(spec_a_bc_d)

matplot(x = waves_s, y = t(spec_a_bc_d), ylim=c(-0.1,1.5), type = "l", lty = 1,
        main = "Absorbance", xlab = "Wavelength (nm)",
        ylab = "Absorbance", col = rep(palette(), each = 3))


### Merge the Spectral ID and merge 
UniqueID <- spec_s$UniqueID

spec_a_bc_d_6000_600_f <- cbind(UniqueID, spec_a_bc_d)
#names(spec_a_bc_d_6000_600_f)


### Select soil property and merge with spectral dataset
head(data)

## Example: clay 
analytical <- dplyr::select(data, UniqueID, Clay)
head(analytical)

## Merge data set
final_d <- dplyr::left_join(analytical, spec_a_bc_d_6000_600_f, by='UniqueID')
dim(final_d)

######################################################################################################################################################
### Divide the data set as model calibration and validation
## 80 % cal and 20 % val
## Based on KS

set.seed(1981)
ken_mahal <- kenStone(X = final_d[,3:1403], k = 160, metric = "mahal", pc= 10)
ken_mahal$model

# The pc components in the output list stores the pc scores
plot(ken_mahal$pc[,1], 
ken_mahal$pc[,2], 
col = rgb(0, 0, 0, 0.3), 
pch = 19, 
xlab = "PC1",
ylab = "PC2",
main = "Kennard-Stone") 
grid()

# This is the selected points in the pc space
points(ken_mahal$pc[ken_mahal$model, 1],
ken_mahal$pc[ken_mahal$model,2],
pch = 19, col = "red")

## Index for selected samples
index <- ken_mahal$model

## Cal data set
data_cal <- final_d[index,]
dim(data_cal)

## Val data set
data_val <- final_d[-index,]
dim(data_val)
#############################################################################################################################################
### Fitting a PLSR model
data_cal$UniqueID <- NULL
ncomp <- 10
fit_plsr <- plsr(Clay ~ ., data=data_cal, ncomp, validation = "CV")
plot(fit_plsr)
abline(0,1)

### See the components returned lowest RMSE and select the best one
plot(fit_plsr, "val") 
ncomp.onesigma <- selectNcomp(fit_plsr, method = "onesigma", plot = TRUE,ylim = c(0, .5))
ncomp.onesigma


### Model validation
data_val$UniqueID <- NULL
val_pred <- predict(fit_plsr, newdata = data_val, ncomp=ncomp.onesigma)
plot(data_val$Clay, val_pred, xlim=c(0,60), ylim=c(0,60), xlab='Measured', ylab='Predicted')
abline(0,1, col='blue')


### Validation statistics 
gof(as.data.frame(val_pred), as.data.frame(data_val$Clay))


