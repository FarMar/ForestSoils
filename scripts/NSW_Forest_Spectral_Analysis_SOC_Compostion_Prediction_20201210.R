##############################################################################################################################################################################
########### Prediction of Soil Carbon Fractions ##############################################################################################################################
###############################################################################################################################################################################

##################################################################################################################################
## Authors: 
## Senani Karunaratne | senanni.karunaratne@csiro.au
## Mark Farrell | mark.farrell@csiro.au 
## Jeff Baldock | jeff.baldock@csiro.au
## Start date: 20201210
## End date: 20201210
## ################################################################################################################################


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
library(hyperSpec)

#######################################################################################################################################
### Preparation of data sets

### Set working directory
#setwd('V:\\_Projects\\2019 Projects\\2019_Farrell_NSW_Forestry\\MIR\\Analysis')
getwd()
dir()


### Read analytical data
data <- read.csv("data/working/MasterFieldDataFC_NSW - Data.csv")
head(data)
dim(data)
summary(data$Clay)

### Read spectral data
spec <- read.csv("data/working/MasterFieldDataFC_NSW - MIR_raw.csv")
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
head(spec_s)

### Quick plot to see the spectra
waves <- seq(7999.28, 401.1211, by= -3.8569)
colnames(spec_s[,2:1972]) <- waves

matplot(x = waves, y = t(spec_s[,2:1972]), ylim=c(0,3.5), type = "l", lty = 1,
        main = "Raw spectra", xlab = "Wavelength (nm)",
        ylab = "Absorbance", col = rep(palette(), each = 3))

### Convert absorbance to reflectance - Not needed for Thermo Nicolet files
##spec_a <- log(1/spec_s[,2:1972])

spec_a <- spec_s[,2:1972] # Just so rest of code runs easily


#matplot(x = waves, y = t(spec_a), ylim=c(0,3.5), type = "l", lty = 1,
#        main = "Absorbance", xlab = "Wavelength (nm)",
#        ylab = "Absorbance", col = rep(palette(), each = 3))


#colnames(spec_a) <- waves_s

#matplot(x = waves_s, y = t(spec_a_d_6000_600), ylim=c(0,3), type = "l", lty = 1,
#        main = "Absorbance - 600 to 6000", xlab = "Wavelength (nm)",
#        ylab = "Absorbance", col = rep(palette(), each = 3))




#######################################################################################################################################
###  Re-sample / Interpolation of the data set ################################################################################################################
#######################################################################################################################################
## Senani noted: Here used the entire data matrix rather trimmed it to 6100 to 500 

mirinterp <- spec_a  ##[,2:1402] #assigns the file to be processed to the object name "mirinterp".  This object will be processed and interpolated.
sample <- spec


#mirinterp[,1] <- NULL  #remove first column with sample names
#mirinterp[,2] <- NULL  #remove second column with Project name

colnames(mirinterp) <- sub("X","",colnames(mirinterp))  #remove the 'X' from headers (MIR wavenumbers)

# brings file into hyperspec
mirinterp1 <- new("hyperSpec", spc=mirinterp[,grep('[[:digit:]]',colnames(mirinterp))], wavelength= as.numeric(colnames(mirinterp)[grep	('[[:digit:]]',colnames(mirinterp))]),label=list(.wavelength="Wavenumber",spc="Intensity"))
#head(mirinterp1[[]])

#mirinterp2 <- (1/10^mirinterp1)  # converts absorbance spectra into trasnmittance using 1/10^absorb

mirinterp3 <- hyperSpec::spc.loess(mirinterp1, c(seq(6000,600,-4)))  #interpolates spectra between 6000 and 600 cm-1 with a data spacing of 4
plot(mirinterp3, "spc", wl.reverse =F) # plot the transmittance spectra produced between 6000 and 600 cm-1 with a data spacing of 4

output <- mirinterp3[[]]  #generates the mir spectra with correct headings (removes "spc." from the front of wavenumbers)


waves_l <- seq(6000, 600, by= -4)
colnames(output) <- waves_l

#puts sample names back into the data as the first column
ID <- as.data.frame(sample$UniqueID)
#project <- as.data.frame(sample$Project_name)
final <- cbind(ID, output)
#colnames(final)[1] <- "Project_name"
#colnames(final)[2] <- "ID"
head(final)
dim(final)

matplot(x = waves_l, y = t(final[,2:1352]), ylim=c(0,3), type = "l", lty = 1,
        main = "Absorbance - 600 to 6000 & reample with resolution of 4", xlab = "Wavelength (nm)",
        ylab = "Absorbance", col = rep(palette(), each = 3))


#######################################################################################################################################
###  Baseline correction ################################################################################################################
#######################################################################################################################################

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

### Due to error in the function need to delete the first raw 
#spec_a_bc_d <- spec_a_bc_d[-1,] 
#dim(spec_a_bc_d)


## Wave numbers are in reverse order
waves_ss <- seq(600, 6000, by=4)


matplot(x = waves_ss, y = t(spec_a_bc_d), ylim=c(0,3), type = "l", lty = 1,
        main = "Absorbance - baseline corrected", xlab = "Wavelength (nm)",
        ylab = "Absorbance", col = rep(palette(), each = 3))




#################################################################################################################################################################
##################### Prediction of SOC fractions ###############################################################################################################
#################################################################################################################################################################
# All R object models are inside the folder called 
# V:\_Projects\2019 Projects\2019_Farrell_NSW_Forestry\MIR\Analysis\SCaRP_models_20201210
# Copy those objects to your working folder and before you run the analysis i.e. V:\\_Projects\\2019 Projects\\2019_Farrell_NSW_Forestry\\MIR\\Analysis

### When working on Mac, the *.rds files are in the `data/working` folder. Ensure the .gitignore inclides *.rds before committing and pushing

# load the relevant model - in this example its for OC
dir()
model <- readRDS("ROC_model_20201210.rds")



### Carryout predictions
pred <- predict(model, newdata = spec_a_bc_d)
str(pred)

### Notes
## Relevant optimum factors
# sqrt OC = 9
# sqrt POC = 9
# sqrt HOC = 6
# sqrt ROC = 10

pred_d <- as.data.frame(pred)


### Save results
# Remember to extract the predictions associated with the optimum number of factors included in the relevant PLSR model (see above notes)
results <- as.data.frame(cbind(final$`sample$UniqueID`, pred_d$`sqrt_ROC.10 comps`))
head(results)

# Rename the column headings
# change the second heading accordingly 
colnames(results)<- c("UniqueID","sqrtROC")

# Write your results 
write.csv(results, 'sqrtROC_Predictions_20201211.csv')






