## https://github.com/VictimOfMaths/DeathsOfDespair/blob/master/ONSRidgePlots.R
install.packages("devtools")
install.packages("data.table")

devtools::install_github("coolbutuseless/lofi")      # Colour encoding
devtools::install_github("coolbutuseless/minisvg")   # SVG creation
devtools::install_github("coolbutuseless/devout")    # Device interface
devtools::install_github("coolbutuseless/devoutsvg") # This package
devtools::install_github("coolbutuseless/poissoned") # This package
devtools::install_github("coolbutuseless/svgpatternsimple") # This package

library(data.table)
library(tidyr)
library(ggplot2)
library(ggridges)
library(dplyr)
library(lofi)
library(minisvg)
library(devout)
library(devoutsvg)
library(svgpatternsimple)
library(poissoned)


read.url <- function(url, ...){
  tmpFile <- tempfile()
  download.file(url, destfile = tmpFile, method = "curl")
  url.data <- read.csv(tmpFile, ...)
  return(url.data)
}
data <- read.url("https://www.ons.gov.uk/visualisations/dvc661/drugs/datadownload.csv")

colnames(data) <- c("Age", 1993:2017)
data$Age <- as.character(data$Age)
data$Age <- case_when(
  data$Age=="<10" ~ "9",
  data$Age=="90+" ~ "90",
  TRUE ~ data$Age)

data_long <- melt(as.data.table(data), id.vars = "Age", variable.name = "Year", value.name = "Deaths")
str(data_long)
data_long[ , Year := as.integer(as.vector(Year))]
str(data_long)
data_long[,flag:=1]
str(data_long)

#Create pattern
gradgreen <- svgpatternsimple::create_pattern_gradient(id="p1", angle=90, colour1="Brown", colour2="Green")
gradgreen$show()
