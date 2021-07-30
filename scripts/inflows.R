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
gradgreen <- create_pattern_gradient(id="p1", angle=90, colour1="Brown", colour2="Green")
gradgreen$show()

gradgreenlist <- list(gradgreen)

myrect <- stag$rect(x=0, y=0, width=400, height=100, fill=gradgreen)

#encode it
gradgreenlist <- SVGPatternList_to_svg(gradgreenlist)

svgout(filename = "data/working/example-manual.svg", width = 8, height = 4,
       pattern_list = gradgreenlist)

ggplot(data_long, aes(x=Age, y=-Year, height=Deaths, group=Year, fill=as.factor(flag)))+
  geom_density_ridges(stat="identity", alpha=1, scale=3, colour=alpha(0.0001))+
  theme_classic()+
  scale_y_continuous(breaks=c(-1993:-2017), labels=c(1993:2017), name="", position="right")+
  scale_x_discrete(breaks=c("10", "20", "30", "40", "50", "60", "70", "80", "90"))+
  theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), text=element_text(family="Georgia"))+
  labs(title="Trends in deaths from drug poisoning", 
       subtitle="Data from England and Wales 1993-2017", 
       caption="Source: Office for National Statistics\nPlot by @VictimOfMaths")+
  scale_fill_manual(values=names(gradgreenlist), guide=FALSE)

########################################
#Go from here: https://coolbutuseless.github.io/package/devoutsvg/articles/svg-with-gradient-fill.html
########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create SVG gradient pattern definition
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gradient_pattern <- svgpatternsimple::create_pattern_gradient(
  id      = "p1",      # HTML/SVG id to assign to this pattern
  angle   = 90,        # Direction of the gradient
  colour1 = "wheat",   # Starting colour
  colour2 = "darkgreen"  # Final colour
)

# Contents of 'gradient_pattern'
#> <linearGradient id="p1" x1="0%" y1="100%" x2="0%" y2="0%">
#>   <stop style="stop-color:White;stop-opacity:1" offset="0%" />
#>   <stop style="stop-color:#0570b0;stop-opacity:1" offset="100%" />
#> </linearGradient>

# Visualise in viewer in Rstudio
# gradient_pattern$show()

my_pattern_list <- list(
  `#000001` = list(fill = gradient_pattern)
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Render the graph to the 'svgout' device and nominate any patterns to be 
# rendered by the 'svgpatternsimple' package
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
svgout(filename = "svg/test-gradient.svg", pattern_list = my_pattern_list)
ggplot(iris, aes(x=Sepal.Width, y=Species)) +
  geom_density_ridges(alpha=0.66, scale=2, fill='#000001', colour=alpha(0.5)) +
  theme_classic()
#> Picking joint bandwidth of 0.13
invisible(dev.off())    

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read (and cache) the data from the ONS 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 drugs <- readr::read_csv("https://www.ons.gov.uk/visualisations/dvc661/drugs/datadownload.csv"   )
 saveRDS(drugs, "data/drugs.rds")
drugs <- readRDS("data/drugs.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tidy + reshape data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
drugs <- drugs %>%
  mutate(
    Age = case_when(
      Age   == "<10" ~ "9",
      Age   == "90+" ~ "90",
      TRUE  ~  Age
    )
  ) %>%
  tidyr::gather("Year", "Deaths", -Age) %>%
  mutate(
    Age  = as.integer(Age),
    Year = as.integer(Year)
  )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Render the graph to the 'svgout' device and nominate any fill colours to be 
# rendered by the 'svgpatternsimple' package
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
svgout(filename = "svg/DrugDeaths.svg", pattern_list = my_pattern_list, width=4, height=6)
ggplot(drugs, aes(Age, Year, height=Deaths, group=Year)) +
  geom_density_ridges(stat='identity', scale = 3, colour=NA, fill='#000001') + 
  scale_y_reverse(position = 'right', breaks = sort(unique(drugs$Year))) +
  scale_x_continuous(breaks = seq(10, 90, 10)) +
  theme_classic() +
  theme(
    axis.line.y  = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.title.y = element_blank(),
    text         = element_text(family="Georgia")
  ) +
  labs(
    title    = "Trends in deaths from drug poisoning",
    subtitle = "Data from England and Wales 1993-2017",
    caption  = "Source: Office for National Statistics\nPlot by @VictimOfMaths"
  )
invisible(dev.off())

svgout(filename = "svg/inflows.svg", pattern_list = my_pattern_list, width=8, height=12)
ggplot(newDF, aes(x = Dates, y = year, height = FitData, group = year)) +
  geom_ridgeline(stat = "identity", fill = '#000001', alpha = 0.8, scale = 0.003, min_height = 10, size = 0.2, show.legend = FALSE) +
  theme_classic() +
  scale_y_reverse(breaks = c(1895, 1900, 1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010, 2019), 
                  minor_breaks = seq(1895, 2019, 5),
                  expand = c(0,0), name = "", position = "right") +
  scale_x_date(date_breaks = "1 month", minor_breaks = "1 week", labels=date_format("%b"), expand = c(0,0.1), name = "") +
  theme(axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_text(hjust = -1.5), 
        panel.grid.major.y = element_line(color = "black", size = 0.2, linetype = "dotted"),
        panel.grid.minor.y = element_line(color = "black", size = 0.2, linetype = "dotted"))# +
  #scale_fill_manual(values = inflow_col)
invisible(dev.off())
