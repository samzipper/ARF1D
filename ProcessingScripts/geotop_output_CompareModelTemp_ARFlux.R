## geotop_output_CompareModelTemp_ARFlux.R
#' This is intended to compare modeled and measured temperature at different
#' depths using data from the ARFlux sensors themselves.
#' 

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

require(lubridate)
require(ggplot2)
require(reshape2)
require(dplyr)
require(gridExtra)

# version name
version <- "20170214-10yrSpinUp-5mWTD"
fire <- "Unburned"  # options are: Unburned, Moderate, Severe

## paths
# modeled data path
path.mod.temp <- paste0(git.dir, "geotop/output-tabs/soiltemp0001.txt")
path.obs <- paste0(git.dir, "data/ARFlux/ARFlux_2008-2012_Daily.csv")

# save plot path
path.plot <- paste0(git.dir, "geotop/output-plots/CompareModelTemp_ARFlux_", version, "-", fire, ".png")

## read in data
df.mod.temp <- read.csv(path.mod.temp, stringsAsFactors=F)
df.obs <- read.csv(path.obs, stringsAsFactors=F)

# trim observations to only unburned
df.obs <- df.obs[df.obs$fire==fire,]

## pre-processing of observed data
# assign depths based on site 
# (depths are listed on repository, e.g. http://arc-lter.ecosystems.mbl.edu/2012arfluxunburned )
depths <- c(20, 60)  # depth is average of 2 cm and 6 cm sensor

# calculate mean depth
depth.min <- min(depths)
depth.mean <- mean(depths)
depth.max <- max(depths)

# trim data frame
df.obs <- df.obs[,c("Date", "Tsoil.C", "Tsoil.C.min", "Tsoil.C.max")]

## Process model output data
# make Date column
df.mod.temp$Date <- format(dmy_hm(df.mod.temp[,1]), "%m/%d/%Y")

# keep column for Date and anything beginning with "X" (these are the depths)
all.cols <- colnames(df.mod.temp)
cols.keep <- c("Date", all.cols[startsWith(all.cols, "X")])
df.mod.temp <- df.mod.temp[,cols.keep]

# figure out depth for each column
cols.depth <- as.numeric(sub("X", "", cols.keep))

# determine which columns are within the depth range of the observations
cols.compare <- which(cols.depth > depth.min & cols.depth < depth.max)

# summarize model for each day
Dates.all <- unique(df.mod.temp$Date)
df.mod.day <- data.frame(Date = Dates.all,
                         Temp.mean = NaN,
                         Temp.std = NaN,
                         Temp.min = NaN,
                         Temp.max = NaN)
for (d in Dates.all){
  # figure out index for this date
  i.d <- which(df.mod.day$Date==d)
  
  # summarize VWC and temperature
  df.mod.day$Temp.mean[i.d] <- mean(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
  df.mod.day$Temp.std[i.d] <- sd(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
  df.mod.day$Temp.min[i.d] <- min(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
  df.mod.day$Temp.max[i.d] <- max(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
}

## make plots
# convert dates
df.mod.day$Date <- mdy(df.mod.day$Date)
df.obs$Date <- ymd(df.obs$Date)

# trim model to 2010-2013
df.mod.day <- subset(df.mod.day, year(Date)>=2008 & year(Date)<=2012)

# temperature plot
p.temp.compare <-
  ggplot() +
  geom_line(data=df.mod.day, aes(x=Date, y=Temp.mean), color="red") +
  geom_point(data=df.obs, aes(x=Date, y=Tsoil.C)) +
  scale_y_continuous(name="Daily Soil Temperature [C]") +
  theme_bw() +
  theme(panel.grid=element_blank())
ggsave(path.plot, p.temp.compare,
       width=12, height=4, units="in")