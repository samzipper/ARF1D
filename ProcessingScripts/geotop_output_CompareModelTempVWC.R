## geotop_output_CompareModelTempVWC.R
#' This is intended to compare modeled and measured temperature at different
#' depths from a set of soil moisture/temperature sensors.
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
path.mod.VWC <- paste0(git.dir, "geotop/output-tabs/thetaliq0001.txt")
path.obs <- paste0(git.dir, "data/ARFlux/ARFlux_SoilTemp-Moisture_", fire, "_2010-2013.csv")

# save plot path
path.plot <- paste0(git.dir, "geotop/output-plots/CompareModelTempVWC_", version, "-", fire, ".png")

## read in data
df.mod.temp <- read.csv(path.mod.temp, stringsAsFactors=F)
df.mod.VWC <- read.csv(path.mod.VWC, stringsAsFactors=F)
df.obs <- read.csv(path.obs, stringsAsFactors=F)

## pre-processing of observed data
# rename columns to standardized form
colnames(df.obs) <- c("Date", "Year", "Day", "Hour_Min", 
                      "Temp.1", "Temp.2", "Temp.3", "Temp.4", "Temp.5", "Temp.6",
                      "PA.1", "PA.2", "PA.3", "PA.4", "PA.5", "PA.6",
                      "VWC.1", "VWC.2", "VWC.3", "VWC.4", "VWC.5", "VWC.6")

# assign depths based on site 
# (depths are listed on repository, e.g. http://arc-lter.ecosystems.mbl.edu/arfunburned-soiltemp-moisture)
# note that depths may change between years and were not continuously monitored
# depth units are [mm]
if (fire=="Unburned"){
  depths <- c("1"=80, "2"=40, "3"=90, "4"=70, "5"=70, "6"=60)
} else if (fire=="Moderate"){
  depths <- c("1"=40, "2"=70, "3"=60, "4"=100, "5"=50, "6"=50)
} else if (fire=="Severe"){
  depths <- c("1"=70, "2"=70, "3"=50, "4"=55, "5"=90, "6"=45)
}

# calculate mean depth
depth.min <- min(depths)
depth.mean <- mean(depths)
depth.max <- max(depths)

## calculate mean temperature and VWC across all sensors for each day
VWC.cols <- c("VWC.1", "VWC.2", "VWC.3", "VWC.4", "VWC.5", "VWC.6")
Temp.cols <- c("Temp.1", "Temp.2", "Temp.3", "Temp.4", "Temp.5", "Temp.6")

# calculate just Date (without time)
df.obs$Date <- format(mdy_hm(df.obs$Date), "%m/%d/%Y")
Dates.all <- unique(df.obs$Date)
df.obs.day <- data.frame(Date = Dates.all, 
                         VWC.mean = NaN,
                         VWC.std = NaN,
                         VWC.min = NaN,
                         VWC.max = NaN,
                         Temp.mean = NaN,
                         Temp.std = NaN,
                         Temp.min = NaN,
                         Temp.max = NaN)
for (d in Dates.all){
  # figure out index for this date
  i.d <- which(df.obs.day$Date==d)
  
  # summarize VWC and temperature
  df.obs.day$VWC.mean[i.d] <- mean(unlist(df.obs[df.obs$Date==d, VWC.cols]), na.rm=T)/100
  df.obs.day$VWC.std[i.d] <- sd(unlist(df.obs[df.obs$Date==d, VWC.cols]), na.rm=T)/100
  df.obs.day$VWC.min[i.d] <- min(unlist(df.obs[df.obs$Date==d, VWC.cols]), na.rm=T)/100
  df.obs.day$VWC.max[i.d] <- max(unlist(df.obs[df.obs$Date==d, VWC.cols]), na.rm=T)/100
  
  df.obs.day$Temp.mean[i.d] <- mean(unlist(df.obs[df.obs$Date==d, Temp.cols]), na.rm=T)
  df.obs.day$Temp.std[i.d] <- sd(unlist(df.obs[df.obs$Date==d, Temp.cols]), na.rm=T)
  df.obs.day$Temp.min[i.d] <- min(unlist(df.obs[df.obs$Date==d, Temp.cols]), na.rm=T)
  df.obs.day$Temp.max[i.d] <- max(unlist(df.obs[df.obs$Date==d, Temp.cols]), na.rm=T)
}

## Process model output data
# make Date column
df.mod.temp$Date <- format(dmy_hm(df.mod.temp[,1]), "%m/%d/%Y")
df.mod.VWC$Date <- format(dmy_hm(df.mod.VWC[,1]), "%m/%d/%Y")

# keep column for Date and anything beginning with "X" (these are the depths)
all.cols <- colnames(df.mod.temp)
cols.keep <- c("Date", all.cols[startsWith(all.cols, "X")])
df.mod.temp <- df.mod.temp[,cols.keep]
df.mod.VWC <- df.mod.VWC[,cols.keep]

# figure out depth for each column
cols.depth <- as.numeric(sub("X", "", cols.keep))

# determine which columns are within the depth range of the observations
cols.compare <- which(cols.depth > depth.min & cols.depth < depth.max)

# summarize model for each day
Dates.all <- unique(df.mod.temp$Date)
df.mod.day <- data.frame(Date = Dates.all, 
                         VWC.mean = NaN,
                         VWC.std = NaN,
                         VWC.min = NaN,
                         VWC.max = NaN,
                         Temp.mean = NaN,
                         Temp.std = NaN,
                         Temp.min = NaN,
                         Temp.max = NaN)
for (d in Dates.all){
  # figure out index for this date
  i.d <- which(df.mod.day$Date==d)
  
  # summarize VWC and temperature
  df.mod.day$VWC.mean[i.d] <- mean(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare]), na.rm=T)
  df.mod.day$VWC.std[i.d] <- sd(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare]), na.rm=T)
  df.mod.day$VWC.min[i.d] <- min(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare]), na.rm=T)
  df.mod.day$VWC.max[i.d] <- max(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare]), na.rm=T)
  
  df.mod.day$Temp.mean[i.d] <- mean(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
  df.mod.day$Temp.std[i.d] <- sd(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
  df.mod.day$Temp.min[i.d] <- min(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
  df.mod.day$Temp.max[i.d] <- max(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
}

## make plots
# convert dates
df.mod.day$Date <- mdy(df.mod.day$Date)
df.obs.day$Date <- mdy(df.obs.day$Date)

# trim model to 2010-2013
df.mod.day <- subset(df.mod.day, year(Date)>=2010 & year(Date)<=2013)

# temperature plot
p.temp.compare <-
  ggplot() +
  geom_line(data=df.mod.day, aes(x=Date, y=Temp.mean), color="red") +
  geom_point(data=df.obs.day, aes(x=Date, y=Temp.mean)) +
  scale_y_continuous(name="Daily Temperature [C]") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.VWC.compare <-
  ggplot() +
  geom_line(data=df.mod.day, aes(x=Date, y=VWC.mean), color="blue") +
  geom_point(data=df.obs.day, aes(x=Date, y=VWC.mean)) +
  scale_y_continuous(name="Daily VWC [m3/m3]") +
  theme_bw() +
  theme(panel.grid=element_blank())
ggsave(path.plot, arrangeGrob(p.temp.compare, p.VWC.compare, ncol=1),
       width=12, height=8, units="in")