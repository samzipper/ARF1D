## geotop_SoilPropertyCalibration_SelectBest+Plot.R
#' This is intended to compare modeled and measured data:
#'   -Soil temperature
#'   -Soil moisture
#'   -Thaw depth
#'   -Sensible heat flux
#'   -Latent heat flux
#'   -Ground heat flux
#'   -Net radiation

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

require(lubridate)
require(ggplot2)
require(reshape2)
require(dplyr)
require(gridExtra)
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))

#version
version <- "20170228b-Match27n+Porosity"

# function to find closest
which.closest <- function(x, vec){
  # function returns the index (i) of the element in the vector (vec) which is closest to the input value (x).
  # if x is equidistant from multiple elements in x, all are returned.
  # a tolerance parameter is built in to avoid rounding errors.
  # 
  eps <- .Machine$double.eps^0.5  # tolerance parameter
  
  i <- which(abs(vec-x)-min(abs(vec-x), na.rm=T) < eps)
  
  return(i)
}

## start and end year for calibration/validation
cal.start <- 2008
cal.end <- 2010
val.start <- 2011
val.end <- 2012

yr.min <- min(c(cal.start, val.start))
yr.max <- max(c(cal.end, val.end))

# burn severity site
fire <- "Severe"  # options are: Unburned, Moderate, Severe

## read in observed data and preprocess
# paths
path.obs <- paste0(git.dir, "data/ARFlux/ARFlux_SoilTemp-Moisture_", fire, "_2010-2013.csv")
path.obs.ARFlux <- paste0(git.dir, "data/ARFlux/ARFlux-Merged_2008-2012_Daily.csv")
path.thaw <- paste0(git.dir, "data/ARFlux/ARFlux_ThawDepths_2008-2014.csv")

# read in data
df.obs <- read.csv(path.obs, stringsAsFactors=F)
df.obs.ARFlux <- read.csv(path.obs.ARFlux, stringsAsFactors=F)
df.thaw <- read.csv(path.thaw, stringsAsFactors=F)

# add date column
df.obs$Date <- mdy_hm(df.obs$Date)
df.obs.ARFlux$Date <- ymd(df.obs.ARFlux$Date)
df.thaw$Date <- dmy(df.thaw$Date)

# subset to only calibration/validation period
df.obs <- subset(df.obs, year(Date)>=yr.min & year(Date)<=yr.max)
df.obs.ARFlux <- subset(df.obs.ARFlux, year(Date)>=yr.min & year(Date)<=yr.max)
df.thaw <- subset(df.thaw, year(Date)>=yr.min & year(Date)<=yr.max)

## pre-processing of observed data
# trim observations to only this fire
df.obs.ARFlux <- df.obs.ARFlux[df.obs.ARFlux$fire==fire,]
df.thaw <- data.frame(Date = df.thaw$Date,
                      ThawDepth.mm = 10*as.numeric(df.thaw[,which(colnames(df.thaw)==paste0("Thaw.depth.", fire))]))
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

# ARFlux depths are listed on repository, e.g. http://arc-lter.ecosystems.mbl.edu/2012arfluxunburned
depths.ARFlux <- c(20, 60)  # depth is average of 2 cm and 6 cm sensor
depths.ARFlux.VWC <- c(25)  # VWC data is at 2.5 cm

# calculate mean depth
depth.ARFlux.min <- min(depths.ARFlux)
depth.ARFlux.mean <- mean(depths.ARFlux)
depth.ARFlux.max <- max(depths.ARFlux)

# trim ARFlux data frame
df.obs.ARFlux <- df.obs.ARFlux[,c("Date", "Tsoil.C", "Tsoil.C.min", "Tsoil.C.max", "VWC", "VWC.min", "VWC.max", "SWin.W.m2", "SWout.W.m2", "LWin.W.m2","LWout.W.m2", "Rnet.W.m2", "LE.W.m2", "H.W.m2", "G.W.m2")]

## calculate mean temperature and VWC across all sensors for each day
VWC.cols <- c("VWC.1", "VWC.2", "VWC.3", "VWC.4", "VWC.5", "VWC.6")
Temp.cols <- c("Temp.1", "Temp.2", "Temp.3", "Temp.4", "Temp.5", "Temp.6")

# calculate just Date (without time)
df.obs$Date <- format(df.obs$Date, "%m/%d/%Y")
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
# parameter set number
number <- "0001"

# modeled data path
#path.mod.temp <- paste0(git.dir, "geotop/SoilPropertyCalibration/output_", number, "_", fire, "_soiltemp0001.txt")
#path.mod.VWC <- paste0(git.dir, "geotop/SoilPropertyCalibration/output_", number, "_", fire, "_thetaliq0001.txt")
path.mod.temp <- paste0(git.dir, "geotop/output-tabs/soiltemp0001.txt")
path.mod.VWC <- paste0(git.dir, "geotop/output-tabs/thetaliq0001.txt")
path.mod.point <- paste0(git.dir, "geotop/output-tabs/point0001.txt")

# read in data
df.mod.temp <- read.csv(path.mod.temp, stringsAsFactors=F)
df.mod.VWC <- read.csv(path.mod.VWC, stringsAsFactors=F)
df.mod.point <- read.csv(path.mod.point, stringsAsFactors=F)

# make Date column
df.mod.temp$Date <- dmy_hm(df.mod.temp$Date12.DDMMYYYYhhmm.)
df.mod.VWC$Date <- dmy_hm(df.mod.VWC$Date12.DDMMYYYYhhmm.)
df.mod.point$Date <- dmy_hm(df.mod.point$Date12.DDMMYYYYhhmm.)

# subset
df.mod.temp <- subset(df.mod.temp, year(Date)>=yr.min & year(Date)<=yr.max)
df.mod.VWC <- subset(df.mod.VWC, year(Date)>=yr.min & year(Date)<=yr.max)
df.mod.point <- subset(df.mod.point, year(Date)>=yr.min & year(Date)<=yr.max)

# format Date column
df.mod.temp$Date <- format(df.mod.temp$Date, "%m/%d/%Y")
df.mod.VWC$Date <- format(df.mod.VWC$Date, "%m/%d/%Y")
df.mod.point$Date <- format(df.mod.point$Date, "%m/%d/%Y")

# keep column for Date and anything beginning with "X" (these are the depths)
all.cols <- colnames(df.mod.temp)
cols.keep <- c("Date", all.cols[startsWith(all.cols, "X")])
df.mod.temp <- df.mod.temp[,cols.keep]
df.mod.VWC <- df.mod.VWC[,cols.keep]

# figure out depth for each column
cols.depth <- as.numeric(sub("X", "", cols.keep))

# determine which columns are within the depth range of the observations
cols.compare <- unlist(lapply(depths, FUN=which.closest, vec=cols.depth))
cols.compare.ARFlux <- unlist(lapply(depths.ARFlux, FUN=which.closest, vec=cols.depth))
cols.compare.ARFlux.VWC <- unlist(lapply(depths.ARFlux.VWC, FUN=which.closest, vec=cols.depth))

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
                         Temp.max = NaN,
                         Temp.ARFlux.mean = NaN,
                         Temp.ARFlux.std = NaN,
                         Temp.ARFlux.min = NaN,
                         Temp.ARFlux.max = NaN,
                         ThawDepth.mm = NaN)
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
  
  df.mod.day$Temp.ARFlux.mean[i.d] <- mean(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
  df.mod.day$Temp.ARFlux.std[i.d] <- sd(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
  df.mod.day$Temp.ARFlux.min[i.d] <- min(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
  df.mod.day$Temp.ARFlux.max[i.d] <- max(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
  
  df.mod.day$VWC.ARFlux.mean[i.d] <- mean(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)
  df.mod.day$VWC.ARFlux.std[i.d] <- sd(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)
  df.mod.day$VWC.ARFlux.min[i.d] <- min(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)
  df.mod.day$VWC.ARFlux.max[i.d] <- max(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)
  
  ## figure out thaw depth: shallowest depth at which temperature > 0
  # extract all temperatures for that date
  df.mod.temp.profile <- data.frame(depth = cols.depth[is.finite(cols.depth)], 
                                    Temp = unlist(df.mod.temp[df.mod.temp$Date==d, is.finite(cols.depth)]))
  
  # find the shallowest temperature that is frozen
  df.mod.day$ThawDepth.mm[i.d] <- min(df.mod.temp.profile$depth[which(df.mod.temp.profile$Temp<0)])
}

## make plots
# convert dates
df.mod.day$Date <- mdy(df.mod.day$Date)
df.mod.point$Date <- mdy(df.mod.point$Date)
df.obs.day$Date <- mdy(df.obs.day$Date)
df.obs.ARFlux$Date <- ymd(df.obs.ARFlux$Date)

# for each thaw date, summarize to mean and range
df.thaw.day <- summarize(group_by(df.thaw, Date),
                         ThawDepth.mm.mean = mean(ThawDepth.mm, na.rm=T),
                         ThawDepth.mm.std = sd(ThawDepth.mm, na.rm=T),
                         ThawDepth.mm.min = min(ThawDepth.mm, na.rm=T),
                         ThawDepth.mm.max = max(ThawDepth.mm, na.rm=T))


## calculate fit metrics - calibration period
# ARFlux
df.fit.temp.ARFlux <- df.obs.ARFlux[is.finite(df.obs.ARFlux$Tsoil.C),]
df.fit.temp.ARFlux$Temp.mod <- df.mod.day$Temp.ARFlux.mean[match(df.fit.temp.ARFlux$Date, df.mod.day$Date)]

df.fit.VWC.ARFlux <- df.obs.ARFlux[is.finite(df.obs.ARFlux$VWC),]
df.fit.VWC.ARFlux$VWC.mod <- df.mod.day$VWC.ARFlux.mean[match(df.fit.VWC.ARFlux$Date, df.mod.day$Date)]

df.thaw.day$thaw.mod <- df.mod.day$ThawDepth.mm[match(df.thaw.day$Date, df.mod.day$Date)]

df.fit.LE <- df.obs.ARFlux[is.finite(df.obs.ARFlux$LE.W.m2),]
df.fit.LE$LE.mod <- df.mod.point$LE.W.m2.[match(df.fit.LE$Date, df.mod.point$Date)]

df.fit.H <- df.obs.ARFlux[is.finite(df.obs.ARFlux$H.W.m2),]
df.fit.H$H.mod <- df.mod.point$H.W.m2.[match(df.fit.H$Date, df.mod.point$Date)]

df.fit.Rnet <- df.obs.ARFlux[is.finite(df.obs.ARFlux$Rnet.W.m2),]
df.fit.Rnet$Rnet.mod <- df.mod.point$SWnet.W.m2.[match(df.fit.Rnet$Date, df.mod.point$Date)]+df.mod.point$LWnet.W.m2.[match(df.fit.Rnet$Date, df.mod.point$Date)]

df.fit.SWin <- df.obs.ARFlux[is.finite(df.obs.ARFlux$SWin.W.m2),]
df.fit.SWin$SWin.mod <- df.mod.point$SWin.W.m2.[match(df.fit.SWin$Date, df.mod.point$Date)]

df.fit.LWin <- df.obs.ARFlux[is.finite(df.obs.ARFlux$LWin.W.m2),]
df.fit.LWin$LWin.mod <- df.mod.point$LWin.W.m2.[match(df.fit.LWin$Date, df.mod.point$Date)]

df.fit.SWout <- df.obs.ARFlux[is.finite(df.obs.ARFlux$SWout.W.m2),]
df.fit.SWout$SWout.mod <- df.mod.point$SWup.W.m2.[match(df.fit.SWout$Date, df.mod.point$Date)]

df.fit.LWout <- df.obs.ARFlux[is.finite(df.obs.ARFlux$LWout.W.m2),]
df.fit.LWout$LWout.mod <- df.mod.point$LWup.W.m2.[match(df.fit.LWout$Date, df.mod.point$Date)]

df.fit.G <- df.obs.ARFlux[is.finite(df.obs.ARFlux$G.W.m2),]
df.fit.G$G.mod <- df.mod.point$Soil_heat_flux.W.m2.[match(df.fit.G$Date, df.mod.point$Date)]

# temp + VWC
df.fit.temp <- df.obs.day[is.finite(df.obs.day$Temp.mean),]
df.fit.temp$Temp.mod <- df.mod.day$Temp.mean[match(df.fit.temp$Date, df.mod.day$Date)]
df.fit.VWC <- df.obs.day[is.finite(df.obs.day$VWC.mean),]
df.fit.VWC$VWC.mod <- df.mod.day$VWC.mean[match(df.fit.VWC$Date, df.mod.day$Date)]

# make calibration/validation column
df.fit.temp.ARFlux$period <- "cal"
df.fit.temp.ARFlux$period[year(df.fit.temp.ARFlux$Date) >= val.start & year(df.fit.temp.ARFlux$Date) <= val.end] <- "val"

df.fit.VWC.ARFlux$period <- "cal"
df.fit.VWC.ARFlux$period[year(df.fit.VWC.ARFlux$Date) >= val.start & year(df.fit.VWC.ARFlux$Date) <= val.end] <- "val"

df.fit.temp$period <- "cal"
df.fit.temp$period[year(df.fit.temp$Date) >= val.start & year(df.fit.temp$Date) <= val.end] <- "val"

df.fit.VWC$period <- "cal"
df.fit.VWC$period[year(df.fit.VWC$Date) >= val.start & year(df.fit.VWC$Date) <= val.end] <- "val"

df.thaw.day$period <- "cal"
df.thaw.day$period[year(df.thaw.day$Date) >= val.start & year(df.thaw.day$Date) <= val.end] <- "val"

df.fit.LE$period <- "cal"
df.fit.LE$period[year(df.fit.LE$Date) >= val.start & year(df.fit.LE$Date) <= val.end] <- "val"

df.fit.H$period <- "cal"
df.fit.H$period[year(df.fit.H$Date) >= val.start & year(df.fit.H$Date) <= val.end] <- "val"

df.fit.Rnet$period <- "cal"
df.fit.Rnet$period[year(df.fit.Rnet$Date) >= val.start & year(df.fit.Rnet$Date) <= val.end] <- "val"

df.fit.SWin$period <- "cal"
df.fit.SWin$period[year(df.fit.SWin$Date) >= val.start & year(df.fit.SWin$Date) <= val.end] <- "val"

df.fit.LWin$period <- "cal"
df.fit.LWin$period[year(df.fit.LWin$Date) >= val.start & year(df.fit.LWin$Date) <= val.end] <- "val"

df.fit.SWout$period <- "cal"
df.fit.SWout$period[year(df.fit.SWout$Date) >= val.start & year(df.fit.SWout$Date) <= val.end] <- "val"

df.fit.LWout$period <- "cal"
df.fit.LWout$period[year(df.fit.LWout$Date) >= val.start & year(df.fit.LWout$Date) <= val.end] <- "val"

df.fit.G$period <- "cal"
df.fit.G$period[year(df.fit.G$Date) >= val.start & year(df.fit.G$Date) <= val.end] <- "val"

# make a matrix with fit statistics
df.fit.table.cal <- data.frame(RMSE=numeric(5), NRMSE=numeric(5), NSE=numeric(5), R2=numeric(5),
                               row.names=c("Thaw Depth","Tower Soil Temp", "Tower VWC", "Nest Soil Temp", "Nest VWC"))
df.fit.table.cal["Tower Soil Temp",] <- c(RMSE(subset(df.fit.temp.ARFlux, period=="cal")$Temp.mod, subset(df.fit.temp.ARFlux, period=="cal")$Tsoil.C),
                                          NRMSE(subset(df.fit.temp.ARFlux, period=="cal")$Temp.mod, subset(df.fit.temp.ARFlux, period=="cal")$Tsoil.C),
                                          NashSutcliffe(subset(df.fit.temp.ARFlux, period=="cal")$Temp.mod, subset(df.fit.temp.ARFlux, period=="cal")$Tsoil.C),
                                          R2(subset(df.fit.temp.ARFlux, period=="cal")$Temp.mod, subset(df.fit.temp.ARFlux, period=="cal")$Tsoil.C))
df.fit.table.cal["Tower VWC",] <- c(RMSE(subset(df.fit.VWC.ARFlux, period=="cal")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="cal")$VWC),
                                    NRMSE(subset(df.fit.VWC.ARFlux, period=="cal")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="cal")$VWC),
                                    NashSutcliffe(subset(df.fit.VWC.ARFlux, period=="cal")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="cal")$VWC),
                                    R2(subset(df.fit.VWC.ARFlux, period=="cal")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="cal")$VWC))
df.fit.table.cal["Thaw Depth",] <- c(RMSE(subset(df.thaw.day, period=="cal")$thaw.mod, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean),
                                     NRMSE(subset(df.thaw.day, period=="cal")$thaw.mod, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean),
                                     NashSutcliffe(subset(df.thaw.day, period=="cal")$thaw.mod, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean),
                                     R2(subset(df.thaw.day, period=="cal")$thaw.mod, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean))
df.fit.table.cal["Nest Soil Temp",] <- c(RMSE(subset(df.fit.temp, period=="cal")$Temp.mod, subset(df.fit.temp, period=="cal")$Temp.mean),
                                         NRMSE(subset(df.fit.temp, period=="cal")$Temp.mod, subset(df.fit.temp, period=="cal")$Temp.mean),
                                         NashSutcliffe(subset(df.fit.temp, period=="cal")$Temp.mod, subset(df.fit.temp, period=="cal")$Temp.mean),
                                         R2(subset(df.fit.temp, period=="cal")$Temp.mod, subset(df.fit.temp, period=="cal")$Temp.mean))
df.fit.table.cal["Nest VWC",] <- c(RMSE(subset(df.fit.VWC, period=="cal")$VWC.mod, subset(df.fit.VWC, period=="cal")$VWC.mean),
                                   NRMSE(subset(df.fit.VWC, period=="cal")$VWC.mod, subset(df.fit.VWC, period=="cal")$VWC.mean),
                                   NashSutcliffe(subset(df.fit.VWC, period=="cal")$VWC.mod, subset(df.fit.VWC, period=="cal")$VWC.mean),
                                   R2(subset(df.fit.VWC, period=="cal")$VWC.mod, subset(df.fit.VWC, period=="cal")$VWC.mean))

df.fit.table.val <- data.frame(RMSE=numeric(5), NRMSE=numeric(5), NSE=numeric(5), R2=numeric(5),
                               row.names=c("Thaw Depth","Tower Soil Temp", "Tower VWC", "Nest Soil Temp", "Nest VWC"))
df.fit.table.val["Tower Soil Temp",] <- c(RMSE(subset(df.fit.temp.ARFlux, period=="val")$Temp.mod, subset(df.fit.temp.ARFlux, period=="val")$Tsoil.C),
                                          NRMSE(subset(df.fit.temp.ARFlux, period=="val")$Temp.mod, subset(df.fit.temp.ARFlux, period=="val")$Tsoil.C),
                                          NashSutcliffe(subset(df.fit.temp.ARFlux, period=="val")$Temp.mod, subset(df.fit.temp.ARFlux, period=="val")$Tsoil.C),
                                          R2(subset(df.fit.temp.ARFlux, period=="val")$Temp.mod, subset(df.fit.temp.ARFlux, period=="val")$Tsoil.C))
df.fit.table.val["Tower VWC",] <- c(RMSE(subset(df.fit.VWC.ARFlux, period=="val")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="val")$VWC),
                                    NRMSE(subset(df.fit.VWC.ARFlux, period=="val")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="val")$VWC),
                                    NashSutcliffe(subset(df.fit.VWC.ARFlux, period=="val")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="val")$VWC),
                                    R2(subset(df.fit.VWC.ARFlux, period=="val")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="val")$VWC))
df.fit.table.val["Thaw Depth",] <- c(RMSE(subset(df.thaw.day, period=="val")$thaw.mod, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean),
                                     NRMSE(subset(df.thaw.day, period=="val")$thaw.mod, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean),
                                     NashSutcliffe(subset(df.thaw.day, period=="val")$thaw.mod, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean),
                                     R2(subset(df.thaw.day, period=="val")$thaw.mod, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean))
df.fit.table.val["Nest Soil Temp",] <- c(RMSE(subset(df.fit.temp, period=="val")$Temp.mod, subset(df.fit.temp, period=="val")$Temp.mean),
                                         NRMSE(subset(df.fit.temp, period=="val")$Temp.mod, subset(df.fit.temp, period=="val")$Temp.mean),
                                         NashSutcliffe(subset(df.fit.temp, period=="val")$Temp.mod, subset(df.fit.temp, period=="val")$Temp.mean),
                                         R2(subset(df.fit.temp, period=="val")$Temp.mod, subset(df.fit.temp, period=="val")$Temp.mean))
df.fit.table.val["Nest VWC",] <- c(RMSE(subset(df.fit.VWC, period=="val")$VWC.mod, subset(df.fit.VWC, period=="val")$VWC.mean),
                                   NRMSE(subset(df.fit.VWC, period=="val")$VWC.mod, subset(df.fit.VWC, period=="val")$VWC.mean),
                                   NashSutcliffe(subset(df.fit.VWC, period=="val")$VWC.mod, subset(df.fit.VWC, period=="val")$VWC.mean),
                                   R2(subset(df.fit.VWC, period=="val")$VWC.mod, subset(df.fit.VWC, period=="val")$VWC.mean))

df.fit.energy.table.cal <- data.frame(RMSE=numeric(8), NRMSE=numeric(8), NSE=numeric(8), R2=numeric(8),
                                      row.names=c("Rnet","SWin", "SWout", "LWin", "LWout", "LE", "H", "G"))
df.fit.energy.table.cal["Rnet",] <- c(RMSE(subset(df.fit.Rnet, period=="cal")$Rnet.mod, subset(df.fit.Rnet, period=="cal")$Rnet.W.m2),
                                      NRMSE(subset(df.fit.Rnet, period=="cal")$Rnet.mod, subset(df.fit.Rnet, period=="cal")$Rnet.W.m2),
                                      NashSutcliffe(subset(df.fit.Rnet, period=="cal")$Rnet.mod, subset(df.fit.Rnet, period=="cal")$Rnet.W.m2),
                                      R2(subset(df.fit.Rnet, period=="cal")$Rnet.mod, subset(df.fit.Rnet, period=="cal")$Rnet.W.m2))
df.fit.energy.table.cal["SWin",] <- c(RMSE(subset(df.fit.SWin, period=="cal")$SWin.mod, subset(df.fit.SWin, period=="cal")$SWin.W.m2),
                                      NRMSE(subset(df.fit.SWin, period=="cal")$SWin.mod, subset(df.fit.SWin, period=="cal")$SWin.W.m2),
                                      NashSutcliffe(subset(df.fit.SWin, period=="cal")$SWin.mod, subset(df.fit.SWin, period=="cal")$SWin.W.m2),
                                      R2(subset(df.fit.SWin, period=="cal")$SWin.mod, subset(df.fit.SWin, period=="cal")$SWin.W.m2))
df.fit.energy.table.cal["SWout",] <- c(RMSE(subset(df.fit.SWout, period=="cal")$SWout.mod, subset(df.fit.SWout, period=="cal")$SWout.W.m2),
                                       NRMSE(subset(df.fit.SWout, period=="cal")$SWout.mod, subset(df.fit.SWout, period=="cal")$SWout.W.m2),
                                       NashSutcliffe(subset(df.fit.SWout, period=="cal")$SWout.mod, subset(df.fit.SWout, period=="cal")$SWout.W.m2),
                                       R2(subset(df.fit.SWout, period=="cal")$SWout.mod, subset(df.fit.SWout, period=="cal")$SWout.W.m2))
df.fit.energy.table.cal["LWin",] <- c(RMSE(subset(df.fit.LWin, period=="cal")$LWin.mod, subset(df.fit.LWin, period=="cal")$LWin.W.m2),
                                      NRMSE(subset(df.fit.LWin, period=="cal")$LWin.mod, subset(df.fit.LWin, period=="cal")$LWin.W.m2),
                                      NashSutcliffe(subset(df.fit.LWin, period=="cal")$LWin.mod, subset(df.fit.LWin, period=="cal")$LWin.W.m2),
                                      R2(subset(df.fit.LWin, period=="cal")$LWin.mod, subset(df.fit.LWin, period=="cal")$LWin.W.m2))
df.fit.energy.table.cal["LWout",] <- c(RMSE(subset(df.fit.LWout, period=="cal")$LWout.mod, subset(df.fit.LWout, period=="cal")$LWout.W.m2),
                                       NRMSE(subset(df.fit.LWout, period=="cal")$LWout.mod, subset(df.fit.LWout, period=="cal")$LWout.W.m2),
                                       NashSutcliffe(subset(df.fit.LWout, period=="cal")$LWout.mod, subset(df.fit.LWout, period=="cal")$LWout.W.m2),
                                       R2(subset(df.fit.LWout, period=="cal")$LWout.mod, subset(df.fit.LWout, period=="cal")$LWout.W.m2))
df.fit.energy.table.cal["LE",] <- c(RMSE(subset(df.fit.LE, period=="cal")$LE.mod, subset(df.fit.LE, period=="cal")$LE.W.m2),
                                    NRMSE(subset(df.fit.LE, period=="cal")$LE.mod, subset(df.fit.LE, period=="cal")$LE.W.m2),
                                    NashSutcliffe(subset(df.fit.LE, period=="cal")$LE.mod, subset(df.fit.LE, period=="cal")$LE.W.m2),
                                    R2(subset(df.fit.LE, period=="cal")$LE.mod, subset(df.fit.LE, period=="cal")$LE.W.m2))
df.fit.energy.table.cal["H",] <- c(RMSE(subset(df.fit.H, period=="cal")$H.mod, subset(df.fit.H, period=="cal")$H.W.m2),
                                   NRMSE(subset(df.fit.H, period=="cal")$H.mod, subset(df.fit.H, period=="cal")$H.W.m2),
                                   NashSutcliffe(subset(df.fit.H, period=="cal")$H.mod, subset(df.fit.H, period=="cal")$H.W.m2),
                                   R2(subset(df.fit.H, period=="cal")$H.mod, subset(df.fit.H, period=="cal")$H.W.m2))
df.fit.energy.table.cal["G",] <- c(RMSE(subset(df.fit.G, period=="cal")$G.mod, subset(df.fit.G, period=="cal")$G.W.m2),
                                   NRMSE(subset(df.fit.G, period=="cal")$G.mod, subset(df.fit.G, period=="cal")$G.W.m2),
                                   NashSutcliffe(subset(df.fit.G, period=="cal")$G.mod, subset(df.fit.G, period=="cal")$G.W.m2),
                                   R2(subset(df.fit.G, period=="cal")$G.mod, subset(df.fit.G, period=="cal")$G.W.m2))

df.fit.energy.table.val <- data.frame(RMSE=numeric(8), NRMSE=numeric(8), NSE=numeric(8), R2=numeric(8),
                                      row.names=c("Rnet","SWin", "SWout", "LWin", "LWout", "LE", "H", "G"))
df.fit.energy.table.val["Rnet",] <- c(RMSE(subset(df.fit.Rnet, period=="val")$Rnet.mod, subset(df.fit.Rnet, period=="val")$Rnet.W.m2),
                                      NRMSE(subset(df.fit.Rnet, period=="val")$Rnet.mod, subset(df.fit.Rnet, period=="val")$Rnet.W.m2),
                                      NashSutcliffe(subset(df.fit.Rnet, period=="val")$Rnet.mod, subset(df.fit.Rnet, period=="val")$Rnet.W.m2),
                                      R2(subset(df.fit.Rnet, period=="val")$Rnet.mod, subset(df.fit.Rnet, period=="val")$Rnet.W.m2))
df.fit.energy.table.val["SWin",] <- c(RMSE(subset(df.fit.SWin, period=="val")$SWin.mod, subset(df.fit.SWin, period=="val")$SWin.W.m2),
                                      NRMSE(subset(df.fit.SWin, period=="val")$SWin.mod, subset(df.fit.SWin, period=="val")$SWin.W.m2),
                                      NashSutcliffe(subset(df.fit.SWin, period=="val")$SWin.mod, subset(df.fit.SWin, period=="val")$SWin.W.m2),
                                      R2(subset(df.fit.SWin, period=="val")$SWin.mod, subset(df.fit.SWin, period=="val")$SWin.W.m2))
df.fit.energy.table.val["SWout",] <- c(RMSE(subset(df.fit.SWout, period=="val")$SWout.mod, subset(df.fit.SWout, period=="val")$SWout.W.m2),
                                       NRMSE(subset(df.fit.SWout, period=="val")$SWout.mod, subset(df.fit.SWout, period=="val")$SWout.W.m2),
                                       NashSutcliffe(subset(df.fit.SWout, period=="val")$SWout.mod, subset(df.fit.SWout, period=="val")$SWout.W.m2),
                                       R2(subset(df.fit.SWout, period=="val")$SWout.mod, subset(df.fit.SWout, period=="val")$SWout.W.m2))
df.fit.energy.table.val["LWin",] <- c(RMSE(subset(df.fit.LWin, period=="val")$LWin.mod, subset(df.fit.LWin, period=="val")$LWin.W.m2),
                                      NRMSE(subset(df.fit.LWin, period=="val")$LWin.mod, subset(df.fit.LWin, period=="val")$LWin.W.m2),
                                      NashSutcliffe(subset(df.fit.LWin, period=="val")$LWin.mod, subset(df.fit.LWin, period=="val")$LWin.W.m2),
                                      R2(subset(df.fit.LWin, period=="val")$LWin.mod, subset(df.fit.LWin, period=="val")$LWin.W.m2))
df.fit.energy.table.val["LWout",] <- c(RMSE(subset(df.fit.LWout, period=="val")$LWout.mod, subset(df.fit.LWout, period=="val")$LWout.W.m2),
                                       NRMSE(subset(df.fit.LWout, period=="val")$LWout.mod, subset(df.fit.LWout, period=="val")$LWout.W.m2),
                                       NashSutcliffe(subset(df.fit.LWout, period=="val")$LWout.mod, subset(df.fit.LWout, period=="val")$LWout.W.m2),
                                       R2(subset(df.fit.LWout, period=="val")$LWout.mod, subset(df.fit.LWout, period=="val")$LWout.W.m2))
df.fit.energy.table.val["LE",] <- c(RMSE(subset(df.fit.LE, period=="val")$LE.mod, subset(df.fit.LE, period=="val")$LE.W.m2),
                                    NRMSE(subset(df.fit.LE, period=="val")$LE.mod, subset(df.fit.LE, period=="val")$LE.W.m2),
                                    NashSutcliffe(subset(df.fit.LE, period=="val")$LE.mod, subset(df.fit.LE, period=="val")$LE.W.m2),
                                    R2(subset(df.fit.LE, period=="val")$LE.mod, subset(df.fit.LE, period=="val")$LE.W.m2))
df.fit.energy.table.val["H",] <- c(RMSE(subset(df.fit.H, period=="val")$H.mod, subset(df.fit.H, period=="val")$H.W.m2),
                                   NRMSE(subset(df.fit.H, period=="val")$H.mod, subset(df.fit.H, period=="val")$H.W.m2),
                                   NashSutcliffe(subset(df.fit.H, period=="val")$H.mod, subset(df.fit.H, period=="val")$H.W.m2),
                                   R2(subset(df.fit.H, period=="val")$H.mod, subset(df.fit.H, period=="val")$H.W.m2))
df.fit.energy.table.val["G",] <- c(RMSE(subset(df.fit.G, period=="val")$G.mod, subset(df.fit.G, period=="val")$G.W.m2),
                                   NRMSE(subset(df.fit.G, period=="val")$G.mod, subset(df.fit.G, period=="val")$G.W.m2),
                                   NashSutcliffe(subset(df.fit.G, period=="val")$G.mod, subset(df.fit.G, period=="val")$G.W.m2),
                                   R2(subset(df.fit.G, period=="val")$G.mod, subset(df.fit.G, period=="val")$G.W.m2))

## make plots
# path to save plots
#path.plot.val <- paste0(git.dir, "geotop/SoilPropertyCalibration/Plot_", number, "_", fire, ".png")
path.plot.val.sub <- paste0(git.dir, "geotop/output-plots/CalVal_Subsurface_", version, "_", fire, ".png")
path.plot.val.sur <- paste0(git.dir, "geotop/output-plots/CalVal_Surface_", version, "_", fire, ".png")

# plot subsurface temperature and VWC
p.thaw.compare.ARFlux <- 
  ggplot() +
  geom_ribbon(data=subset(df.mod.day, year(Date)>=cal.start & year(Date)<=cal.end), aes(x=Date, ymin=-Inf, ymax=Inf), fill="gray90") +
  geom_line(data=subset(df.mod.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=ThawDepth.mm), color="brown") +
  geom_segment(data=subset(df.thaw.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, xend=Date, y=ThawDepth.mm.mean-ThawDepth.mm.std, yend=ThawDepth.mm.mean+ThawDepth.mm.std)) +
  geom_point(data=subset(df.thaw.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=ThawDepth.mm.mean)) +
  scale_x_date(expand=c(0,0)) +
  scale_y_reverse(name="Tower Thaw Depth [mm]") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.temp.compare.ARFlux <-
  ggplot() +
  geom_ribbon(data=subset(df.mod.day, year(Date)>=cal.start & year(Date)<=cal.end), aes(x=Date, ymin=-Inf, ymax=Inf), fill="gray90") +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(data=subset(df.mod.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=Temp.mean), color="red") +
  geom_point(data=df.obs.ARFlux, aes(x=Date, y=Tsoil.C)) +
  scale_x_date(expand=c(0,0)) +
  scale_y_continuous(name=paste0("Tower Soil Temp [C]: ", depth.ARFlux.min, "-", depth.ARFlux.max, " mm")) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.VWC.compare.ARFlux <-
  ggplot() +
  geom_ribbon(data=subset(df.mod.day, year(Date)>=cal.start & year(Date)<=cal.end), aes(x=Date, ymin=-Inf, ymax=Inf), fill="gray90") +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(data=subset(df.mod.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=VWC.ARFlux.mean), color="blue") +
  geom_point(data=df.obs.ARFlux, aes(x=Date, y=VWC)) +
  scale_x_date(expand=c(0,0)) +
  scale_y_continuous(name=paste0("Tower VWC [m3/m3]: ", depths.ARFlux.VWC, " mm")) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.temp.compare <-
  ggplot() +
  geom_ribbon(data=subset(df.mod.day, year(Date)>=cal.start & year(Date)<=cal.end), aes(x=Date, ymin=-Inf, ymax=Inf), fill="gray90") +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(data=subset(df.mod.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=Temp.mean), color="red") +
  geom_point(data=df.obs.day, aes(x=Date, y=Temp.mean)) +
  scale_y_continuous(name=paste0("Nest Soil Temp [C]: ", depth.min, "-", depth.max, " mm")) +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.VWC.compare <-
  ggplot() +
  geom_ribbon(data=subset(df.mod.day, year(Date)>=cal.start & year(Date)<=cal.end), aes(x=Date, ymin=-Inf, ymax=Inf), fill="gray90") +
  geom_line(data=subset(df.mod.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=VWC.mean), color="blue") +
  geom_point(data=df.obs.day, aes(x=Date, y=VWC.mean)) +
  scale_y_continuous(name=paste0("Nest VWC [m3/m3]: ", depth.min, "-", depth.max, " mm")) +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

ggsave(path.plot.val.sub, arrangeGrob(p.thaw.compare.ARFlux, p.temp.compare.ARFlux, p.VWC.compare.ARFlux, p.temp.compare, p.VWC.compare, 
                                  arrangeGrob(tableGrob(round(df.fit.table.cal, 3)), tableGrob(round(df.fit.table.val, 3)), ncol=2), 
                                  ncol=1, heights=c(1,1,1,1,1,0.75)),
       width=12, height=12, units="in")

# surface fluxes
p.LE.compare <-
  ggplot() +
  geom_ribbon(data=subset(df.mod.point, year(Date)>=cal.start & year(Date)<=cal.end), aes(x=Date, ymin=-Inf, ymax=Inf), fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=LE.W.m2.), color="green") +
  geom_point(data=df.fit.LE, aes(x=Date, y=LE.W.m2)) +
  scale_y_continuous(name="Tower LE Flux [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.H.compare <-
  ggplot() +
  geom_ribbon(data=subset(df.mod.point, year(Date)>=cal.start & year(Date)<=cal.end), aes(x=Date, ymin=-Inf, ymax=Inf), fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=H.W.m2.), color="red") +
  geom_point(data=df.fit.H, aes(x=Date, y=H.W.m2)) +
  scale_y_continuous(name="Tower H Flux [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.G.compare <-
  ggplot() +
  geom_ribbon(data=subset(df.mod.point, year(Date)>=cal.start & year(Date)<=cal.end), aes(x=Date, ymin=-Inf, ymax=Inf), fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=Soil_heat_flux.W.m2.), color="brown") +
  geom_point(data=df.fit.G, aes(x=Date, y=G.W.m2)) +
  scale_y_continuous(name="Tower G Flux [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.Rnet.compare <-
  ggplot() +
  geom_ribbon(data=subset(df.mod.point, year(Date)>=cal.start & year(Date)<=cal.end), aes(x=Date, ymin=-Inf, ymax=Inf), fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=SWnet.W.m2.+LWnet.W.m2.), color="brown") +
  geom_point(data=df.fit.Rnet, aes(x=Date, y=Rnet.W.m2)) +
  scale_y_continuous(name="Tower Rnet [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.SWin.compare <-
  ggplot() +
  geom_ribbon(data=subset(df.mod.point, year(Date)>=cal.start & year(Date)<=cal.end), aes(x=Date, ymin=-Inf, ymax=Inf), fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=SWin.W.m2.), color="brown") +
  geom_point(data=df.obs.ARFlux, aes(x=Date, y=SWin.W.m2)) +
  scale_y_continuous(name="Tower SWin [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.SWout.compare <-
  ggplot() +
  geom_ribbon(data=subset(df.mod.point, year(Date)>=cal.start & year(Date)<=cal.end), aes(x=Date, ymin=-Inf, ymax=Inf), fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=SWup.W.m2.), color="brown") +
  geom_point(data=df.obs.ARFlux, aes(x=Date, y=SWout.W.m2)) +
  scale_y_continuous(name="Tower SWout [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.LWin.compare <-
  ggplot() +
  geom_ribbon(data=subset(df.mod.point, year(Date)>=cal.start & year(Date)<=cal.end), aes(x=Date, ymin=-Inf, ymax=Inf), fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=LWin.W.m2.), color="brown") +
  geom_point(data=df.obs.ARFlux, aes(x=Date, y=LWin.W.m2)) +
  scale_y_continuous(name="Tower LWin [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.LWout.compare <-
  ggplot() +
  geom_ribbon(data=subset(df.mod.point, year(Date)>=cal.start & year(Date)<=cal.end), aes(x=Date, ymin=-Inf, ymax=Inf), fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=LWup.W.m2.), color="brown") +
  geom_point(data=df.obs.ARFlux, aes(x=Date, y=LWout.W.m2)) +
  scale_y_continuous(name="Tower LWout [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

ggsave(path.plot.val.sur, arrangeGrob(p.Rnet.compare, p.SWin.compare, p.SWout.compare, p.LWin.compare, p.LWout.compare, p.LE.compare, p.H.compare, p.G.compare, 
                                      arrangeGrob(tableGrob(round(df.fit.energy.table.cal, 3)), tableGrob(round(df.fit.energy.table.val, 3)), ncol=2), 
                                      ncol=1),
       width=12, height=12, units="in")
