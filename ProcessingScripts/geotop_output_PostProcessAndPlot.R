## geotop_PostProcessAndPlot.R
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
#git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"
git.dir <- "C:/Users/Sam/WorkGits/ARF1D/"

require(lubridate)
require(ggplot2)
require(reshape2)
require(dplyr)
require(gridExtra)
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))

#version
version <- "20170327-1hr"

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

## calibrate with odd years, validate with even
# this is done because the only winter with observational data
# is 2008-2009, and I want to have winter in both calibration
# and validation
cal <- c(2008, 2010, 2012)
val <- c(2009, 2011, 2013)

yr.min <- min(c(cal, val))
yr.max <- max(c(cal, val))

# burn severity site
fire <- "Unburned"  # options are: Unburned, Moderate, Severe

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

# divide VWC to fractional
df.obs[,VWC.cols] <- df.obs[,VWC.cols]/100

# melt to long-form data frame
df.obs.melt <- melt(df.obs[,c("Date", VWC.cols, Temp.cols)], id=c("Date"))

# make column for sensor number
var.strsep <- matrix(unlist(strsplit(as.character(df.obs.melt$variable), ".", fixed=T)), 
                     ncol=2, nrow=dim(df.obs.melt)[1], byrow=T)
df.obs.melt$var <- var.strsep[,1]
df.obs.melt$sensor <- var.strsep[,2]
df.obs.melt$depth <- depths[match(df.obs.melt$sensor, names(depths))]

# summarize to daily mean
df.obs.melt.d <- summarize(group_by(df.obs.melt, Date, variable, var, sensor, depth),
                           obs.mean = mean(value),
                           obs.std = sd(value),
                           obs.min = min(value),
                           obs.max = max(value))

# calculate albedo
df.obs.ARFlux$albedo <- df.obs.ARFlux$SWout.W.m2/df.obs.ARFlux$SWin.W.m2

# any days where albedo is <0 or >1, set to SWout, SWin, and albedo to NaN
df.obs.ARFlux$SWout.W.m2[df.obs.ARFlux$albedo < 0 | df.obs.ARFlux$albedo > 1] <- NaN
df.obs.ARFlux$SWin.W.m2[df.obs.ARFlux$albedo < 0 | df.obs.ARFlux$albedo > 1] <- NaN
df.obs.ARFlux$albedo[df.obs.ARFlux$albedo < 0 | df.obs.ARFlux$albedo > 1] <- NaN

## Process model output data
# modeled data path
path.mod.temp <- paste0(git.dir, "geotop/output-tabs/soiltemp0001.txt")
path.mod.VWC <- paste0(git.dir, "geotop/output-tabs/thetaliq0001.txt")
path.mod.point <- paste0(git.dir, "geotop/output-tabs/point0001.txt")
path.mod.discharge <- paste0(git.dir, "geotop/output-tabs/discharge.txt")

# read in data
df.mod.temp <- read.csv(path.mod.temp, stringsAsFactors=F)
df.mod.VWC <- read.csv(path.mod.VWC, stringsAsFactors=F)
df.mod.point <- read.csv(path.mod.point, stringsAsFactors=F)
df.mod.discharge <- read.csv(path.mod.discharge, stringsAsFactors=F)

# make Date column
df.mod.temp$Date <- dmy_hm(df.mod.temp$Date12.DDMMYYYYhhmm.)
df.mod.VWC$Date <- dmy_hm(df.mod.VWC$Date12.DDMMYYYYhhmm.)
df.mod.point$Date <- dmy_hm(df.mod.point$Date12.DDMMYYYYhhmm.)
df.mod.discharge$Date <- dmy_hm(df.mod.discharge$DATE.day.month.year.hour.min.)

# subset
df.mod.temp <- subset(df.mod.temp, year(Date)>=yr.min & year(Date)<=yr.max)
df.mod.VWC <- subset(df.mod.VWC, year(Date)>=yr.min & year(Date)<=yr.max)
df.mod.point <- subset(df.mod.point, year(Date)>=yr.min & year(Date)<=yr.max)
df.mod.discharge <- subset(df.mod.discharge, year(Date)>=yr.min & year(Date)<=yr.max)

# format Date column
df.mod.temp$Date <- format(df.mod.temp$Date, "%m/%d/%Y")
df.mod.VWC$Date <- format(df.mod.VWC$Date, "%m/%d/%Y")
df.mod.point$Date <- format(df.mod.point$Date, "%m/%d/%Y")
df.mod.discharge$Date <- format(df.mod.discharge$Date, "%m/%d/%Y")

# merge discharge with point file
df.mod.point$Qoutbottom.mm <- df.mod.discharge$Qoutbottom.m3.s.*(24*60*60)*1000  # convert to mm/day

# for profile files, keep column for Date and anything beginning with "X" (these are the depths)
all.cols <- colnames(df.mod.temp)
cols.keep <- c("Date", all.cols[startsWith(all.cols, "X")])
df.mod.temp <- df.mod.temp[,cols.keep]
df.mod.VWC <- df.mod.VWC[,cols.keep]

# reformat as long-form
df.mod.temp.melt <- melt(df.mod.temp, id="Date", value.name="value.mod")
df.mod.temp.melt$depth <- as.numeric(sub("X", "", df.mod.temp.melt$variable))
df.mod.temp.melt.d <- summarize(group_by(df.mod.temp.melt, Date),
                                value.mod.all = mean(value.mod))

df.mod.VWC.melt <- melt(df.mod.VWC, id="Date", value.name="value.mod")
df.mod.VWC.melt$depth <- as.numeric(sub("X", "", df.mod.VWC.melt$variable))
df.mod.VWC.melt.d <- summarize(group_by(df.mod.VWC.melt, Date),
                                value.mod.all = mean(value.mod))

# add column names to match df.obs.melt.d
df.mod.temp.melt$var <- "Temp"
df.mod.VWC.melt$var <- "VWC"

# combine
df.mod.melt <- rbind(df.mod.temp.melt, df.mod.VWC.melt)
df.mod.melt$variable <- NULL

# merge with obs data frame
df.nest.melt.d <- merge(df.obs.melt.d, df.mod.melt, by=c("Date", "var", "depth"), all.x=T)
df.nest.melt.d$variable <- NULL

# get rid of missing data
df.nest.melt.d <- df.nest.melt.d[complete.cases(df.nest.melt.d),]

# summarize to average of all sensors for each day
df.nest.melt.d.all <- summarize(group_by(df.nest.melt.d, Date, var),
                                obs.mean.all=mean(obs.mean),
                                obs.std.all=sd(obs.std),
                                mod.mean.all=mean(value.mod),
                                mod.std.all=sd(value.mod))

# figure out depth for each column
cols.depth <- as.numeric(sub("X", "", cols.keep))

# determine which columns are within the depth range of the observations
cols.compare.ARFlux <- unlist(lapply(depths.ARFlux, FUN=which.closest, vec=cols.depth))
cols.compare.ARFlux.VWC <- unlist(lapply(depths.ARFlux.VWC, FUN=which.closest, vec=cols.depth))

# summarize model for each day
Dates.all <- unique(df.mod.temp$Date)
df.mod.day <- data.frame(Date = Dates.all, 
                         Temp.ARFlux.mean = NaN,
                         Temp.ARFlux.std = NaN,
                         Temp.ARFlux.min = NaN,
                         Temp.ARFlux.max = NaN,
                         VWC.ARFlux.mean = NaN,
                         VWC.ARFlux.std = NaN,
                         VWC.ARFlux.min = NaN,
                         VWC.ARFlux.max = NaN)
for (d in Dates.all){
  # figure out index for this date
  i.d <- which(df.mod.day$Date==d)
  
  # summarize VWC and temperature
  df.mod.day$Temp.ARFlux.mean[i.d] <- mean(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
  df.mod.day$Temp.ARFlux.std[i.d] <- sd(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
  df.mod.day$Temp.ARFlux.min[i.d] <- min(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
  df.mod.day$Temp.ARFlux.max[i.d] <- max(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
  
  df.mod.day$VWC.ARFlux.mean[i.d] <- mean(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)
  df.mod.day$VWC.ARFlux.std[i.d] <- sd(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)
  df.mod.day$VWC.ARFlux.min[i.d] <- min(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)
  df.mod.day$VWC.ARFlux.max[i.d] <- max(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)

}

# calculate albedo
df.mod.point$albedo <- df.mod.point$SWup.W.m2./df.mod.point$SWin.W.m2.

## make plots
# convert dates
df.mod.day$Date <- mdy(df.mod.day$Date)
df.mod.point$Date <- mdy(df.mod.point$Date)
df.mod.temp.melt.d$Date <- mdy(df.mod.temp.melt.d$Date)
df.mod.VWC.melt.d$Date <- mdy(df.mod.VWC.melt.d$Date)
df.nest.melt.d$Date <- mdy(df.nest.melt.d$Date)
df.nest.melt.d.all$Date <- mdy(df.nest.melt.d.all$Date)
df.obs.ARFlux$Date <- ymd(df.obs.ARFlux$Date)

# for each thaw date, summarize observations to mean and range
df.thaw.day <- summarize(group_by(df.thaw, Date),
                         ThawDepth.mm.mean = mean(ThawDepth.mm, na.rm=T),
                         ThawDepth.mm.std = sd(ThawDepth.mm, na.rm=T),
                         ThawDepth.mm.min = min(ThawDepth.mm, na.rm=T),
                         ThawDepth.mm.max = max(ThawDepth.mm, na.rm=T))

# calculate LE based on evaporation and transpiration
df.mod.point$LEfromET.W.m2 <- (df.mod.point$Evap_surface.mm. + df.mod.point$Trasp_canopy.mm.)*2.45*1e6/86400
df.mod.point$LE.weighted.W.m2 <- (df.mod.point$LEg_unveg.W.m2.*(1-df.mod.point$Canopy_fraction...))+(df.mod.point$LEv.W.m2.*df.mod.point$Canopy_fraction...)
df.mod.point$H.weighted.W.m2 <- (df.mod.point$Hg_unveg.W.m2.*(1-df.mod.point$Canopy_fraction...))+(df.mod.point$Hv.W.m2.*df.mod.point$Canopy_fraction...)

# calculate water balance components
df.mod.point$ET.mm <- df.mod.point$Evap_surface.mm. + df.mod.point$Trasp_canopy.mm.
df.mod.point$precip.mm <- df.mod.point$Psnow_over_canopy.mm. + df.mod.point$Prain_over_canopy.mm.
df.mod.point$precip.effective.mm <- df.mod.point$Prain_under_canopy.mm. + df.mod.point$snow_melted.mm.

# summarize to annual water balance
df.mod.point$year <- year(df.mod.point$Date)
df.mod.point.yr <- dplyr::summarize(group_by(df.mod.point, year),
                                    ET.mm = sum(ET.mm),
                                    precip.mm = sum(precip.mm),
                                    precip.effective.mm = sum(precip.effective.mm),
                                    Qoutbottom.mm = sum(Qoutbottom.mm))

## calculate fit metrics - calibration period
# ARFlux
df.fit.temp.ARFlux <- df.obs.ARFlux[is.finite(df.obs.ARFlux$Tsoil.C),]
df.fit.temp.ARFlux$Temp.mod <- df.mod.day$Temp.ARFlux.mean[match(df.fit.temp.ARFlux$Date, df.mod.day$Date)]

df.fit.VWC.ARFlux <- df.obs.ARFlux[is.finite(df.obs.ARFlux$VWC),]
df.fit.VWC.ARFlux$VWC.mod <- df.mod.day$VWC.ARFlux.mean[match(df.fit.VWC.ARFlux$Date, df.mod.day$Date)]

df.thaw.day$thaw.mod.point <- df.mod.point$lowest_thawed_soil_depth.mm.[match(df.thaw.day$Date, df.mod.point$Date)]

df.fit.LE <- df.obs.ARFlux[is.finite(df.obs.ARFlux$LE.W.m2),]
df.fit.LE$LE.mod <- df.mod.point$LE.weighted.W.m2[match(df.fit.LE$Date, df.mod.point$Date)]

df.fit.H <- df.obs.ARFlux[is.finite(df.obs.ARFlux$H.W.m2),]
df.fit.H$H.mod <- df.mod.point$H.weighted.W.m2[match(df.fit.H$Date, df.mod.point$Date)]

df.fit.Rnet <- df.obs.ARFlux[is.finite(df.obs.ARFlux$Rnet.W.m2),]
df.fit.Rnet$Rnet.mod <- df.mod.point$SWnet.W.m2.[match(df.fit.Rnet$Date, df.mod.point$Date)]+df.mod.point$LWnet.W.m2.[match(df.fit.Rnet$Date, df.mod.point$Date)]

df.fit.SWin <- df.obs.ARFlux[is.finite(df.obs.ARFlux$SWin.W.m2),]
df.fit.SWin$SWin.mod <- df.mod.point$SWin.W.m2.[match(df.fit.SWin$Date, df.mod.point$Date)]

df.fit.LWin <- df.obs.ARFlux[is.finite(df.obs.ARFlux$LWin.W.m2),]
df.fit.LWin$LWin.mod <- df.mod.point$LWin.W.m2.[match(df.fit.LWin$Date, df.mod.point$Date)]

df.fit.SWout <- df.obs.ARFlux[is.finite(df.obs.ARFlux$SWout.W.m2),]
df.fit.SWout$SWout.mod <- df.mod.point$SWup.W.m2.[match(df.fit.SWout$Date, df.mod.point$Date)]

df.fit.albedo <- df.obs.ARFlux[is.finite(df.obs.ARFlux$albedo),]
df.fit.albedo$albedo.mod <- df.mod.point$albedo[match(df.fit.albedo$Date, df.mod.point$Date)]

df.fit.LWout <- df.obs.ARFlux[is.finite(df.obs.ARFlux$LWout.W.m2),]
df.fit.LWout$LWout.mod <- df.mod.point$LWup.W.m2.[match(df.fit.LWout$Date, df.mod.point$Date)]

df.fit.G <- df.obs.ARFlux[is.finite(df.obs.ARFlux$G.W.m2),]
df.fit.G$G.mod <- df.mod.point$Soil_heat_flux.W.m2.[match(df.fit.G$Date, df.mod.point$Date)]

# make calibration/validation column
df.nest.melt.d$period <- "cal"
df.nest.melt.d$period[year(df.nest.melt.d$Date) %in% val] <- "val"

df.nest.melt.d.all$period <- "cal"
df.nest.melt.d.all$period[year(df.nest.melt.d.all$Date) %in% val] <- "val"

df.fit.temp.ARFlux$period <- "cal"
df.fit.temp.ARFlux$period[year(df.fit.temp.ARFlux$Date) %in% val] <- "val"

df.fit.VWC.ARFlux$period <- "cal"
df.fit.VWC.ARFlux$period[year(df.fit.VWC.ARFlux$Date)  %in% val] <- "val"

df.fit.temp$period <- "cal"
df.fit.temp$period[year(df.fit.temp$Date)  %in% val] <- "val"

df.fit.VWC$period <- "cal"
df.fit.VWC$period[year(df.fit.VWC$Date)  %in% val] <- "val"

df.thaw.day$period <- "cal"
df.thaw.day$period[year(df.thaw.day$Date)  %in% val] <- "val"

df.fit.LE$period <- "cal"
df.fit.LE$period[year(df.fit.LE$Date)  %in% val] <- "val"

df.fit.H$period <- "cal"
df.fit.H$period[year(df.fit.H$Date)  %in% val] <- "val"

df.fit.Rnet$period <- "cal"
df.fit.Rnet$period[year(df.fit.Rnet$Date)  %in% val] <- "val"

df.fit.SWin$period <- "cal"
df.fit.SWin$period[year(df.fit.SWin$Date)  %in% val] <- "val"

df.fit.LWin$period <- "cal"
df.fit.LWin$period[year(df.fit.LWin$Date)  %in% val] <- "val"

df.fit.SWout$period <- "cal"
df.fit.SWout$period[year(df.fit.SWout$Date)  %in% val] <- "val"

df.fit.albedo$period <- "cal"
df.fit.albedo$period[year(df.fit.albedo$Date)  %in% val] <- "val"

df.fit.LWout$period <- "cal"
df.fit.LWout$period[year(df.fit.LWout$Date)  %in% val] <- "val"

df.fit.G$period <- "cal"
df.fit.G$period[year(df.fit.G$Date)  %in% val] <- "val"

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
df.fit.table.cal["Thaw Depth",] <- c(RMSE(subset(df.thaw.day, period=="cal")$thaw.mod.point, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean),
                                     NRMSE(subset(df.thaw.day, period=="cal")$thaw.mod.point, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean),
                                     NashSutcliffe(subset(df.thaw.day, period=="cal")$thaw.mod.point, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean),
                                     R2(subset(df.thaw.day, period=="cal")$thaw.mod.point, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean))
df.fit.table.cal["Nest Soil Temp",] <- c(RMSE(subset(df.nest.melt.d.all, var=="Temp" & period=="cal")$mod.mean.all, subset(df.nest.melt.d.all, var=="Temp" & period=="cal")$obs.mean.all),
                                         NRMSE(subset(df.nest.melt.d.all, var=="Temp" & period=="cal")$mod.mean.all, subset(df.nest.melt.d.all, var=="Temp" & period=="cal")$obs.mean.all),
                                         NashSutcliffe(subset(df.nest.melt.d.all, var=="Temp" & period=="cal")$mod.mean.all, subset(df.nest.melt.d.all, var=="Temp" & period=="cal")$obs.mean.all),
                                         R2(subset(df.nest.melt.d.all, var=="Temp" & period=="cal")$mod.mean.all, subset(df.nest.melt.d.all, var=="Temp" & period=="cal")$obs.mean.all))
df.fit.table.cal["Nest VWC",] <- c(RMSE(subset(df.nest.melt.d.all, var=="VWC" & period=="cal")$mod.mean.all, subset(df.nest.melt.d.all, var=="VWC" & period=="cal")$obs.mean.all),
                                   NRMSE(subset(df.nest.melt.d.all, var=="VWC" & period=="cal")$mod.mean.all, subset(df.nest.melt.d.all, var=="VWC" & period=="cal")$obs.mean.all),
                                   NashSutcliffe(subset(df.nest.melt.d.all, var=="VWC" & period=="cal")$mod.mean.all, subset(df.nest.melt.d.all, var=="VWC" & period=="cal")$obs.mean.all),
                                   R2(subset(df.nest.melt.d.all, var=="VWC" & period=="cal")$mod.mean.all, subset(df.nest.melt.d.all, var=="VWC" & period=="cal")$obs.mean.all))

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
df.fit.table.val["Thaw Depth",] <- c(RMSE(subset(df.thaw.day, period=="val")$thaw.mod.point, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean),
                                     NRMSE(subset(df.thaw.day, period=="val")$thaw.mod.point, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean),
                                     NashSutcliffe(subset(df.thaw.day, period=="val")$thaw.mod.point, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean),
                                     R2(subset(df.thaw.day, period=="val")$thaw.mod.point, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean))
df.fit.table.val["Nest Soil Temp",] <- c(RMSE(subset(df.nest.melt.d.all, var=="Temp" & period=="val")$mod.mean.all, subset(df.nest.melt.d.all, var=="Temp" & period=="val")$obs.mean.all),
                                         NRMSE(subset(df.nest.melt.d.all, var=="Temp" & period=="val")$mod.mean.all, subset(df.nest.melt.d.all, var=="Temp" & period=="val")$obs.mean.all),
                                         NashSutcliffe(subset(df.nest.melt.d.all, var=="Temp" & period=="val")$mod.mean.all, subset(df.nest.melt.d.all, var=="Temp" & period=="val")$obs.mean.all),
                                         R2(subset(df.nest.melt.d.all, var=="Temp" & period=="val")$mod.mean.all, subset(df.nest.melt.d.all, var=="Temp" & period=="val")$obs.mean.all))
df.fit.table.val["Nest VWC",] <- c(RMSE(subset(df.nest.melt.d.all, var=="VWC" & period=="val")$mod.mean.all, subset(df.nest.melt.d.all, var=="VWC" & period=="val")$obs.mean.all),
                                   NRMSE(subset(df.nest.melt.d.all, var=="VWC" & period=="val")$mod.mean.all, subset(df.nest.melt.d.all, var=="VWC" & period=="val")$obs.mean.all),
                                   NashSutcliffe(subset(df.nest.melt.d.all, var=="VWC" & period=="val")$mod.mean.all, subset(df.nest.melt.d.all, var=="VWC" & period=="val")$obs.mean.all),
                                   R2(subset(df.nest.melt.d.all, var=="VWC" & period=="val")$mod.mean.all, subset(df.nest.melt.d.all, var=="VWC" & period=="val")$obs.mean.all))

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
path.plot.val.sub <- paste0(git.dir, "geotop/output-plots/Plots_CalVal_Subsurface_", version, "_", fire, ".png")
path.plot.val.sur <- paste0(git.dir, "geotop/output-plots/Plots_CalVal_Surface_", version, "_", fire, ".png")
path.plot.val.diag <- paste0(git.dir, "geotop/output-plots/Plots_CalVal_Diagnostics_", version, "_", fire, ".png")
path.plot.val.nest <- paste0(git.dir, "geotop/output-plots/Plots_CalVal_Nest_", version, "_", fire, ".png")
path.plot.val.wb <- paste0(git.dir, "geotop/output-plots/Plots_CalVal_WB_", version, "_", fire, ".png")

# water balance plot
p.WB <- 
  ggplot(df.mod.point.yr, aes(x=year, y=precip.effective.mm-ET.mm)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_bar(stat="identity") +
  scale_x_continuous(name="Year") +
  scale_y_continuous(name="Precip - ET [mm]") +
  theme_bw() +
  theme(panel.grid=element_blank())
ggsave(path.plot.val.wb, p.WB, width=6, height=6, units="in")

# plot subsurface temperature and VWC
p.thaw.compare.ARFlux <- 
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=lowest_thawed_soil_depth.mm.), color="brown") +
  geom_segment(data=subset(df.thaw.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, xend=Date, y=ThawDepth.mm.mean-ThawDepth.mm.std, yend=ThawDepth.mm.mean+ThawDepth.mm.std)) +
  geom_point(data=subset(df.thaw.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=ThawDepth.mm.mean), shape=21) +
  scale_x_date(expand=c(0,0)) +
  scale_y_reverse(name="Tower Thaw Depth [mm]") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.temp.compare.ARFlux <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(data=subset(df.mod.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=Temp.ARFlux.mean), color="red") +
  geom_point(data=df.obs.ARFlux, aes(x=Date, y=Tsoil.C), shape=21) +
  scale_x_date(expand=c(0,0)) +
  scale_y_continuous(name=paste0("Tower Soil Temp [C]: Mean ", depth.ARFlux.min, " & ", depth.ARFlux.max, " mm")) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.VWC.compare.ARFlux <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(data=subset(df.mod.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=VWC.ARFlux.mean), color="blue") +
  geom_point(data=df.obs.ARFlux, aes(x=Date, y=VWC), shape=21) +
  scale_x_date(expand=c(0,0)) +
  scale_y_continuous(name=paste0("Tower VWC [m3/m3]: ", depths.ARFlux.VWC, " mm")) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.temp.compare <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(data=df.mod.temp.melt.d, aes(x=Date, y=value.mod.all), color="red") +
  geom_point(data=subset(df.nest.melt.d.all, var=="Temp"), aes(x=Date, y=obs.mean.all), shape=21) +
  scale_y_continuous(name=paste0("Nest Soil Temp [C]: ", depth.min, "-", depth.max, " mm")) +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.VWC.compare <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_line(data=df.mod.VWC.melt.d, aes(x=Date, y=value.mod.all), color="blue") +
  geom_point(data=subset(df.nest.melt.d.all, var=="VWC"), aes(x=Date, y=obs.mean.all), shape=21) +
  scale_y_continuous(name=paste0("Nest VWC [m3/m3]: ", depth.min, "-", depth.max, " mm")) +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

ggsave(path.plot.val.sub, arrangeGrob(p.thaw.compare.ARFlux, p.temp.compare.ARFlux, p.VWC.compare.ARFlux, p.temp.compare, p.VWC.compare, 
                                      arrangeGrob(tableGrob(round(df.fit.table.cal, 3)), tableGrob(round(df.fit.table.val, 3)), ncol=2), 
                                      ncol=1, heights=c(1,1,1,1,1,0.75)),
       width=12, height=12, units="in")

# nest plots, split out by depth
df.mod.melt$Date <- mdy(df.mod.melt$Date)
df.mod.obs.melt <- merge(df.mod.melt, df.nest.melt.d[,c("Date", "depth", "var", "obs.mean")], by=c("Date", "depth", "var"), all.x=T)
df.nest.fit <- summarize(group_by(df.nest.melt.d, var, depth),
                         RMSE=round(RMSE(value.mod, obs.mean), 3),
                         NRMSE=round(NRMSE(value.mod, obs.mean), 3),
                         NSE=round(NashSutcliffe(value.mod, obs.mean), 3),
                         R2=round(R2(value.mod, obs.mean), 3))

p.nest.temp <- 
  ggplot(subset(df.mod.obs.melt, depth %in% depths & var=="Temp"), aes(x=Date)) +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(aes(y=value.mod), color="red") +
  geom_point(aes(y=obs.mean), shape=21) +
  scale_y_continuous(name="Temperature [C]") +
  scale_x_date(name="Date", expand=c(0,0)) +
  facet_grid(depth~.) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.nest.VWC <- 
  ggplot(subset(df.mod.obs.melt, depth %in% depths & var=="VWC"), aes(x=Date)) +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(aes(y=value.mod), color="blue") +
  geom_point(aes(y=obs.mean), shape=21) +
  scale_y_continuous(name="VWC [m3 m-3]") +
  scale_x_date(name="Date", expand=c(0,0)) +
  facet_grid(depth~.) +
  theme_bw() +
  theme(panel.grid=element_blank())

ggsave(path.plot.val.nest, arrangeGrob(arrangeGrob(p.nest.temp, p.nest.VWC, ncol=2), tableGrob(df.nest.fit), 
                                       ncol=1, heights=c(1,0.5)),
       width=12, height=12, units="in")

# surface fluxes
p.LE.compare <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=LE.weighted.W.m2), color="green") +
  geom_point(data=df.fit.LE, aes(x=Date, y=LE.W.m2), shape=21) +
  scale_y_continuous(name="Tower LE Flux [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.LE.scatter <-
  ggplot(df.fit.LE, aes(x=LE.W.m2, y=LE.mod, color=period)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous("Tower LE Flux [W m-2]") +
  scale_y_continuous("Modeled LE Flux [W m-2]") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.H.compare <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=H.weighted.W.m2), color="red") +
  geom_point(data=df.fit.H, aes(x=Date, y=H.W.m2), shape=21) +
  scale_y_continuous(name="Tower H Flux [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.H.scatter <-
  ggplot(df.fit.H, aes(x=H.W.m2, y=H.mod, color=period)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous("Tower H Flux [W m-2]") +
  scale_y_continuous("Modeled H Flux [W m-2]") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.G.compare <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=Soil_heat_flux.W.m2.), color="brown") +
  geom_point(data=df.fit.G, aes(x=Date, y=G.W.m2), shape=21) +
  scale_y_continuous(name="Tower G Flux [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.G.scatter <-
  ggplot(df.fit.G, aes(x=G.W.m2, y=G.mod, color=period)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous("Tower G Flux [W m-2]") +
  scale_y_continuous("Modeled G Flux [W m-2]") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.Rnet.compare <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=SWnet.W.m2.+LWnet.W.m2.), color="brown") +
  geom_point(data=df.fit.Rnet, aes(x=Date, y=Rnet.W.m2), shape=21) +
  scale_y_continuous(name="Tower Rnet [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.Rnet.scatter <-
  ggplot(df.fit.Rnet, aes(x=Rnet.W.m2, y=Rnet.mod, color=period)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous("Tower Rnet Flux [W m-2]") +
  scale_y_continuous("Modeled Rnet Flux [W m-2]") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.SWin.compare <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=SWin.W.m2.), color="brown") +
  geom_point(data=df.obs.ARFlux, aes(x=Date, y=SWin.W.m2), shape=21) +
  scale_y_continuous(name="Tower SWin [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.SWin.scatter <-
  ggplot(df.fit.SWin, aes(x=SWin.W.m2, y=SWin.mod, color=period)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous("Tower SWin Flux [W m-2]") +
  scale_y_continuous("Modeled SWin Flux [W m-2]") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.SWout.compare <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=SWup.W.m2.), color="brown") +
  geom_point(data=df.obs.ARFlux, aes(x=Date, y=SWout.W.m2), shape=21) +
  scale_y_continuous(name="Tower SWout [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.SWout.scatter <-
  ggplot(df.fit.SWout, aes(x=SWout.W.m2, y=SWout.mod, color=period)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous("Tower SWout Flux [W m-2]") +
  scale_y_continuous("Modeled SWout Flux [W m-2]") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.albedo.compare <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=albedo), color="brown") +
  geom_point(data=df.obs.ARFlux, aes(x=Date, y=albedo), shape=21) +
  scale_y_continuous(name="Tower Albedo [-]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.albedo.scatter <-
  ggplot(df.fit.albedo, aes(x=albedo, y=albedo.mod, color=period)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous("Tower Albedo") +
  scale_y_continuous("Modeled Albedo") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.LWin.compare <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=LWin.W.m2.), color="brown") +
  geom_point(data=df.obs.ARFlux, aes(x=Date, y=LWin.W.m2), shape=21) +
  scale_y_continuous(name="Tower LWin [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.LWin.scatter <-
  ggplot(df.fit.LWin, aes(x=LWin.W.m2, y=LWin.mod, color=period)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous("Tower LWin Flux [W m-2]") +
  scale_y_continuous("Modeled LWin Flux [W m-2]") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.LWout.compare <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=LWup.W.m2.), color="brown") +
  geom_point(data=df.obs.ARFlux, aes(x=Date, y=LWout.W.m2), shape=21) +
  scale_y_continuous(name="Tower LWout [W m-2]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.LWout.scatter <-
  ggplot(df.fit.LWout, aes(x=LWout.W.m2, y=LWout.mod, color=period)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous("Tower LWout Flux [W m-2]") +
  scale_y_continuous("Modeled LWout Flux [W m-2]") +
  theme_bw() +
  theme(panel.grid=element_blank())

ggsave(path.plot.val.sur, arrangeGrob(p.Rnet.compare, p.SWin.compare, p.SWout.compare, p.albedo.compare, p.LWin.compare, p.LWout.compare, p.LE.compare, p.H.compare, p.G.compare, 
                                      arrangeGrob(tableGrob(round(df.fit.energy.table.cal, 3)), tableGrob(round(df.fit.energy.table.val, 3)), ncol=2), 
                                      ncol=1),
       width=12, height=16, units="in")

ggsave(path.plot.val.sur, arrangeGrob(p.Rnet.compare, p.Rnet.scatter,
                                      p.SWin.compare, p.SWin.scatter,
                                      p.SWout.compare, p.SWout.scatter,
                                      p.albedo.compare, p.albedo.scatter,
                                      p.LWin.compare, p.LWin.scatter,
                                      p.LWout.compare, p.LWout.scatter,
                                      p.LE.compare, p.LE.scatter,
                                      p.H.compare, p.H.scatter,
                                      p.G.compare, p.G.scatter,
                                      tableGrob(round(df.fit.energy.table.cal, 3)), tableGrob(round(df.fit.energy.table.val, 3)), 
                                      ncol=2, widths=c(3,1)),
       width=16, height=16, units="in")

# other diagnostic plots
p.mod.LSAI <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=LSAI.m2.m2.), color="green") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=Canopy_fraction...), color="brown") +
  scale_y_continuous(name="LSAI (green), Canopy fraction (brown)") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.mod.ET <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_ribbon(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, ymin=0, ymax=Evap_surface.mm.), fill="brown", alpha=0.5) +
  geom_ribbon(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, ymin=Evap_surface.mm., ymax=Evap_surface.mm.+Trasp_canopy.mm.), fill="green", alpha=0.5) +
  scale_y_continuous(name="Transpiration (green), Soil Evap (brown)") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.mod.Tair <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=Tair.C.), color="red") +
#  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=Tvegetation.C.), color="green") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=Tsurface.C.), color="orange") +
  scale_y_continuous(name="Tair (red), Tsurface (orange)") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.mod.precip <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=Prain_over_canopy.mm.+Psnow_over_canopy.mm.), color="blue") +
  scale_y_continuous(name="Daily Precip [mm]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.mod.snow <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_hline(yintercept=0, color="gray65") +
  geom_ribbon(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, ymin=0, ymax=snow_depth.mm.), fill="deepskyblue1") +
  scale_y_continuous(name="Snow Depth [mm]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.mod.WTD <-
  ggplot() +
  annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(data=subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=highest_water_table_depth.mm.), color="blue") +
  scale_y_reverse(name="Water Table Depth [mm]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

ggsave(path.plot.val.diag, arrangeGrob(p.mod.precip, p.mod.Tair, p.mod.ET, p.mod.LSAI, p.mod.snow, p.mod.WTD, ncol=1),
       width=12, height=12, units="in")

# save zoomed-in nice plots
path.plot.nice <- paste0(git.dir, "geotop/output-plots/Plots_CalVal_Zoom_", version, "_", fire, ".png")
ggsave(path.plot.nice, 
       arrangeGrob(
         arrangeGrob(p.mod.precip+scale_x_date(limits=c(as.Date("2007-12-31"), as.Date("2010-01-01")), expand=c(0,0)),
                     p.mod.Tair+scale_x_date(limits=c(as.Date("2007-12-31"), as.Date("2010-01-01")), expand=c(0,0)),
                     p.mod.ET+scale_x_date(limits=c(as.Date("2007-12-31"), as.Date("2010-01-01")), expand=c(0,0)), ncol=1),
         arrangeGrob(p.thaw.compare.ARFlux+scale_x_date(limits=c(as.Date("2007-12-31"), as.Date("2010-01-01")), expand=c(0,0)),
                     p.temp.compare.ARFlux+scale_x_date(limits=c(as.Date("2007-12-31"), as.Date("2010-01-01")), expand=c(0,0)),
                     p.VWC.compare.ARFlux+scale_x_date(limits=c(as.Date("2007-12-31"), as.Date("2010-01-01")), expand=c(0,0)), ncol=1),
         arrangeGrob(p.Rnet.compare+scale_x_date(limits=c(as.Date("2007-12-31"), as.Date("2010-01-01")), expand=c(0,0)),
                     p.SWin.compare+scale_x_date(limits=c(as.Date("2007-12-31"), as.Date("2010-01-01")), expand=c(0,0)),
                     p.SWout.compare+scale_x_date(limits=c(as.Date("2007-12-31"), as.Date("2010-01-01")), expand=c(0,0)), ncol=1),
         arrangeGrob(p.LE.compare+scale_x_date(limits=c(as.Date("2007-12-31"), as.Date("2010-01-01")), expand=c(0,0)),
                     p.H.compare+scale_x_date(limits=c(as.Date("2007-12-31"), as.Date("2010-01-01")), expand=c(0,0)),
                     p.G.compare+scale_x_date(limits=c(as.Date("2007-12-31"), as.Date("2010-01-01")), expand=c(0,0)), ncol=1),
       ncol=4), width=16, height=12, units="in")
