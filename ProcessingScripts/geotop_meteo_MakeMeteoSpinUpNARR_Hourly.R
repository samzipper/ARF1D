## geotop_meteo_MakeMeteoSpinUpNARR_Hourly.R
#' This script is intended to take daily meteo data from NARR (generated with
#' script geotop_meteo_MakeMeteoSpinUpNARR_Daily.R). Two output files will be generated:
#'   -fname.hourly has a simple constant value the entire day based on daily data
#'   -fname.hourly.dynamic statistically downscales to estimated hourly values with diel patterns
#' 
#' NARR link: https://www.esrl.noaa.gov/psd/data/gridded/data.narr.monolevel.html

rm(list=ls())

# git directory for relative paths
#git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"
git.dir <- "C:/Users/Sam/WorkGits/ARF1D/"

# which geotop version
f.geotop <- "geotop_NRCS/"
#f.geotop <- "geotop/"

require(lubridate)
require(ggplot2)
require(dplyr)
source(paste0(git.dir, "ProcessingScripts/Tair_HourlyFromDaily.R"))

# filenames
fname.daily <- paste0(git.dir, f.geotop, "meteo/meteoNARRdailyWithSpinUp0001.txt")   # daily input data
fname.stats <- paste0(git.dir, "data/ARFlux/ARFlux-Merged_MonthlyStats.csv")      # file with monthly statistics for downscaling
fname.hourly<- paste0(git.dir, f.geotop, "meteo/meteoNARRhourlyWithSpinUp0001.txt")  # where to save hourly output data
fname.3hourly<- paste0(git.dir, f.geotop, "meteo/meteoNARR3hourlyWithSpinUp0001.txt")  # where to save 3 hourly output data
fname.hourly.dynamic <- paste0(git.dir, f.geotop, "meteo/meteoNARRhourlyDynamicWithSpinUp0001.txt")  # where to save hourly output data
fname.3hourly.dynamic<- paste0(git.dir, f.geotop, "meteo/meteoNARR3hourlyDynamicWithSpinUp0001.txt")  # where to save 3 hourly output data

# year to start/end hourly data
yr.start <- 1978
yr.end <- 2011

# coordinates of site
site.lat <- 68.93      # unburned=68.93, moderate=68.95, severe=68.99
site.lon <- -150.27    # unburned=-150.27, moderate=-150.21, severe=-150.28
site.elev<- 760        # [m] - site elevation
site.lat.rad <- (pi/180)*site.lat        # convert latitude to radians
Lm.center <- -135           # longitude of center of time zone [deg E of Greenwich] 75, 90, 105, 120, 135 for Eastern, Central, Rocky Mountain, Pacific, Alaska time zones

## read in data
df.in <- read.csv(fname.daily, stringsAsFactors=F)
df.stats <- read.csv(fname.stats, stringsAsFactors=F)

# date processing
df.in$Date <- as.Date(dmy_hm(df.in$POSIX))  # strip hours
df.in$DOY <- yday(df.in$Date)
df.in$Month <- month(df.in$Date)

# subset to start/end data
df.in <- subset(df.in, year(Date)>=yr.start & year(Date)<=yr.end)

## solar time calculations
# equation of time (Campbell & Norman Eq 11.4)
f <- (279.575+0.9856*df.in$DOY)*pi/180
EqTime <- (-104.7*sin(f) + 596.2*sin(2*f)+4.3*sin(3*f)-12.7*sin(4*f)-429.3*cos(f)-2*cos(2*f)+19.3*cos(3*f))/3600
df.in$SolarNoon <- 12 - (site.lon-Lm.center)/15 - EqTime

## calculate total precipitation [mm/d], start time [hr], precipitation duration [hrs], and precipitation intensity [mm/hr]
df.in$prec.mm_d <- df.in$Iprec*24
df.in$prec.duration.hrs <- sample(seq(4,24), dim(df.in)[1], replace=T)
df.in$prec.intensity.mm_hr <- df.in$prec.mm_d/df.in$prec.duration.hrs

sample.cap <- function(lower, duration, cap){
  # this function takes a random sample while ensuring that the sample 
  # plus some duration do not exceed a cap
  #
  # required inputs are:
  #  -lower = lower bound of sample
  #  -duration = how long the event lasts
  #  -cap = upper bound of event + duration
  #
  #lower <- 0
  #cap <- 23
  #duration <- 6
  
  sample(seq(lower,cap-duration), 1)
}

df.in$prec.start.hr <- unlist(lapply(X=df.in$prec.duration.hrs, MARGIN=1, FUN=function(x,...) {sample.cap(duration=x, lower=0, cap=23)}))
df.in$prec.end.hr <- df.in$prec.start.hr + df.in$prec.duration.hrs-1

## merge in monthly diurnal temperature range, calculate max/min
df.in <- merge(df.in, df.stats, by="Month")
df.in$Tair.C.min <- df.in$AirT - (df.in$Tair.C.range/2)
df.in$Tair.C.max <- df.in$AirT + (df.in$Tair.C.range/2)

# also get yesterday's Tmax and tomorrow Tmin, for hourly scaling
df.in <- df.in[order(df.in$Date),]  # order df.in chronologically

df.in$Tair.C.max.yesterday[1] <- df.in$Tair.C.max[1]
df.in$Tair.C.max.yesterday[2:length(df.in$Tair.C.max)] <- df.in$Tair.C.max[1:(length(df.in$Tair.C.max)-1)]

df.in$Tair.C.min.tomorrow[length(df.in$Tair.C.min)] <- df.in$Tair.C.min[length(df.in$Tair.C.min)]
df.in$Tair.C.min.tomorrow[1:(length(df.in$Tair.C.min)-1)] <- df.in$Tair.C.min[2:length(df.in$Tair.C.min)]

# Hourly data -------------------------------------------------------------

## make hourly data
hours <- data.frame(DOY=rep(seq(1,366), each=24),
                    hr = rep(seq(0,23), 366))

# make data frame with hourly data, year, and DOY
df.h <- merge(hours, df.in[,c("DOY","Month","Date","SolarNoon","Swglob",
                              "Iprec", "prec.intensity.mm_hr","prec.start.hr", "prec.end.hr",
                              "WindSp","RH","P",
                              "AirT", "Tair.C.min", "Tair.C.max", "Tair.C.max.yesterday", "Tair.C.min.tomorrow")], by=c("DOY"), all.x=T)
df.h <- df.h[order(df.h$Date, df.h$hr),]

## air temperature
# calculate solar time (12=solar noon)
df.h$hr.solar <- 12 + (df.h$hr - df.h$SolarNoon)
df.h$hr.solar[df.h$hr.solar >= 24] <- df.h$hr.solar[df.h$hr.solar >= 24]-24
df.h$hr.solar[df.h$hr.solar < 0] <- df.h$hr.solar[df.h$hr.solar < 0]+24

# calculate hourly temperature
df.h$Tair.C <- Tair_HourlyFromDaily(solar.time=df.h$hr.solar,
                                    Tmax=df.h$Tair.C.max,
                                    Tmin=df.h$Tair.C.min,
                                    Tmax.yesterday=df.h$Tair.C.max.yesterday,
                                    Tmin.tomorrow=df.h$Tair.C.min.tomorrow)

# summarize to double-check against daily
df.h.d <- summarize(group_by(df.h, Date),
                    Tair.hr = mean(Tair.C))
df.comp <- data.frame(Tair.d = df.in$AirT,
                      Tair.hr = df.h.d$Tair.hr)

# calculate bias and correct
Tair.bias <- coef(lm(Tair.hr ~ Tair.d, data=df.comp))[1]
df.h$Tair.C <- df.h$Tair.C - Tair.bias

## precip
# turn off precip when it's not raining
df.h$prec.intensity.mm_hr[df.h$hr < df.h$prec.start.hr | df.h$hr > df.h$prec.end.hr] <- 0

## wind speed
# randomly distribute wind speed (note: this is not used in output)
df.h$WindSp <- 1.13989*df.h$WindSp*(-log(runif(1,0,1)))^0.30   # this is the way that EPIC model does it, according to AgroIBIS
df.h$WindSp[df.h$WindSp<0] <- 0

## solar radiation - copied from AgroIBIS toa_radiation.f and weather.f
# cosine of zenith angle
orbit <- 2.0 * pi * df.h$DOY / 365.2425
angle <- 2.0 * pi * (df.h$hr.solar - 12.0) / 24.0
xdecl <-  0.006918 - 0.399912 * cos(orbit) + 
  0.070257 * sin(orbit) - 
  0.006758 * cos(2.0 * orbit) + 
  0.000907 * sin(2.0 * orbit) - 
  0.002697 * cos(3.0 * orbit) + 
  0.001480 * sin(3.0 * orbit)
calc_sw <- 1370 * (1.000110 + 0.034221 * cos(orbit) + 0.001280 * sin(orbit) + 0.000719 * cos(2.0 * orbit)  + 0.000077 * sin(2.0 * orbit))
df.h$coszen <- sin(site.lat.rad)*sin(xdecl) + cos(site.lat.rad)*cos(xdecl)*cos(angle)
df.h$coszen[df.h$coszen<0] <- 0

# summarize to daily (will be same every year)
df.coszen.d <- summarize(group_by(df.h, DOY),
                         coszen.mean = mean(coszen))
df.coszen.d$coszen.mean[df.coszen.d$coszen.mean < 0.0001] <- 0.0001

# merge
df.h <- merge(df.h, df.coszen.d, by=c("DOY"), all.x=T)

# solar radiation
df.h$SWin <- df.h$Swglob * df.h$coszen/df.h$coszen.mean
df.h$SWin[df.h$SWin<0] <- 0

# order df.h
df.h <- df.h[order(df.h$Date, df.h$hr),]

# make output data frame
df.out.dynamic <- data.frame(POSIX = format(df.h$Date + hours(df.h$hr), "%d/%m/%Y %H:%M"),
                             Iprec = df.h$prec.intensity.mm_hr,
                             WindSp = df.h$WindSp,
                             AirT = df.h$Tair.C,
                             RH = df.h$RH,
                             #P = df.h$P,
                             Swglob = df.h$SWin,
                             CloudTrans = 1.0)
df.out <- data.frame(POSIX = format(df.h$Date + hours(df.h$hr), "%d/%m/%Y %H:%M"),
                             Iprec = df.h$Iprec,
                             WindSp = df.h$WindSp,
                             AirT = df.h$AirT,
                             RH = df.h$RH,
                             #P = df.h$P,
                             Swglob = df.h$Swglob,
                             CloudTrans = 1.0)
#df.out.3hr <- df.out[seq(1,dim(df.out)[1],3), ]
#df.out.3hr.dynamic <- df.out.dynamic[seq(1,dim(df.out.dynamic)[1],3), ]

# write output
f.hr.dynamic <- file(fname.hourly.dynamic, open="wb")
write.table(df.out.dynamic, file=f.hr.dynamic, quote=F, sep=",", na="-9999.0", row.names=F, eol="\n")
close(f.hr.dynamic)

f.hr <- file(fname.hourly, open="wb")
write.table(df.out, file=f.hr, quote=F, sep=",", na="-9999.0", row.names=F)
close(f.hr)
