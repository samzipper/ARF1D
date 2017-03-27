## geotop_meteo_MakeMeteoSpinUpNARRrepeated_Hourly.R
#' This script is intended to take hourly meteo data from NARR (generated with
#' script geotop_meteo_MakeMeteoSpinUpNARR_Hourly.R), determine what year has an
#' "average" amount of precipitation, and generated a new output file that repeats
#' that year over and over again.
#' 
#' NARR link: https://www.esrl.noaa.gov/psd/data/gridded/data.narr.monolevel.html

rm(list=ls())

# git directory for relative paths
#git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"
git.dir <- "C:/Users/Sam/WorkGits/ARF1D/"

require(lubridate)
require(ggplot2)
require(dplyr)

# filenames
fname.hourly.in <- paste0(git.dir, "geotop/meteo/meteoNARRhourlyDynamicWithSpinUp0001.txt")  # where to save hourly output data
fname.hourly.out <- paste0(git.dir, "geotop/meteo/meteoNARRhourlyDynamicRepeatedWithSpinUp0001.txt")  # where to save hourly output data

# year to start/end hourly data
yr.start <- 2000
yr.end <- 2015

## read in data
df.in <- read.csv(fname.hourly.in, stringsAsFactors=F)

# date processing
df.in$datetime <- dmy_hm(df.in$POSIX)
df.in$Date <- as.Date(dmy_hm(df.in$POSIX))  # strip hours
df.in$DOY <- yday(df.in$datetime)
df.in$Month <- month(df.in$datetime)
df.in$year <- year(df.in$datetime)
df.in$hour <- hour(df.in$datetime)

# summarize by year
df.in.yr <- summarize(group_by(df.in, year),
                      prec.mm = sum(Iprec),
                      Tair.C = mean(AirT))

# make plot
p.annual.precip.Tair <-
  ggplot(df.in.yr, aes(x=prec.mm, y=Tair.C)) +
  geom_point()

# make z-score
df.in.yr$prec.z <- as.numeric(scale(df.in.yr$prec.mm))
df.in.yr$Tair.z <- as.numeric(scale(df.in.yr$Tair.C))

# find year with combined z-score closest to 0
yr.avg <- df.in.yr$year[which.min(abs(df.in.yr$prec.z) + abs(df.in.yr$Tair.z))]

# make a data frame that contains data from start to end
df.avg <- subset(df.in, year==yr.avg)

df.yrs <- data.frame(datetime = seq(ymd_hm(paste0(yr.start-1, "-12-31 23:00")), 
                                    ymd_hm(paste0(yr.end, "-12-31 23:00")),
                                    by="hour"))
df.yrs$DOY <- yday(df.yrs$datetime)
df.yrs$hour <- hour(df.yrs$datetime)

# set any leap year DOY equal to 365
df.yrs$DOY[df.yrs$DOY==366] <- 365

# add data from average year
df.yrs <- merge(df.yrs, df.avg[, !(colnames(df.avg) %in% c("datetime"))], by=c("DOY", "hour"), all.x=T)

# put in order
df.yrs <- df.yrs[order(df.yrs$datetime), ]

# make output data frame
df.out <- data.frame(POSIX=format(df.yrs$datetime, "%d/%m/%Y %H:%M"),
                     Iprec=df.yrs$Iprec,
                     WindSp=df.yrs$WindSp,
                     AirT=df.yrs$AirT,
                     RH=df.yrs$RH,
                     Swglob=df.yrs$Swglob,
                     CloudTrans=1.0)

# write output
f.out <- file(fname.hourly.out, open="wb")
write.table(df.out, file=f.out, quote=F, sep=",", na="-9999.0", row.names=F, eol="\n")
close(f.out)
