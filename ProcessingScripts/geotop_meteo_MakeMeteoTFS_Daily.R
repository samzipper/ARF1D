## geotop_meteo_MakeMeteoTFS_Daily.R
#' This script is intended to take a gap-filled meteorological
#' file and convert it to the necessary format for GEOtop.

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

require(lubridate)

es <- function(Tair.C){
  0.6108*exp((Tair.C*17.27)/(Tair.C+237.3))
}

# filename of baseline met data
fname.in <- paste0(git.dir, "data/meteo/Daily_1989-2013_TFS_GapFill.csv")
fname.out <- paste0(git.dir, "geotop/meteo/meteoTFSdaily0001.txt")

# read in met data
df.in <- read.csv(fname.in)
df.in$Date <- ymd(df.in$Date)

# calculate RH
df.in$RH <- 100*df.in$ea.kPa/((es(df.in$Tair.C.min)+es(df.in$Tair.C.max))/2)
df.in$RH[df.in$RH>100] <- 100
df.in$RH[df.in$RH<0] <- 0

# assemble output data frame
df.out <- data.frame(POSIX=format(df.in$Date, "%d/%m/%Y %H:%M"),
                     Iprec = df.in$precip.mm/24,  # convert to mm/hr
                     WindSp = df.in$wind.m_s,
                     AirT = df.in$Tair.C.mean,
                     RH = df.in$RH,
                     P = df.in$P.kPa*10,  # convert kPa to mbar
                     Swglob = df.in$rad.W_m2)

# write output
write.table(df.out, file=fname.out, quote=F, sep=",", na="-9999.0", row.names=F)
