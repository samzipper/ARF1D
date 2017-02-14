## meteo_MakeGEOtopSpinUp.R
#' This script is intended to take a gap-filled meteorological
#' file, generate a long record by randomly selecting from existing
#' years, and convert it to the necessary format for GEOtop.

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

require(lubridate)

es <- function(Tair.C){
  0.6108*exp((Tair.C*17.27)/(Tair.C+237.3))
}

# filename of baseline met data
fname.in <- paste0(git.dir, "data/meteo/Daily_1989-2013_TFS_GapFill.csv")
fname.out <- paste0(git.dir, "geotop/meteo/meteoTFSdailyWithSpinUp0001.txt")

# number of spin-up years
yr.n <- 100

# read in met data
df.in <- read.csv(fname.in)
df.in$Date <- ymd(df.in$Date)

# figure out start year
yr.start <- min(df.in$Year)-yr.n

# calculate RH
df.in$RH <- 100*df.in$ea.kPa/((es(df.in$Tair.C.min)+es(df.in$Tair.C.max))/2)
df.in$RH[df.in$RH>100] <- 100
df.in$RH[df.in$RH<0] <- 0

# list of years that exist
yrs.exist.leap <- unique(df.in$Year[leap_year(df.in$Year)])
yrs.exist.nonleap <- unique(df.in$Year[!leap_year(df.in$Year)])

# cycle through years and build output dataset
set.seed(1)
for (yr in yr.start:(yr.start+yr.n-1)){
  # select synthetic year, depending on whether or not it is a leap year
  if (leap_year(yr)){
    yr.syn <- base::sample(yrs.exist.leap,1)
  } else {
    yr.syn <- base::sample(yrs.exist.nonleap,1)
  }
  
  # assemble output data frame for that year
  df.out.syn <- data.frame(POSIX=format((ymd(paste0(yr, "-01-01")) + days(subset(df.in, Year==yr.syn)$DOY - 1)), "%d/%m/%Y %H:%M"),
                           Iprec = subset(df.in, Year==yr.syn)$precip.mm/24,  # convert to mm/hr
                           WindSp = subset(df.in, Year==yr.syn)$wind.m_s,
                           AirT = subset(df.in, Year==yr.syn)$Tair.C.mean,
                           RH = subset(df.in, Year==yr.syn)$RH,
                           P = subset(df.in, Year==yr.syn)$P.kPa*10,  # convert kPa to mbar
                           Swglob = subset(df.in, Year==yr.syn)$rad.W_m2)
  
  if (yr==yr.start){
    df.out <- df.out.syn
  } else {
    df.out <- rbind(df.out, df.out.syn)
  }
  
}

# now, get df.in into the same format
df.in.syn <- data.frame(POSIX=format(df.in$Date, "%d/%m/%Y %H:%M"),
                        Iprec = df.in$precip.mm/24,
                        WindSp = df.in$wind.m_s,
                        AirT = df.in$Tair.C.mean,
                        RH = df.in$RH,
                        P = df.in$P.kPa*10,
                        Swglob = df.in$rad.W_m2)

# combin df.in and df.out
df.out.combo <- rbind(df.out, df.in.syn)

# write output
write.table(df.out.combo, file=fname.out, quote=F, sep=",", na="-9999.0", row.names=F)
