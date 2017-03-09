## geotop_meteo_MakeMeteoSpinUpTFS+ARFlux_Daily.R
#' This script is intended to take a gap-filled meteorological
#' file, generate a long record by randomly selecting from existing
#' years, and convert it to the necessary format for GEOtop. A separate
#' file is made for each burn site, in which 2008-2013 air temperature 
#' data are taken from the relevant ARFlux site.

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

require(lubridate)

# saturation vapor pressure function
es <- function(Tair.C){
  0.6108*exp((Tair.C*17.27)/(Tair.C+237.3))
}

# filename of baseline met data
fname.in <- paste0(git.dir, "data/meteo/Daily_1989-2015_TFS_GapFill.csv")
fname.out <- paste0(git.dir, "geotop/meteo/meteoTFSdailyWithSpinUp")
fname.ARFlux <- paste0(git.dir, "data/ARFlux/ARFlux_2008-2012_Daily.csv")

# number of spin-up years
yr.n <- 100

# read in met data
df.in <- read.csv(fname.in)
df.in$Date <- ymd(df.in$Date)

# figure out start year
yr.start <- min(df.in$Year)-yr.n

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
write.table(df.out.combo, file=paste0(fname.out, "0001.txt"), quote=F, sep=",", na="-9999.0", row.names=F)

## now, make a separate file for each ARFlux site
# read in data
df.ARFlux <- read.csv(fname.ARFlux, stringsAsFactors=F)

# format date to match df.in.syn
df.ARFlux$Date <- format(ymd(df.ARFlux$Date), "%d/%m/%Y %H:%M")

# calculate relative humidity
df.ARFlux$RH <- 100*df.ARFlux$ea.kPa/es(df.ARFlux$Tair.C.mean)
df.ARFlux$RH[df.ARFlux$RH<0] <- 0
df.ARFlux$RH[df.ARFlux$RH>100] <- 100

# scroll through each burn site
for (fire in c("Unburned", "Moderate", "Severe")){
  # copy df.in.syn
  df.base <- df.out.combo
  
  # subset to that burn site only
  df.fire <- df.ARFlux[df.ARFlux$fire==fire, ]
  
  # find match for each df.ARFlux date
  i.fire <- match(df.fire$Date, df.base$POSIX)
  
  # replace each variable where it exists
  df.base$AirT[i.fire[is.finite(df.fire$Tair.C.mean)]] <- df.fire$Tair.C.mean[is.finite(df.fire$Tair.C.mean)]
  df.base$RH[i.fire[is.finite(df.fire$RH)]] <- df.fire$RH[is.finite(df.fire$RH)]
  df.base$WindSp[i.fire[is.finite(df.fire$windSpeed.m2.s)]] <- df.fire$windSpeed.m2.s[is.finite(df.fire$windSpeed.m2.s)]
  df.base$Swglob[i.fire[is.finite(df.fire$SWin.W.m2)]] <- df.fire$SWin.W.m2[is.finite(df.fire$SWin.W.m2)]
  df.base$P[i.fire[is.finite(df.fire$P.kPa)]] <- df.fire$P.kPa[is.finite(df.fire$P.kPa)]
  
  # write output
  write.table(df.base, file=paste0(fname.out, fire, "0001.txt"), quote=F, sep=",", na="-9999.0", row.names=F)
}

