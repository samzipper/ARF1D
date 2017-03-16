## geotop_input_MakeVegFile.R
#' This script is intended to make a time-dependent vegetation file.
#' A daily meteorological data file is used to determine the start/end
#' of the growing season, and also the timestep of the output file.

rm(list=ls())

# git directory for relative paths
#git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"
git.dir <- "C:/Users/Sam/WorkGits/ARF1D/"

require(lubridate)
require(dplyr)

# some vegetation parameters to set
VegHeight.min <- 200  # [mm] - min veg height (50 mm is min allowed in GEOtop)
VegHeight.max <- 500  # [mm] - max veg height (based on Google Image Search for Anaktuvuk River Fire, looked like it was ~knee high at flowering)
LSAI.min <- 0.05      # [m2/m2] - min LAI+SAI used in non-growing season
LSAI.max <- 2.0       # [m2/m2] - max LAI+SAI, based on Rocha & Shaver (2009) AFM
RootDepth.max <- 500  # [mm] - max root depth; defined based on ~annual thaw depth (Iversen et al., 2015, NP)
ThresVeg1 <- 50       # [mm] - at snow depths >ThresVeg1, vegetation completely buried (=ThresSnowVegUp)
ThresVeg2 <- 10       # [mm] - at snow depths >ThresVeg2, vegetation completely exposed (=ThresSnowVegDown)

# year to start/end data
yr.start <- 1999
yr.end <- 2015

# time to ramp up/down after SOS and before EOS
days.ramp <- 30

# filename to save output
fname.out <- paste0(git.dir, "geotop/veg/VegParams0001.txt")

# filename of meteo
fname.meteo <- paste0(git.dir, "geotop/meteo/meteoNARRdailyWithSpinUp0001.txt")

# read in meteo file
df.meteo <- read.csv(fname.meteo, stringsAsFactors=F)

# extract Date and Year
df.meteo$Date <- dmy_hm(df.meteo$POSIX)
df.meteo$year <- year(df.meteo$Date)
df.meteo$DOY <- yday(df.meteo$Date)

# subset to output years
df.meteo <- subset(df.meteo, year >= yr.start & year <= yr.end)

## find annual SOS and EOS based on when 10-day moving average air temperature exceeds 0 for first time (SOS) and last time (EOS)
# make 10-day moving average temperature
df.meteo$AirT.avg <- as.numeric(stats::filter(df.meteo$AirT, rep(0.1,10), sides=1))
df.meteo$GS <- df.meteo$AirT.avg>0  # growing season logical: T means during growing season

# for each year, find DOY of min and max GS
df.meteo.gs <- summarize(group_by(df.meteo[df.meteo$GS,], year),
                         SOS = min(DOY, na.rm=T),
                         EOS = max(DOY, na.rm=T))
df.meteo.gs <- df.meteo.gs[complete.cases(df.meteo.gs),]

# scroll through years and calculate time-dependent veg
for (yr in df.meteo.gs$year){
  # get SOS and EOS
  SOS <- df.meteo.gs$SOS[df.meteo.gs$year==yr]
  EOS <- df.meteo.gs$EOS[df.meteo.gs$year==yr]
  
  # make data frame
  df.yr <- data.frame(year = yr, 
                      DOY = df.meteo$DOY[df.meteo$year==yr],
                      VegHeight = VegHeight.min,
                      LSAI = LSAI.min,
                      CanopyFraction = 0.0,
                      RootDepth = RootDepth.max)
  
  # scale linearly between min and max for VegHeight and LSAI
  df.yr$VegHeight[df.yr$DOY>=SOS & df.yr$DOY<SOS+days.ramp] <- VegHeight.min + (VegHeight.max-VegHeight.min)*
    (df.yr$DOY[df.yr$DOY>=SOS & df.yr$DOY<SOS+days.ramp]-SOS)/days.ramp
  df.yr$VegHeight[df.yr$DOY>=SOS+days.ramp & df.yr$DOY<=EOS-days.ramp] <- VegHeight.max
  df.yr$VegHeight[df.yr$DOY>=(EOS-days.ramp) & df.yr$DOY<EOS] <- VegHeight.max - (VegHeight.max-VegHeight.min)*
    (df.yr$DOY[df.yr$DOY>=(EOS-days.ramp) & df.yr$DOY<EOS]-(EOS-days.ramp))/days.ramp
  
  df.yr$LSAI[df.yr$DOY>=SOS & df.yr$DOY<SOS+days.ramp] <- LSAI.min + (LSAI.max-LSAI.min)*
    (df.yr$DOY[df.yr$DOY>=SOS & df.yr$DOY<SOS+days.ramp]-SOS)/days.ramp
  df.yr$LSAI[df.yr$DOY>=SOS+days.ramp & df.yr$DOY<=EOS-days.ramp] <- LSAI.max
  df.yr$LSAI[df.yr$DOY>=(EOS-days.ramp) & df.yr$DOY<EOS] <- LSAI.max - (LSAI.max-LSAI.min)*
    (df.yr$DOY[df.yr$DOY>=(EOS-days.ramp) & df.yr$DOY<EOS]-(EOS-days.ramp))/days.ramp
  
# for canopy fraction, use Norman et al. (1995) AFM Eq. 3
df.yr$CanopyFraction <- 1-exp(-0.5*df.yr$LSAI)
  
  # make output data frame
  if (exists("df.yr.all")){
    df.yr.all <- rbind(df.yr.all, df.yr)
  } else {
    df.yr.all <- df.yr
  }
  
}

# make output data frame
df.out.combo <- data.frame(POSIX = df.meteo$POSIX,
                           VegHeight = df.yr.all$VegHeight,
                           ThresVeg1 = ThresVeg1,
                           ThresVeg2 = ThresVeg2,
                           LSAI = df.yr.all$LSAI,
                           CanopyFraction = df.yr.all$CanopyFraction,
                           DecayCoeff=2.5,
                           SnowBurialCoeff=1.0,
                           RootDepth = df.yr.all$RootDepth,
                           MinStomatalRes=30)

# write output
write.table(df.out.combo, file=fname.out, quote=F, sep=",", na="-9999.0", row.names=F)
