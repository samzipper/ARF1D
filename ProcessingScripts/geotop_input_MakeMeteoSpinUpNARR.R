## geotop_input_MakeMeteoSpinUpNARR.R
#' This script is intended to extract data for a given coordinate from NARR
#' NetCDF files generate a long-term record, and then build spin-up data by
#' randomly selecting years.
#' 
#' NARR link: https://www.esrl.noaa.gov/psd/data/gridded/data.narr.monolevel.html

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"
NARR.dir <- "D:/Dropbox/GIS_GeneralFiles/Meteorology_Gridded/NARR/NARR-DAILIES-MONOLEVEL/"

require(lubridate)
require(ncdf4)

# filename to save output
fname.NARR.meteo.out <- paste0(git.dir, "data/meteo/Raw/NARR_Daily_1979-2015.csv")
fname.out <- paste0(git.dir, "geotop/meteo/meteoNARRdailyWithSpinUp0001.txt")

# number of spin-up years
yr.n <- 100

# coordinates of site
site.lat <- 68.93      # unburned=68.93, moderate=68.95, severe=68.99
site.lon <- -150.27    # unburned=-150.27, moderate=-150.21, severe=-150.28

# First, build NARR input dataset for 1979-2015 ---------------------------
yr.list <- seq(1979,2015)

# get some metadata that will be needed for all stations and shouldn't change between data products or years
NARR.nc <- nc_open(paste0(NARR.dir, "air.2m.1979.nc"))
NARR.lat <- ncvar_get(NARR.nc, "lat")    # latitude [degN]
NARR.lon <- ncvar_get(NARR.nc, "lon")    # longitude [degE]
nc_close(NARR.nc)

# find index of your point based on difference between lat/long
# get total lat and lon difference and find index of this point
site.diff <- abs(NARR.lat-site.lat)+abs(NARR.lon-site.lon)
i.NARR <- which(site.diff==min(site.diff), arr.ind=T)

# make empty data.frame
df.out  <- data.frame(year = numeric(0),
                      DOY = numeric(0),
                      Tmean.C = numeric(0),     # mean air temperature at 2 m [K]
                      Tdew.C = numeric(0),      # mean dewpoint temperature at 2 m [K]
                      precip.mm = numeric(0),   # total precipitation [mm]
                      SWin = numeric(0),        # mean downward shortwave radiation [W/m2]
                      SWout = numeric(0),       # mean upward shortwave radiation [W/m2]
                      LWin = numeric(0),        # mean downward longwave radiation [W/m2]
                      LWout = numeric(0),       # mean upward longwave radiation [W/m2]
                      P = numeric(0),           # mean air pressure [Pa]
                      RH = numeric(0),          # mean relative humidity at 2 m [%]
                      wind.u = numeric(0),      # mean u-component of wind at 10 m [m/s]
                      wind.v = numeric(0))      # mean v-component of wind at 10 m [m/s]

# scroll through years and grab data
for (yr in yr.list){

  # determine number of days in year
  yr.days <- yday(ymd(paste0(yr, "-12-31")))

  # open NetCDFs for this year
  air.2m.nc <- nc_open(paste0(NARR.dir, "air.2m.", yr, ".nc"))
  apcp.nc <- nc_open(paste0(NARR.dir, "apcp.", yr, ".nc"))
  dlwrf.nc <- nc_open(paste0(NARR.dir, "dlwrf.", yr, ".nc"))
  ulwrf.sfc.nc <- nc_open(paste0(NARR.dir, "ulwrf.sfc.", yr, ".nc"))
  dswrf.nc <- nc_open(paste0(NARR.dir, "dswrf.", yr, ".nc"))
  uswrf.sfc.nc <- nc_open(paste0(NARR.dir, "uswrf.sfc.", yr, ".nc"))
  dpt.2m.nc <- nc_open(paste0(NARR.dir, "dpt.2m.", yr, ".nc"))
  pres.sfc.nc <- nc_open(paste0(NARR.dir, "pres.sfc.", yr, ".nc"))
  rhum.2m.nc <- nc_open(paste0(NARR.dir, "rhum.2m.", yr, ".nc"))
  uwnd.10m.nc <- nc_open(paste0(NARR.dir, "uwnd.10m.", yr, ".nc"))
  vwnd.10m.nc <- nc_open(paste0(NARR.dir, "vwnd.10m.", yr, ".nc"))

  # grab data for this year
  df.year <- data.frame(year = yr,
                        DOY = seq(1,yr.days),
                        Tmean.C = ncvar_get(air.2m.nc, start=c(i.NARR, 1), count=c(1, 1, yr.days))-273.15,
                        Tdew.C = ncvar_get(dpt.2m.nc, start=c(i.NARR, 1), count=c(1, 1, yr.days))-273.15,
                        precip.mm = ncvar_get(apcp.nc, start=c(i.NARR, 1), count=c(1, 1, yr.days)),
                        SWin = ncvar_get(dswrf.nc, start=c(i.NARR, 1), count=c(1, 1, yr.days)),
                        SWout = ncvar_get(uswrf.sfc.nc, start=c(i.NARR, 1), count=c(1, 1, yr.days)),
                        LWin = ncvar_get(dlwrf.nc, start=c(i.NARR, 1), count=c(1, 1, yr.days)),
                        LWout = ncvar_get(ulwrf.sfc.nc, start=c(i.NARR, 1), count=c(1, 1, yr.days)),
                        P = ncvar_get(pres.sfc.nc, start=c(i.NARR, 1), count=c(1, 1, yr.days)),
                        RH = ncvar_get(rhum.2m.nc, start=c(i.NARR, 1), count=c(1, 1, yr.days)),
                        wind.u = ncvar_get(uwnd.10m.nc, start=c(i.NARR, 1), count=c(1, 1, yr.days)),
                        wind.v = ncvar_get(vwnd.10m.nc, start=c(i.NARR, 1), count=c(1, 1, yr.days)))

  # bind to overall data frame
  df.out <- rbind(df.out, df.year)

  # close NetCDFs
  nc_close(air.2m.nc)
  nc_close(apcp.nc)
  nc_close(dlwrf.nc)
  nc_close(ulwrf.sfc.nc)
  nc_close(dswrf.nc)
  nc_close(uswrf.sfc.nc)
  nc_close(dpt.2m.nc)
  nc_close(pres.sfc.nc)
  nc_close(rhum.2m.nc)
  nc_close(uwnd.10m.nc)
  nc_close(vwnd.10m.nc)

  # status update
  print(paste0(yr, " complete"))
}

# save as CSV
write.csv(df.out, fname.NARR.meteo.out, row.names=F, quote=F)

# Build data file formatted to GEOtop standards with spin-up --------------

# read in NARR data
df.in <- read.csv(fname.NARR.meteo.out, stringsAsFactors=F)
colnames(df.in)[colnames(df.in)=="year"] <- "Year"
df.in$Date <- ymd(paste0(df.in$Year, "-01-01"))+days(df.in$DOY-1)

# figure out start year
yr.start <- min(df.in$Year)-yr.n

# list of years that exist
yrs.exist.leap <- unique(df.in$Year[leap_year(df.in$Year)])
yrs.exist.nonleap <- unique(df.in$Year[!leap_year(df.in$Year)])

# calculate wind speed in m/s
df.in$wind.m_s <- ((df.in$wind.u^2)+(df.in$wind.v^2))^0.5

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
                           AirT = subset(df.in, Year==yr.syn)$Tmean.C,
                           RH = subset(df.in, Year==yr.syn)$RH,
                           P = subset(df.in, Year==yr.syn)$P/1000,  # convert Pa to kPa
                           Swglob = subset(df.in, Year==yr.syn)$SWin)
  
  if (yr==yr.start){
    df.out <- df.out.syn
  } else {
    df.out <- rbind(df.out, df.out.syn)
  }
  
  # status update
  print(paste0(yr, " complete"))
  
}

# now, get df.in into the same format
df.in.syn <- data.frame(POSIX=format(df.in$Date, "%d/%m/%Y %H:%M"),
                        Iprec = df.in$precip.mm/24,
                        WindSp = df.in$wind.m_s,
                        AirT = df.in$Tmean.C,
                        RH = df.in$RH,
                        P = df.in$P/1000,
                        Swglob = df.in$SWin)

# combin df.in and df.out
df.out.combo <- rbind(df.out, df.in.syn)

# order
df.out.combo <- df.out.combo[order(dmy_hm(df.out.combo$POSIX)), ]

# write output
write.table(df.out.combo, file=fname.out, quote=F, sep=",", na="-9999.0", row.names=F)
