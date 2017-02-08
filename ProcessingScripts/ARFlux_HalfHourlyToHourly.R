## ARFlux_HalfHourlyToHourly.R
# This aggregates half-hourly ARFlux data to hourly means/sums.

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

require(lubridate)
require(ggplot2)
require(dplyr)
require(reshape2)
require(gridExtra)
source(paste0(git.dir, "ProcessingScripts/FindMissingDates.R"))

setwd(paste0(git.dir, "data/ARFlux"))

for (fire in c("Moderate", "Severe", "Unburned")){
  # read in CSV
  df <- read.csv(paste0("Raw/ARFlux_2008-2012_", fire, ".csv"))
  
  # make date column
  df$DOY.dec <- df$DOY
  df$Date <- as.Date(ymd(paste0(df$Year, "-01-01")) + minutes(round((df$DOY.dec*24*60))))
  df$DOY <- yday(df$Date)
  df$hour <- floor(round((1+df$DOY.dec - df$DOY)*24,1))
  
  # summarize all to daily values
  df.h <- summarize(group_by(df, Year, DOY, hour),
                    fire = fire,
                    P.kPa = mean(Pressure, na.rm=T),
                    Tair.C = mean(Air.Temperature, na.rm=T),
                    Tsurf.K = mean(Surface.Temperature..Apogee., na.rm=T),
                    Tsoil.C = mean(Soil.Temperature, na.rm=T),
                    Tsoil.C.min = min(Soil.Temperature, na.rm=T),
                    Tsoil.C.max = max(Soil.Temperature, na.rm=T),
                    ea.kPa = mean(Ambient.Vapor.Pressure, na.rm=T),
                    SWin.W.m2 = mean(Incoming.Shortwave, na.rm=T),
                    SWout.W.m2 = mean(Outgoing.Shortwave, na.rm=T),
                    LWin.mm.m2.s = mean(Incoming.Longwave, na.rm=T),
                    LWout.mm.m2.s = mean(Outgoing.Longwave, na.rm=T),
                    Rnet.W.m2 = mean(Net.Radiation, na.rm=T),
                    frictionVeloc.m2.s = mean(Friction.Velocity, na.rm=T),
                    windDir.deg = mean(Wind.Direction, na.rm=T),
                    windSpeed.m2.s = mean(Wind.Speed, na.rm=T),
                    LE.W.m2 = mean(Latent.Heat.Flux, na.rm=T),
                    H.W.m2 = mean(Sensible.Heat.Flux, na.rm=T),
                    G.W.m2 = mean(Ground.Heat.Flux, na.rm=T),
                    NEE.mm.m2.s = mean(Net.Ecosystem.Exchange.of.CO2, na.rm=T))
  
  if (exists("df.h.all")){
    df.h.all <- rbind(df.h.all, df.h)
  } else {
    df.h.all <- df.h
  }
  
}

# save as CSV
write.csv(df.h.all, "ARFlux_2008-2012_Hourly.csv", row.names=F)
