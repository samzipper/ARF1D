## ARFlux_HalfHourlyToDaily.R
# This aggregates half-hourly ARFlux data to daily means/sums.

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
  
  # summarize all to daily values
  df.d <- summarize(group_by(df, Date),
                    fire = fire,
                    P.kPa = mean(Pressure, na.rm=T),
                    Tair.C.mean = mean(Air.Temperature, na.rm=T),
                    Tair.C.min = min(Air.Temperature, na.rm=T),
                    Tair.C.max = max(Air.Temperature, na.rm=T),
                    Tsurf.K = mean(Surface.Temperature..Apogee., na.rm=T),
                    Tsoil.C = mean(Soil.Temperature, na.rm=T),
                    Tsoil.C.min = min(Soil.Temperature, na.rm=T),
                    Tsoil.C.max = max(Soil.Temperature, na.rm=T),
                    ea.kPa = mean(Ambient.Vapor.Pressure, na.rm=T),
                    SWin.W.m2 = mean(Incoming.Shortwave, na.rm=T),
                    SWout.W.m2 = mean(Outgoing.Shortwave, na.rm=T),
                    LWin.mm.m2.s = mean(Incoming.Longwave, na.rm=T),
                    LWout.mm.m2.s = mean(Outgoing.Longwave, na.rm=T),
                    Rnet.W.m2 = mean(Net.Radiation),
                    frictionVeloc.m2.s = mean(Friction.Velocity, na.rm=T),
                    windDir.deg = mean(Wind.Direction, na.rm=T),
                    windSpeed.m2.s = mean(Wind.Speed, na.rm=T),
                    LE.W.m2 = mean(Latent.Heat.Flux),
                    H.W.m2 = mean(Sensible.Heat.Flux),
                    G.W.m2 = mean(Ground.Heat.Flux),
                    NEE.mm.m2.s = mean(Net.Ecosystem.Exchange.of.CO2))

  # make some plots
  p.Tair.C <- 
    ggplot(df.d, aes(x=Date, y=Tair.C.mean)) +
    geom_hline(yintercept=0, color="gray65") +
    geom_point() +
    scale_y_continuous(limits=c(-35,20), breaks=seq(-30,20,10)) +
    theme_bw() +
    theme(panel.grid=element_blank())
  
  p.Tsoil.C <- 
    ggplot(df.d, aes(x=Date, y=Tsoil.C)) +
    geom_hline(yintercept=0, color="gray65") +
    geom_point() +
    scale_y_continuous(limits=c(-35,20), breaks=seq(-30,20,10)) +
    theme_bw() +
    theme(panel.grid=element_blank())
  
  ggsave(paste0("p.Tair.Tsoil_", fire, ".png"),
         arrangeGrob(p.Tair.C, p.Tsoil.C, ncol=1),
         width=12, height=8, units="in")
  
  df.energy <- df.d[,c("Date", "Rnet.W.m2", "LE.W.m2", "H.W.m2", "G.W.m2")]
  
  if (exists("df.d.all")){
    df.d.all <- rbind(df.d.all, df.d)
  } else {
    df.d.all <- df.d
  }
  
}

# save as CSV
write.csv(df.d.all, "ARFlux_2008-2012_Daily.csv", row.names=F)