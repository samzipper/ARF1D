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
require(zoo)
source(paste0(git.dir, "ProcessingScripts/FindMissingDates.R"))
es <- function(Tair.C){
  0.6108*exp((Tair.C*17.27)/(Tair.C+237.3))
}

setwd(paste0(git.dir, "data/ARFlux"))

# gap-filling: max gap to use linear interpolation
max.gap.hrs <- 6 # [hrs]

for (fire in c("Moderate", "Severe", "Unburned")){
  # read in CSV
  df <- read.csv(paste0("Raw/ARFlux_2008-2012_", fire, ".csv"))
  
  # make date column
  df$DOY.dec <- df$DOY
  df$Date <- as.Date(ymd(paste0(df$Year, "-01-01")) + minutes(round((df$DOY.dec*24*60))))
  df$DOY <- yday(df$Date)
  
  # figure out timestep of data
  ts <- minutes(round((df$DOY.dec[2]-df$DOY.dec[1])*24*60))
  max.gap.pts <- hours(max.gap.hrs)/ts
  
  # gap-fill columns that you want to output
  df$Pressure <- na.approx(df$Pressure, na.rm=F, maxgap=max.gap.pts)
  df$Air.Temperature <- na.approx(df$Air.Temperature, na.rm=F, maxgap=max.gap.pts)
  df$Soil.Temperature <- na.approx(df$Soil.Temperature, na.rm=F, maxgap=max.gap.pts)
  df$Surface.Temperature..Apogee. <- na.approx(df$Surface.Temperature..Apogee., na.rm=F, maxgap=max.gap.pts)
  df$Ambient.Vapor.Pressure <- na.approx(df$Ambient.Vapor.Pressure, na.rm=F, maxgap=max.gap.pts)
  df$Incoming.Shortwave <- na.approx(df$Incoming.Shortwave, na.rm=F, maxgap=max.gap.pts)
  df$Outgoing.Shortwave <- na.approx(df$Outgoing.Shortwave, na.rm=F, maxgap=max.gap.pts)
  df$Incoming.Longwave <- na.approx(df$Incoming.Longwave, na.rm=F, maxgap=max.gap.pts)
  df$Outgoing.Longwave <- na.approx(df$Outgoing.Longwave, na.rm=F, maxgap=max.gap.pts)
  df$Net.Radiation <- na.approx(df$Net.Radiation, na.rm=F, maxgap=max.gap.pts)
  df$Friction.Velocity <- na.approx(df$Friction.Velocity, na.rm=F, maxgap=max.gap.pts)
  df$Wind.Direction <- na.approx(df$Wind.Direction, na.rm=F, maxgap=max.gap.pts)
  df$Wind.Speed <- na.approx(df$Wind.Speed, na.rm=F, maxgap=max.gap.pts)
  df$Latent.Heat.Flux <- na.approx(df$Latent.Heat.Flux, na.rm=F, maxgap=max.gap.pts)
  df$Sensible.Heat.Flux <- na.approx(df$Sensible.Heat.Flux, na.rm=F, maxgap=max.gap.pts)
  df$Ground.Heat.Flux <- na.approx(df$Ground.Heat.Flux, na.rm=F, maxgap=max.gap.pts)
  df$Net.Ecosystem.Exchange.of.CO2 <- na.approx(df$Net.Ecosystem.Exchange.of.CO2, na.rm=F, maxgap=max.gap.pts)
  
  # calculate derived variables
  df$Saturation.Vapor.Pressure <- es(df$Air.Temperature)
  df$RH <- 100*df$Ambient.Vapor.Pressure/df$Saturation.Vapor.Pressure
  
  # summarize all to daily values
  df.d <- summarize(group_by(df, Date),
                    fire = fire,
                    P.kPa = mean(Pressure),
                    Tair.C.mean = mean(Air.Temperature),
                    Tair.C.min = min(Air.Temperature),
                    Tair.C.max = max(Air.Temperature),
                    Tsurf.K = mean(Surface.Temperature..Apogee.),
                    Tsoil.C = mean(Soil.Temperature),
                    Tsoil.C.min = min(Soil.Temperature),
                    Tsoil.C.max = max(Soil.Temperature),
                    ea.kPa = mean(Ambient.Vapor.Pressure),
                    es.kPa = mean(Saturation.Vapor.Pressure),
                    RH = mean(RH),
                    SWin.W.m2 = mean(Incoming.Shortwave),
                    SWout.W.m2 = mean(Outgoing.Shortwave),
                    LWin.mm.m2.s = mean(Incoming.Longwave),
                    LWout.mm.m2.s = mean(Outgoing.Longwave),
                    Rnet.W.m2 = mean(Net.Radiation),
                    frictionVeloc.m.s = mean(Friction.Velocity),
                    windDir.deg = mean(Wind.Direction),
                    windSpeed.m.s = mean(Wind.Speed),
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
  
  ggsave(paste0("ARFlux_p.Tair.Tsoil_", fire, ".png"),
         arrangeGrob(p.Tair.C, p.Tsoil.C, ncol=1),
         width=12, height=8, units="in")
  
  if (exists("df.d.all")){
    df.d.all <- rbind(df.d.all, df.d)
  } else {
    df.d.all <- df.d
  }
  
}

# save as CSV
write.csv(df.d.all, "ARFlux_2008-2012_Daily.csv", row.names=F)