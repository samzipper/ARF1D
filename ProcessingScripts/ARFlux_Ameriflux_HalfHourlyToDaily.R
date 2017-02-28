## ARFlux_Ameriflux_HalfHourlyToDaily.R
# This aggregates half-hourly ARFlux Ameriflux data to daily means/sums.

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
  # define path based on site
  if (fire=="Unburned"){
    path.fire <- "Raw-Ameriflux/AMF_US-An3_BASE_HH_1-1.csv"
  } else if (fire=="Moderate"){
    path.fire <- "Raw-Ameriflux/AMF_US-An2_BASE_HH_1-1.csv"
  } else if (fire=="Severe"){
    path.fire <- "Raw-Ameriflux/AMF_US-An1_BASE_HH_1-1.csv"
  }
  
  # read in CSV
  df <- read.csv(path.fire, skip=2, stringsAsFactors=F)
  
  # set nodata
  df[df==-9999] <- NaN
  
  # make a date column (this is the end of the half-hourly measurement period)
  df$Date <- ymd_hm(df$TIMESTAMP_END)
  
  # make date column
  df$Year <- year(df$Date)
  df$DOY <- yday(df$Date)
  
  # figure out timestep of data
  ts <- minutes(df$Date[2]-df$Date[1])
  max.gap.pts <- hours(max.gap.hrs)/ts
  
  # gap-fill columns that you want to output, using the same names as the 
  # data downloaded from ARC LTER so you don't have to mess with rest of script
  df$Pressure <- na.approx(df$PA, na.rm=F, maxgap=max.gap.pts)
  df$Air.Temperature <- na.approx(df$TA, na.rm=F, maxgap=max.gap.pts)
  df$Soil.Temperature <- na.approx(df$TS_1, na.rm=F, maxgap=max.gap.pts)
  df$SWC <- na.approx(df$SWC_1, na.rm=F, maxgap=max.gap.pts)
  df$RH <- na.approx(df$RH, na.rm=F, maxgap=max.gap.pts)
  df$Incoming.Shortwave <- na.approx(df$SW_IN, na.rm=F, maxgap=max.gap.pts)
  df$Outgoing.Shortwave <- na.approx(df$SW_OUT, na.rm=F, maxgap=max.gap.pts)
  df$Incoming.Longwave <- na.approx(df$LW_IN, na.rm=F, maxgap=max.gap.pts)
  df$Outgoing.Longwave <- na.approx(df$LW_OUT, na.rm=F, maxgap=max.gap.pts)
  df$Net.Radiation <- na.approx(df$NETRAD, na.rm=F, maxgap=max.gap.pts)
  df$Friction.Velocity <- na.approx(df$USTAR, na.rm=F, maxgap=max.gap.pts)
  df$Wind.Direction <- na.approx(df$WD, na.rm=F, maxgap=max.gap.pts)
  df$Wind.Speed <- na.approx(df$WS, na.rm=F, maxgap=max.gap.pts)
  df$Latent.Heat.Flux <- na.approx(df$LE, na.rm=F, maxgap=max.gap.pts)
  df$Sensible.Heat.Flux <- na.approx(df$H, na.rm=F, maxgap=max.gap.pts)
  df$Ground.Heat.Flux <- na.approx(df$G, na.rm=F, maxgap=max.gap.pts)
  df$Net.Ecosystem.Exchange.of.CO2 <- na.approx(df$NEE_PI, na.rm=F, maxgap=max.gap.pts)
  
  # calculate derived variables
  df$Saturation.Vapor.Pressure <- es(df$Air.Temperature)
  df$Ambient.Vapor.Pressure <- (df$RH/100)*df$Saturation.Vapor.Pressure
  
  # summarize all to daily values
  df.d <- summarize(group_by(df, Year, DOY),
                    fire = fire,
                    P.kPa = mean(Pressure),
                    Tair.C.mean = mean(Air.Temperature),
                    Tair.C.min = min(Air.Temperature),
                    Tair.C.max = max(Air.Temperature),
                    Tsoil.C = mean(Soil.Temperature),
                    Tsoil.C.min = min(Soil.Temperature),
                    Tsoil.C.max = max(Soil.Temperature),
                    VWC = mean(SWC/100),
                    VWC.min = min(SWC/100),
                    VWC.max = max(SWC/100),
                    RH = mean(RH),
                    es.kPa = mean(Saturation.Vapor.Pressure),
                    ea.kPa = mean(Ambient.Vapor.Pressure),
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

  # add a date column
  df.d$Date <- ymd(paste0(df.d$Year, "-01-01"))+(df.d$DOY-1)
  
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
  
  p.VWC <- 
    ggplot(df.d, aes(x=Date, y=VWC)) +
    #geom_hline(yintercept=0, color="gray65") +
    geom_point() +
    scale_y_continuous() +
    theme_bw() +
    theme(panel.grid=element_blank())
  
  ggsave(paste0("ARFlux-Ameriflux_p.Tair.Tsoil.VWC_", fire, ".png"),
         arrangeGrob(p.Tair.C, p.Tsoil.C, p.VWC, ncol=1),
         width=12, height=8, units="in")
  
  if (exists("df.d.all")){
    df.d.all <- rbind(df.d.all, df.d)
  } else {
    df.d.all <- df.d
  }
  
}

# save as CSV
write.csv(df.d.all, "ARFlux-Ameriflux_2008-2010_Daily.csv", row.names=F)
