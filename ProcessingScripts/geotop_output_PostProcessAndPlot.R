## geotop_output_CompareModelTempVWC.R
#' This is intended to compare modeled and measured temperature at different
#' depths from a set of soil moisture/temperature sensors.
#' 

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

require(lubridate)
require(ggplot2)
require(reshape2)
require(dplyr)
require(gridExtra)
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))

# version name
version <- "20170220-TFSmet-FixCoord-FixDEM"
fire <- "Unburned"  # options are: Unburned, Moderate, Severe

## paths
# modeled data path
path.mod.temp <- paste0(git.dir, "geotop/output-tabs/soiltemp0001.txt")
path.mod.VWC <- paste0(git.dir, "geotop/output-tabs/thetaliq0001.txt")
path.mod.point <- paste0(git.dir, "geotop/output-tabs/point0001.txt")
path.obs <- paste0(git.dir, "data/ARFlux/ARFlux_SoilTemp-Moisture_", fire, "_2010-2013.csv")
path.obs.ARFlux <- paste0(git.dir, "data/ARFlux/ARFlux_2008-2012_Daily.csv")
path.thaw <- paste0(git.dir, "data/ARFlux/ARFlux_ThawDepths_2008-2014.csv")
path.snow <- paste0(git.dir, "data/meteo/Raw/TFS-EDC_3-hour_data.csv")

# save plot path
path.plot.val <- paste0(git.dir, "geotop/output-plots/PostProcessAndPlot_Val_", version, "-", fire, ".png")
path.plot.met <- paste0(git.dir, "geotop/output-plots/PostProcessAndPlot_Met_", version, "-", fire, ".png")
path.plot.snow <- paste0(git.dir, "geotop/output-plots/PostProcessAndPlot_Snow_", version, "-", fire, ".png")

## read in data
df.mod.temp <- read.csv(path.mod.temp, stringsAsFactors=F)
df.mod.VWC <- read.csv(path.mod.VWC, stringsAsFactors=F)
df.mod.point <- read.csv(path.mod.point, stringsAsFactors=F)
df.obs <- read.csv(path.obs, stringsAsFactors=F)
df.obs.ARFlux <- read.csv(path.obs.ARFlux, stringsAsFactors=F)
df.thaw <- read.csv(path.thaw, stringsAsFactors=F)
df.snow <- read.csv(path.snow, stringsAsFactors=F)

## pre-processing of observed data
# trim observations to only unburned
df.obs.ARFlux <- df.obs.ARFlux[df.obs.ARFlux$fire==fire,]
df.thaw <- data.frame(Date = dmy(df.thaw$Date),
                      ThawDepth.mm = 10*as.numeric(df.thaw[,which(colnames(df.thaw)==paste0("Thaw.depth.", fire))]))
# rename columns to standardized form
colnames(df.obs) <- c("Date", "Year", "Day", "Hour_Min", 
                      "Temp.1", "Temp.2", "Temp.3", "Temp.4", "Temp.5", "Temp.6",
                      "PA.1", "PA.2", "PA.3", "PA.4", "PA.5", "PA.6",
                      "VWC.1", "VWC.2", "VWC.3", "VWC.4", "VWC.5", "VWC.6")

# assign depths based on site 
# (depths are listed on repository, e.g. http://arc-lter.ecosystems.mbl.edu/arfunburned-soiltemp-moisture)
# note that depths may change between years and were not continuously monitored
# depth units are [mm]
if (fire=="Unburned"){
  depths <- c("1"=80, "2"=40, "3"=90, "4"=70, "5"=70, "6"=60)
} else if (fire=="Moderate"){
  depths <- c("1"=40, "2"=70, "3"=60, "4"=100, "5"=50, "6"=50)
} else if (fire=="Severe"){
  depths <- c("1"=70, "2"=70, "3"=50, "4"=55, "5"=90, "6"=45)
}

# calculate mean depth
depth.min <- min(depths)
depth.mean <- mean(depths)
depth.max <- max(depths)

# ARFlux depths are listed on repository, e.g. http://arc-lter.ecosystems.mbl.edu/2012arfluxunburned
depths.ARFlux <- c(20, 60)  # depth is average of 2 cm and 6 cm sensor

# calculate mean depth
depth.ARFlux.min <- min(depths.ARFlux)
depth.ARFlux.mean <- mean(depths.ARFlux)
depth.ARFlux.max <- max(depths.ARFlux)

# trim ARFlux data frame
df.obs.ARFlux <- df.obs.ARFlux[,c("Date", "Tsoil.C", "Tsoil.C.min", "Tsoil.C.max")]

## calculate mean temperature and VWC across all sensors for each day
VWC.cols <- c("VWC.1", "VWC.2", "VWC.3", "VWC.4", "VWC.5", "VWC.6")
Temp.cols <- c("Temp.1", "Temp.2", "Temp.3", "Temp.4", "Temp.5", "Temp.6")

# calculate just Date (without time)
df.obs$Date <- format(mdy_hm(df.obs$Date), "%m/%d/%Y")
Dates.all <- unique(df.obs$Date)
df.obs.day <- data.frame(Date = Dates.all, 
                         VWC.mean = NaN,
                         VWC.std = NaN,
                         VWC.min = NaN,
                         VWC.max = NaN,
                         Temp.mean = NaN,
                         Temp.std = NaN,
                         Temp.min = NaN,
                         Temp.max = NaN)
for (d in Dates.all){
  # figure out index for this date
  i.d <- which(df.obs.day$Date==d)
  
  # summarize VWC and temperature
  df.obs.day$VWC.mean[i.d] <- mean(unlist(df.obs[df.obs$Date==d, VWC.cols]), na.rm=T)/100
  df.obs.day$VWC.std[i.d] <- sd(unlist(df.obs[df.obs$Date==d, VWC.cols]), na.rm=T)/100
  df.obs.day$VWC.min[i.d] <- min(unlist(df.obs[df.obs$Date==d, VWC.cols]), na.rm=T)/100
  df.obs.day$VWC.max[i.d] <- max(unlist(df.obs[df.obs$Date==d, VWC.cols]), na.rm=T)/100
  
  df.obs.day$Temp.mean[i.d] <- mean(unlist(df.obs[df.obs$Date==d, Temp.cols]), na.rm=T)
  df.obs.day$Temp.std[i.d] <- sd(unlist(df.obs[df.obs$Date==d, Temp.cols]), na.rm=T)
  df.obs.day$Temp.min[i.d] <- min(unlist(df.obs[df.obs$Date==d, Temp.cols]), na.rm=T)
  df.obs.day$Temp.max[i.d] <- max(unlist(df.obs[df.obs$Date==d, Temp.cols]), na.rm=T)
}

## Process model output data
# make Date column
df.mod.temp$Date <- format(dmy_hm(df.mod.temp[,1]), "%m/%d/%Y")
df.mod.VWC$Date <- format(dmy_hm(df.mod.VWC[,1]), "%m/%d/%Y")
df.mod.point$Date <- format(dmy_hm(df.mod.point[,1]), "%m/%d/%Y")

# keep column for Date and anything beginning with "X" (these are the depths)
all.cols <- colnames(df.mod.temp)
cols.keep <- c("Date", all.cols[startsWith(all.cols, "X")])
df.mod.temp <- df.mod.temp[,cols.keep]
df.mod.VWC <- df.mod.VWC[,cols.keep]

# figure out depth for each column
cols.depth <- as.numeric(sub("X", "", cols.keep))

# determine which columns are within the depth range of the observations
cols.compare <- which(cols.depth > depth.min & cols.depth < depth.max)
cols.compare.ARFlux <- which(cols.depth > depth.ARFlux.min & cols.depth < depth.ARFlux.max)

# summarize model for each day
Dates.all <- unique(df.mod.temp$Date)
df.mod.day <- data.frame(Date = Dates.all, 
                         VWC.mean = NaN,
                         VWC.std = NaN,
                         VWC.min = NaN,
                         VWC.max = NaN,
                         Temp.mean = NaN,
                         Temp.std = NaN,
                         Temp.min = NaN,
                         Temp.max = NaN,
                         Temp.ARFlux.mean = NaN,
                         Temp.ARFlux.std = NaN,
                         Temp.ARFlux.min = NaN,
                         Temp.ARFlux.max = NaN,
                         ThawDepth.mm = NaN)
for (d in Dates.all){
  # figure out index for this date
  i.d <- which(df.mod.day$Date==d)
  
  # summarize VWC and temperature
  df.mod.day$VWC.mean[i.d] <- mean(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare]), na.rm=T)
  df.mod.day$VWC.std[i.d] <- sd(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare]), na.rm=T)
  df.mod.day$VWC.min[i.d] <- min(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare]), na.rm=T)
  df.mod.day$VWC.max[i.d] <- max(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare]), na.rm=T)
  
  df.mod.day$Temp.mean[i.d] <- mean(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
  df.mod.day$Temp.std[i.d] <- sd(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
  df.mod.day$Temp.min[i.d] <- min(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
  df.mod.day$Temp.max[i.d] <- max(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
  
  df.mod.day$Temp.ARFlux.mean[i.d] <- mean(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
  df.mod.day$Temp.ARFlux.std[i.d] <- sd(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
  df.mod.day$Temp.ARFlux.min[i.d] <- min(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
  df.mod.day$Temp.ARFlux.max[i.d] <- max(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
  
  ## figure out thaw depth: shallowest depth at which temperature > 0
  # extract all temperatures for that date
  df.mod.temp.profile <- data.frame(depth = cols.depth[is.finite(cols.depth)], 
                                    Temp = unlist(df.mod.temp[df.mod.temp$Date==d, is.finite(cols.depth)]))
  
  # find the shallowest temperature that is frozen
  df.mod.day$ThawDepth.mm[i.d] <- min(df.mod.temp.profile$depth[which(df.mod.temp.profile$Temp<0)])
}

## make plots
# convert dates
df.mod.day$Date <- mdy(df.mod.day$Date)
df.mod.point$Date <- mdy(df.mod.point$Date)
df.obs.day$Date <- mdy(df.obs.day$Date)
df.obs.ARFlux$Date <- ymd(df.obs.ARFlux$Date)
df.snow$Date <- ymd(df.snow$date)

# aggregate snow to daily
df.snow <- summarize(group_by(df.snow, Date),
                     snow_depth.obs.mm = mean(snow_depth, na.rm=T)*10)  # convert to mm
df.snow <- df.snow[year(df.snow$Date)>=2014, ]  # snow depth data only available 2014-2015
df.snow$snow_depth.obs.mm[df.snow$snow_depth.obs.mm < 0] <- 0
df.snow <- merge(df.snow, df.mod.point[,c("Date", "snow_depth.mm.")], all.x=T)
colnames(df.snow)[colnames(df.snow)=="snow_depth.mm."] <- "snow_depth.mod.mm"

# for each thaw date, summarize to mean and range
df.thaw.day <- summarize(group_by(df.thaw, Date),
                         ThawDepth.mm.mean = mean(ThawDepth.mm, na.rm=T),
                         ThawDepth.mm.std = sd(ThawDepth.mm, na.rm=T),
                         ThawDepth.mm.min = min(ThawDepth.mm, na.rm=T),
                         ThawDepth.mm.max = max(ThawDepth.mm, na.rm=T))


## calculate fit metrics
# ARFlux
df.fit.temp.ARFlux <- df.obs.ARFlux[complete.cases(df.obs.ARFlux),]
df.fit.temp.ARFlux$Temp.mod <- df.mod.day$Temp.ARFlux.mean[match(df.fit.temp.ARFlux$Date, df.mod.day$Date)]
df.thaw.day$thaw.mod <- df.mod.day$ThawDepth.mm[match(df.thaw.day$Date, df.mod.day$Date)]

# temp + VWC
df.fit.temp <- df.obs.day[is.finite(df.obs.day$Temp.mean),]
df.fit.temp$Temp.mod <- df.mod.day$Temp.mean[match(df.fit.temp$Date, df.mod.day$Date)]
df.fit.VWC <- df.obs.day[is.finite(df.obs.day$VWC.mean),]
df.fit.VWC$VWC.mod <- df.mod.day$VWC.mean[match(df.fit.VWC$Date, df.mod.day$Date)]

# make a matrix with fit statistics
df.fit.table <- data.frame(RMSE=numeric(4), NRMSE=numeric(4), NSE=numeric(4), R2=numeric(4),
                           row.names=c("Thaw Depth","Tower Soil Temp", "Nest Soil Temp", "Nest VWC"))
df.fit.table["Tower Soil Temp",] <- c(RMSE(df.fit.temp.ARFlux$Temp.mod, df.fit.temp.ARFlux$Tsoil.C),
                                      NRMSE(df.fit.temp.ARFlux$Temp.mod, df.fit.temp.ARFlux$Tsoil.C),
                                      NashSutcliffe(df.fit.temp.ARFlux$Temp.mod, df.fit.temp.ARFlux$Tsoil.C),
                                      R2(df.fit.temp.ARFlux$Temp.mod, df.fit.temp.ARFlux$Tsoil.C))
df.fit.table["Thaw Depth",] <- c(RMSE(df.thaw.day$thaw.mod, df.thaw.day$ThawDepth.mm.mean),
                                 NRMSE(df.thaw.day$thaw.mod, df.thaw.day$ThawDepth.mm.mean),
                                 NashSutcliffe(df.thaw.day$thaw.mod, df.thaw.day$ThawDepth.mm.mean),
                                 R2(df.thaw.day$thaw.mod, df.thaw.day$ThawDepth.mm.mean))
df.fit.table["Nest Soil Temp",] <- c(RMSE(df.fit.temp$Temp.mod, df.fit.temp$Temp.mean),
                                     NRMSE(df.fit.temp$Temp.mod, df.fit.temp$Temp.mean),
                                     NashSutcliffe(df.fit.temp$Temp.mod, df.fit.temp$Temp.mean),
                                     R2(df.fit.temp$Temp.mod, df.fit.temp$Temp.mean))
df.fit.table["Nest VWC",] <- c(RMSE(df.fit.VWC$VWC.mod, df.fit.VWC$VWC.mean),
                               NRMSE(df.fit.VWC$VWC.mod, df.fit.VWC$VWC.mean),
                               NashSutcliffe(df.fit.VWC$VWC.mod, df.fit.VWC$VWC.mean),
                               R2(df.fit.VWC$VWC.mod, df.fit.VWC$VWC.mean))


## make plots
# ARFlux temperature and thaw plot
p.thaw.compare.ARFlux <- 
  ggplot() +
  geom_line(data=subset(df.mod.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=ThawDepth.mm), color="brown") +
  geom_segment(data=subset(df.thaw.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, xend=Date, y=ThawDepth.mm.mean-ThawDepth.mm.std, yend=ThawDepth.mm.mean+ThawDepth.mm.std)) +
  geom_point(data=subset(df.thaw.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=ThawDepth.mm.mean)) +
  scale_x_date(expand=c(0,0)) +
  scale_y_reverse(name="Tower Thaw Depth [mm]") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.temp.compare.ARFlux <-
  ggplot() +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(data=subset(df.mod.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=Temp.mean), color="red") +
  geom_point(data=df.obs.ARFlux, aes(x=Date, y=Tsoil.C)) +
  scale_x_date(expand=c(0,0)) +
  scale_y_continuous(name=paste0("Tower Soil Temp [C]: ", depth.ARFlux.min, "-", depth.ARFlux.max, " mm")) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.temp.compare <-
  ggplot() +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(data=subset(df.mod.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=Temp.mean), color="red") +
  geom_point(data=df.obs.day, aes(x=Date, y=Temp.mean)) +
  scale_y_continuous(name=paste0("Nest Soil Temp [C]: ", depth.min, "-", depth.max, " mm")) +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.VWC.compare <-
  ggplot() +
  geom_line(data=subset(df.mod.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=VWC.mean), color="blue") +
  geom_point(data=df.obs.day, aes(x=Date, y=VWC.mean)) +
  scale_y_continuous(name=paste0("Nest VWC [m3/m3]: ", depth.min, "-", depth.max, " mm")) +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

ggsave(path.plot.val, arrangeGrob(p.thaw.compare.ARFlux, p.temp.compare.ARFlux, p.temp.compare, p.VWC.compare, tableGrob(round(df.fit.table, 3)), 
                                  ncol=1, heights=c(1,1,1,1,0.5)),
       width=12, height=12, units="in")

p.precip.mm <- 
  ggplot(subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, xend=Date, y=0, yend=Prain_over_canopy.mm.+Psnow_over_canopy.mm.)) +
  geom_segment(color="blue") +
  scale_y_continuous(name="Precipitation (Rain + Snow) [mm]", expand=c(0,0), 
                     limits=c(0,max(subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013)$Prain_over_canopy.mm.+subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013)$Psnow_over_canopy.mm.)*1.05)) +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.snow.depth.mm <- 
  ggplot(subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, ymin=0, ymax=snow_depth.mm.)) +
  geom_ribbon(fill="deepskyblue1") +
  scale_y_continuous(name="Snow Depth [mm]", expand=c(0,0), 
                     limits=c(0,max(subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013)$snow_depth.mm.)*1.05)) +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.Tair.C <- 
  ggplot(subset(df.mod.point, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=Tair.C.)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(color="red") +
  scale_y_continuous(name="Air Temperature [C]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

ggsave(path.plot.met, arrangeGrob(p.Tair.C, p.precip.mm, p.snow.depth.mm, p.temp.compare.ARFlux, 
                                  ncol=1, heights=c(1,1,1,1)),
       width=12, height=12, units="in")

# snow validation fit metrics & plot
snow.RMSE <- RMSE(df.snow$snow_depth.mod.mm[complete.cases(df.snow)], df.snow$snow_depth.obs.mm[complete.cases(df.snow)])
snow.NRMSE <- NRMSE(df.snow$snow_depth.mod.mm[complete.cases(df.snow)], df.snow$snow_depth.obs.mm[complete.cases(df.snow)])
snow.NSE <- NashSutcliffe(df.snow$snow_depth.mod.mm[complete.cases(df.snow)], df.snow$snow_depth.obs.mm[complete.cases(df.snow)])
snow.R2 <- R2(df.snow$snow_depth.mod.mm[complete.cases(df.snow)], df.snow$snow_depth.obs.mm[complete.cases(df.snow)])
p.snow.depth.mm.val <- 
  ggplot(df.snow) +
  geom_ribbon(fill="deepskyblue1", aes(x=Date, ymin=0, ymax=snow_depth.mod.mm)) +
  geom_point(aes(x=Date, y=snow_depth.obs.mm)) +
  ggtitle(paste0("RMSE=", round(snow.RMSE,3), " (", 100*round(snow.NRMSE,3), "%), NSE=", round(snow.NSE, 3), ", R2=", round(snow.R2, 3))) +
  scale_y_continuous(name="Snow Depth [mm]", expand=c(0,0), 
                     limits=c(0,max(c(max(df.snow$snow_depth.obs.mm), max(df.snow$snow_depth.mod.mm))*1.05))) +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.snow.precip.mm <- 
  ggplot(subset(df.mod.point, year(Date)>=2014 & year(Date)<=2015), aes(x=Date, xend=Date, y=0, yend=Prain_over_canopy.mm.+Psnow_over_canopy.mm.)) +
  geom_segment(color="blue") +
  scale_y_continuous(name="Precipitation (Rain + Snow) [mm]", expand=c(0,0), 
                     limits=c(0,max(subset(df.mod.point, 
                                           year(Date)>=2014 & year(Date)<=2014)$Prain_over_canopy.mm.+
                                      subset(df.mod.point, year(Date)>=2014 & year(Date)<=2014)$Psnow_over_canopy.mm.)*1.05)) +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())

p.snow.Tair.C <- 
  ggplot(subset(df.mod.point, year(Date)>=2014 & year(Date)<=2015), aes(x=Date, y=Tair.C.)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(color="red") +
  scale_y_continuous(name="Air Temperature [C]") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())
ggsave(path.plot.snow, arrangeGrob(p.snow.depth.mm.val, p.snow.precip.mm, p.snow.Tair.C, 
                                  ncol=1, heights=c(1,1,1)),
       width=12, height=9, units="in")