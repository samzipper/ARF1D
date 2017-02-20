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
version <- "20170220e-Snow"
fire <- "Unburned"  # options are: Unburned, Moderate, Severe

## paths
# modeled data path
path.mod.point <- paste0(git.dir, "geotop/output-tabs/point0001.txt")
path.snow <- paste0(git.dir, "data/meteo/Raw/TFS-EDC_3-hour_data.csv")

# save plot path
path.plot.snow <- paste0(git.dir, "geotop/output-plots/PostProcessAndPlot_Snow_", version, "-", fire, ".png")

## read in data
df.mod.point <- read.csv(path.mod.point, stringsAsFactors=F)
df.snow <- read.csv(path.snow, stringsAsFactors=F)

## make plots
# convert dates
df.mod.point$Date <- as.Date(dmy_hm(df.mod.point$Date12.DDMMYYYYhhmm.))
df.snow$Date <- ymd(df.snow$date)

# aggregate snow to daily
df.snow <- summarize(group_by(df.snow, Date),
                     snow_depth.obs.mm = mean(snow_depth, na.rm=T)*10)  # convert to mm
df.snow <- df.snow[year(df.snow$Date)>=2014, ]  # snow depth data only available 2014-2015
df.snow$snow_depth.obs.mm[df.snow$snow_depth.obs.mm < 0] <- 0
df.snow <- merge(df.snow, df.mod.point[,c("Date", "snow_depth.mm.")], all.x=T)
colnames(df.snow)[colnames(df.snow)=="snow_depth.mm."] <- "snow_depth.mod.mm"

## make plots
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