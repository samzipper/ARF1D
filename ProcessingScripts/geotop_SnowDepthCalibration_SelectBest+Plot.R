## geotop_SnowDepthCalibration_SelectBest+Plot_Uniorm.R
#' This is intended to go through the output of a bunch of different 
#' simulations, calculate some fit metrics for each, and choose the
#' parameters with the best fit. A plot will be made of the best fit
#' only. 

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

require(lubridate)
require(ggplot2)
require(reshape2)
require(dplyr)
require(tidyr)
require(gridExtra)
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))

# path where GEOtop results are stored, and output should be saved (including plot)
path.out <- paste0(git.dir, "geotop/SnowDepthCalibration/")

# path of validation data
path.snow <- paste0(git.dir, "data/meteo/Raw/TFS-EDC_3-hour_data.csv")

# Process observation data ------------------------------------------------

# load observations
df.snow <- read.csv(path.snow, stringsAsFactors=F)

# format date
df.snow$Date <- ymd(df.snow$date)

# aggregate snow to daily
df.snow <- summarize(group_by(df.snow, Date),
                     snow_depth.obs.mm = mean(snow_depth, na.rm=T)*10)  # convert to mm
df.snow <- df.snow[year(df.snow$Date)>=2014, ]  # snow depth data only available 2014-2015
df.snow$snow_depth.obs.mm[df.snow$snow_depth.obs.mm < 0] <- 0

# Process model data ------------------------------------------------------

# load file that summarize all scenarios
df.in <- read.csv(paste0(path.out, "SnowDepthCalibration_Input.csv"), stringsAsFactors=F)

# change "number" column to character with padded zeros
df.in$number <- sprintf("%04d", df.in$number)

# make empty columns for fit metrics
df.in$RMSE <- NaN
df.in$NRMSE <- NaN
df.in$NSE <- NaN
df.in$R2 <- NaN

# scroll through each scenario
for (i in 1:dim(df.in)[1]){
  # path to point model output file
  path.mod.point <- paste0(path.out, df.in$number[i], "_point0001.txt")
  
  # read in data
  df.mod.point <- read.csv(path.mod.point, stringsAsFactors=F)
  
  # convert date
  df.mod.point$Date <- as.Date(dmy_hm(df.mod.point$Date12.DDMMYYYYhhmm.))
  
  # merge model with observations
  df.snow.i <- merge(df.snow, df.mod.point[,c("Date", "snow_depth.mm.")], all.x=T)
  colnames(df.snow.i)[colnames(df.snow.i)=="snow_depth.mm."] <- "snow_depth.mod.mm"
  
  # calculate fit metrics
  df.in$RMSE[i] <- RMSE(df.snow.i$snow_depth.mod.mm[complete.cases(df.snow.i)], df.snow.i$snow_depth.obs.mm[complete.cases(df.snow.i)])
  df.in$NRMSE[i] <- NRMSE(df.snow.i$snow_depth.mod.mm[complete.cases(df.snow.i)], df.snow.i$snow_depth.obs.mm[complete.cases(df.snow.i)])
  df.in$NSE[i] <- NashSutcliffe(df.snow.i$snow_depth.mod.mm[complete.cases(df.snow.i)], df.snow.i$snow_depth.obs.mm[complete.cases(df.snow.i)])
  df.in$R2[i] <- R2(df.snow.i$snow_depth.mod.mm[complete.cases(df.snow.i)], df.snow.i$snow_depth.obs.mm[complete.cases(df.snow.i)])

  # status update
  print(paste0(df.in$number[i], " complete"))
}

# save results
write.csv(df.in, paste0(path.out, "SnowDepthCalibration_Output.csv"), row.names=F)

# Determine best fit and make plot ----------------------------------------

# figure out best based on minimum RMSE
i.best <- which.min(df.in$RMSE)
n.best <- df.in$number[i.best]

# read in best
df.mod.best <- read.csv(paste0(path.out, n.best, "_point0001.txt"), stringsAsFactors=F)

# convert date
df.mod.best$Date <- as.Date(dmy_hm(df.mod.best$Date12.DDMMYYYYhhmm.))

# merge model with observations
df.snow.best <- merge(df.snow, df.mod.best[,c("Date", "snow_depth.mm.")], all.x=T)
colnames(df.snow.best)[colnames(df.snow.best)=="snow_depth.mm."] <- "snow_depth.mod.mm"

# make plot
p.snow.depth.cal <- 
  ggplot(df.snow.best) +
  geom_ribbon(fill="deepskyblue1", aes(x=Date, ymin=0, ymax=snow_depth.mod.mm)) +
  geom_point(aes(x=Date, y=snow_depth.obs.mm), shape=21) +
  labs(title=paste0("Snow Depth Calibration ", n.best), subtitle=paste0("RMSE=", round(df.in$RMSE[i.best],3), " (", 100*round(df.in$NRMSE[i.best],3), "%), NSE=", round(df.in$NSE[i.best], 3), ", R2=", round(df.in$R2[i.best], 3))) +
  scale_y_continuous(name="Snow Depth [mm]", expand=c(0,0), 
                     limits=c(0,max(c(max(df.snow.best$snow_depth.obs.mm), max(df.snow.best$snow_depth.mod.mm))*1.05))) +
  scale_x_date(expand=c(0,0), date_labels="%Y-%m", date_breaks="3 months", date_minor_breaks="1 month") +
  theme_bw() +
  theme(panel.grid=element_blank())
ggsave(paste0(path.out, "p.snow.depth.cal.png"), p.snow.depth.cal,
       width=6, height=4, units="in")