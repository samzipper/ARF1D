## geotop_output_CompareModelTemp_ARFlux.R
#' This is intended to compare modeled and measured temperature at different
#' depths using data from the ARFlux sensors themselves.
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
version <- "20170216-SoilThermalProps"
fire <- "Unburned"  # options are: Unburned, Moderate, Severe

## paths
# modeled data path
path.mod.temp <- paste0(git.dir, "geotop/output-tabs/soiltemp0001.txt")
path.obs <- paste0(git.dir, "data/ARFlux/ARFlux_2008-2012_Daily.csv")
path.thaw <- paste0(git.dir, "data/ARFlux/ARFlux_ThawDepths_2008-2014.csv")

# save plot path
path.plot <- paste0(git.dir, "geotop/output-plots/CompareModelTemp_ARFlux_", version, "-", fire, ".png")

## read in data
df.mod.temp <- read.csv(path.mod.temp, stringsAsFactors=F)
df.obs <- read.csv(path.obs, stringsAsFactors=F)
df.thaw <- read.csv(path.thaw, stringsAsFactors=F)

# trim observations to only unburned
df.obs <- df.obs[df.obs$fire==fire,]
df.thaw <- data.frame(Date = dmy(df.thaw$Date),
                      ThawDepth.mm = 10*as.numeric(df.thaw[,which(colnames(df.thaw)==paste0("Thaw.depth.", fire))]))

## pre-processing of observed data
# assign depths based on site 
# (depths are listed on repository, e.g. http://arc-lter.ecosystems.mbl.edu/2012arfluxunburned )
depths <- c(20, 60)  # depth is average of 2 cm and 6 cm sensor

# calculate mean depth
depth.min <- min(depths)
depth.mean <- mean(depths)
depth.max <- max(depths)

# trim data frame
df.obs <- df.obs[,c("Date", "Tsoil.C", "Tsoil.C.min", "Tsoil.C.max")]

## Process model output data
# make Date column
df.mod.temp$Date <- format(dmy_hm(df.mod.temp[,1]), "%m/%d/%Y")

# keep column for Date and anything beginning with "X" (these are the depths)
all.cols <- colnames(df.mod.temp)
cols.keep <- c("Date", all.cols[startsWith(all.cols, "X")])
df.mod.temp <- df.mod.temp[,cols.keep]

# figure out depth for each column
cols.depth <- as.numeric(sub("X", "", cols.keep))

# determine which columns are within the depth range of the observations
cols.compare <- which(cols.depth > depth.min & cols.depth < depth.max)

# summarize model for each day
Dates.all <- unique(df.mod.temp$Date)
df.mod.day <- data.frame(Date = Dates.all,
                         Temp.mean = NaN,
                         Temp.std = NaN,
                         Temp.min = NaN,
                         Temp.max = NaN,
                         ThawDepth.mm = NaN)
for (d in Dates.all){
  # figure out index for this date
  i.d <- which(df.mod.day$Date==d)
  
  # summarize temperature
  df.mod.day$Temp.mean[i.d] <- mean(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
  df.mod.day$Temp.std[i.d] <- sd(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
  df.mod.day$Temp.min[i.d] <- min(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
  df.mod.day$Temp.max[i.d] <- max(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare]), na.rm=T)
  
  ## figure out thaw depth: shallowest depth at which temperature > 0
  # extract all temperatures for that date
  df.mod.temp.profile <- data.frame(depth = cols.depth[is.finite(cols.depth)], 
                                    Temp = unlist(df.mod.temp[df.mod.temp$Date==d, is.finite(cols.depth)]))
  
  # find the shallowest temperature that is frozen
  df.mod.day$ThawDepth.mm[i.d] <- min(df.mod.temp.profile$depth[which(df.mod.temp.profile$Temp<0)])
}

# convert dates
df.mod.day$Date <- mdy(df.mod.day$Date)
df.obs$Date <- ymd(df.obs$Date)

# trim model to 2010-2013
df.mod.day <- subset(df.mod.day, year(Date)>=2008 & year(Date)<=2012)
df.thaw <- subset(df.thaw, year(Date)>=2008 & year(Date)<=2012)

# for each thaw date, summarize to mean and range
df.thaw.day <- summarize(group_by(df.thaw, Date),
                         ThawDepth.mm.mean = mean(ThawDepth.mm, na.rm=T),
                         ThawDepth.mm.std = sd(ThawDepth.mm, na.rm=T),
                         ThawDepth.mm.min = min(ThawDepth.mm, na.rm=T),
                         ThawDepth.mm.max = max(ThawDepth.mm, na.rm=T))


## calculate fit metrics
# temperature: make a data frame that has a model value for each date there is an observed value
df.fit.temp <- df.obs[complete.cases(df.obs),]
df.fit.temp$Temp.mod <- df.mod.day$Temp.mean[match(df.fit.temp$Date, df.mod.day$Date)]
df.thaw.day$thaw.mod <- df.mod.day$ThawDepth.mm[match(df.thaw.day$Date, df.mod.day$Date)]

# make a matrix
df.fit.table <- data.frame(RMSE=numeric(2), NRMSE=numeric(2), NSE=numeric(2), R2=numeric(2),
                           row.names=c("SoilTemp", "ThawDepth"))
df.fit.table["SoilTemp",] <- c(RMSE(df.fit.temp$Temp.mod, df.fit.temp$Tsoil.C),
                               NRMSE(df.fit.temp$Temp.mod, df.fit.temp$Tsoil.C),
                               NashSutcliffe(df.fit.temp$Temp.mod, df.fit.temp$Tsoil.C),
                               R2(df.fit.temp$Temp.mod, df.fit.temp$Tsoil.C))
df.fit.table["ThawDepth",] <- c(RMSE(df.thaw.day$thaw.mod, df.thaw.day$ThawDepth.mm.mean),
                                NRMSE(df.thaw.day$thaw.mod, df.thaw.day$ThawDepth.mm.mean),
                                NashSutcliffe(df.thaw.day$thaw.mod, df.thaw.day$ThawDepth.mm.mean),
                                R2(df.thaw.day$thaw.mod, df.thaw.day$ThawDepth.mm.mean))

## make plots
# temperature plot
p.temp.compare <-
  ggplot() +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(data=df.mod.day, aes(x=Date, y=Temp.mean), color="red") +
  geom_point(data=df.obs, aes(x=Date, y=Tsoil.C)) +
  scale_x_date(expand=c(0,0)) +
  scale_y_continuous(name="Daily Soil Temperature [C]") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.thaw.compare <- 
  ggplot() +
  geom_line(data=df.mod.day, aes(x=Date, y=ThawDepth.mm), color="blue") +
  geom_segment(data=df.thaw.day, aes(x=Date, xend=Date, y=ThawDepth.mm.mean-ThawDepth.mm.std, yend=ThawDepth.mm.mean+ThawDepth.mm.std)) +
  geom_point(data=df.thaw.day, aes(x=Date, y=ThawDepth.mm.mean)) +
  scale_x_date(expand=c(0,0)) +
  scale_y_reverse(name="Thaw Depth [mm]") +
  theme_bw() +
  theme(panel.grid=element_blank())

ggsave(path.plot, arrangeGrob(p.temp.compare, p.thaw.compare, tableGrob(round(df.fit.table, 3)), 
                              ncol=1, heights=c(1,1,0.33)),
       width=12, height=8, units="in")
