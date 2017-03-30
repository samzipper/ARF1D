## geotop_output_WaterBalance.R
#' This is intended to check the water balance of a GEOtop simulation.

rm(list=ls())

# git directory for relative paths
#git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"
git.dir <- "C:/Users/Sam/WorkGits/ARF1D/"

# which geotop version
geo.dir <- "geotop_NRCS/"
#geo.dir <- "geotop/"

require(lubridate)
require(ggplot2)
require(reshape2)
require(dplyr)
require(gridExtra)
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))

#version
version <- "20170327-1hr"

## Process model output data
# modeled data path
path.mod.point <- paste0(git.dir, geo.dir, "output-tabs/point0001.txt")
path.mod.discharge <- paste0(git.dir, geo.dir, "output-tabs/discharge.txt")
path.mod.VWC <- paste0(git.dir, geo.dir, "output-tabs/thetaliq0001.txt")
path.mod.ice <- paste0(git.dir, geo.dir, "output-tabs/thetaice0001.txt")
path.soil <- paste0(git.dir, geo.dir, "soil/soilNRCS0001.txt")

# read in data
df.mod.point <- read.csv(path.mod.point, stringsAsFactors=F)
df.mod.discharge <- read.csv(path.mod.discharge, stringsAsFactors=F)
df.mod.VWC <- read.csv(path.mod.VWC, stringsAsFactors=F)
df.mod.ice <- read.csv(path.mod.ice, stringsAsFactors=F)
df.soil <- read.csv(path.soil, stringsAsFactors=F)

# make Date column
df.mod.point$Date <- dmy_hm(df.mod.point$Date12.DDMMYYYYhhmm.)
df.mod.discharge$Date <- dmy_hm(df.mod.discharge$DATE.day.month.year.hour.min.)
df.mod.VWC$Date <- dmy_hm(df.mod.VWC$Date12.DDMMYYYYhhmm.)
df.mod.ice$Date <- dmy_hm(df.mod.ice$Date12.DDMMYYYYhhmm.)

# for profile files, keep column for Date and anything beginning with "X" (these are the depths)
all.cols <- colnames(df.mod.VWC)
cols.keep <- c("Date", all.cols[startsWith(all.cols, "X")])
df.mod.VWC <- df.mod.VWC[,cols.keep]
df.mod.ice <- df.mod.ice[,cols.keep]

# combine VWC and ice
df.mod.water <- df.mod.VWC[,-1]+df.mod.ice[,-1]
df.mod.water$Date <- df.mod.VWC$Date
df.mod.water <- df.mod.water[,cols.keep]

# calculate total storage and storage change
storage.mm <- apply(apply(df.mod.water[,-1], 1, function(x) x*df.soil$Dz), 2, sum)
storage.change.mm <- diff(storage.mm, lag=1)

# calculate water balance components
df.mod.point$ET.mm <- df.mod.point$Evap_surface.mm. + df.mod.point$Trasp_canopy.mm.
df.mod.point$precip.mm <- df.mod.point$Psnow_over_canopy.mm. + df.mod.point$Prain_over_canopy.mm.
df.mod.point$precip.effective.mm <- df.mod.point$Prain_under_canopy.mm. + df.mod.point$snow_melted.mm.

# merge discharge out bottom
df.mod.point$drainage.mm <- df.mod.discharge$Qoutbottom.m3.s.*(24*60*60)*1000  # convert to mm/day

# merge lateral outflow
df.mod.point$lateral.mm <- df.mod.discharge$Qoutlandsub.m3.s.*(24*60*60)*1000  # convert to mm/day

# merge storage change
df.mod.point$storage.change.mm <- storage.change.mm

# calculate storage change based on inflows/outflows
df.mod.point$storage.change.calc.mm <- df.mod.point$precip.effective.mm - df.mod.point$ET.mm - df.mod.point$drainage.mm - df.mod.point$lateral.mm

# summarize to annual water balance
df.mod.point$year <- year(df.mod.point$Date)
df.mod.point.yr <- dplyr::summarize(group_by(df.mod.point, year),
                                    ET.mm = sum(ET.mm),
                                    precip.mm = sum(precip.mm),
                                    precip.effective.mm = sum(precip.effective.mm),
                                    drainage.mm = sum(drainage.mm),
                                    lateral.mm = sum(lateral.mm),
                                    storage.change.mm = sum(storage.change.mm),
                                    storage.change.calc.mm = sum(storage.change.calc.mm))

df.mod.point.yr$storage.change.calc.mm <- df.mod.point.yr$precip.effective.mm - df.mod.point.yr$ET.mm - df.mod.point.yr$drainage.mm

# plots
p.stor.day <-
  ggplot(df.mod.point, aes(x=storage.change.mm, y=storage.change.calc.mm)) +
  geom_point()

p.stor.yr <-
  ggplot(df.mod.point.yr, aes(x=storage.change.mm, y=storage.change.calc.mm)) +
  geom_point()

p.stor.yr.ts <-
  ggplot(df.mod.point.yr, aes(x=year)) +
  geom_line(aes(y=storage.change.mm), color="black") +
  geom_line(aes(y=storage.change.calc.mm), color="red")

# components from each year
df.mod.point.yr.melt <- melt(df.mod.point.yr[,c("year", "ET.mm", "precip.effective.mm", "drainage.mm", "lateral.mm", "storage.change.mm")], id="year")

p.yr <- 
  ggplot(df.mod.point.yr.melt, aes(x=year, y=value, color=variable)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line() +
  geom_point() +
  ggtitle("Annual Water Balance for 10 Years") +
  scale_y_continuous(name="Depth [mm]") +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(name="", labels=c("storage.change.mm"="Soil Storage\nChange",
                                       "precip.effective.mm"="Effective Precip\n(Rain+Snowmelt)",
                                       "drainage.mm"="Drainage",
                                       "lateral.mm"="Lateral\nOutflow",
                                       "ET.mm"="ET"), 
                     values=c("storage.change.mm"="brown",
                              "precip.effective.mm"="blue",
                              "drainage.mm"="orange",
                              "lateral.mm"="purple",
                              "ET.mm"="forestgreen")) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")

# components from a single event
# single rainfall followed by dry days
#df.event <- melt(df.mod.point[1329:1341,c("Date", "storage.change.mm", "storage.change.calc.mm","precip.effective.mm","drainage.mm","ET.mm")], id=c("Date"))

# many rainy days in a row
df.event <- melt(df.mod.point[1248:1268,c("Date", "storage.change.mm", "precip.effective.mm","drainage.mm","lateral.mm","ET.mm")], id=c("Date"))

df.event.sum <- summarize(group_by(df.event, variable),
                          total = sum(value))
p.event <- 
  ggplot(df.event, aes(x=Date, y=value, color=variable)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line() +
  geom_point() +
  ggtitle("Daily Water Balance for a Snowmelt + Rain Event") +
  scale_y_continuous(name="Depth [mm]") +
  scale_x_datetime(expand=c(0,0)) +
  scale_color_manual(name="", labels=c("storage.change.mm"="Soil Storage\nChange",
                                       "precip.effective.mm"="Effective Precip\n(Rain+Snowmelt)",
                                       "drainage.mm"="Drainage",
                                       "lateral.mm"="Lateral\nOutflow",
                                       "ET.mm"="ET"), 
                     values=c("storage.change.mm"="brown",
                              "precip.effective.mm"="blue",
                              "drainage.mm"="orange",
                              "lateral.mm"="purple",
                              "ET.mm"="forestgreen")) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")

ggsave(paste0(git.dir, geo.dir, "output-plots/WaterBudget.png"),
       arrangeGrob(p.yr, p.event, ncol=1), width=6, height=8, units="in")
