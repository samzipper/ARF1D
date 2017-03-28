## geotop_PostProcessAndPlot.R
#' This is intended to compare modeled and measured data:
#'   -Soil temperature
#'   -Soil moisture

rm(list=ls())

# git directory for relative paths
#git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"
git.dir <- "C:/Users/Sam/WorkGits/ARF1D/"

require(lubridate)
require(ggplot2)
require(reshape2)
require(dplyr)
require(gridExtra)
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))

#version
version <- "20170328-1hr-Start1988-RootGrow"

# function to find closest
which.closest <- function(x, vec){
  # function returns the index (i) of the element in the vector (vec) which is closest to the input value (x).
  # if x is equidistant from multiple elements in x, all are returned.
  # a tolerance parameter is built in to avoid rounding errors.
  # 
  eps <- .Machine$double.eps^0.5  # tolerance parameter
  
  i <- which(abs(vec-x)-min(abs(vec-x), na.rm=T) < eps)
  
  return(i)
}

## calibrate with odd years, validate with even
cal <- seq(1998,2003,1)
val <- seq(2005,2011,1)

yr.min <- min(c(cal, val))
yr.max <- max(c(cal, val))

# path to save plots
path.plot <- paste0(git.dir, "geotop_NRCS/output-plots/")

## read in observed data and preprocess
# paths
path.obs <- paste0(git.dir, "data/NRCS/Toolik_AllYears_Clean.csv")

# read in data
df.obs <- read.csv(path.obs, stringsAsFactors=F)

colnames(df.obs)[colnames(df.obs)=="date"] <- "Date"

# add date column
df.obs$Date <- ymd(df.obs$Date)

# subset to only calibration/validation period
df.obs <- subset(df.obs, year(Date)>=yr.min & year(Date)<=yr.max)

# get rid of 2004 - bad observation data
df.obs <- subset(df.obs, year(Date) != 2004)

## Process model output data
# modeled data path
path.mod.temp <- paste0(git.dir, "geotop_NRCS/output-tabs/soiltemp0001.txt")
path.mod.VWC <- paste0(git.dir, "geotop_NRCS/output-tabs/thetaliq0001.txt")

# read in data
df.mod.temp <- read.csv(path.mod.temp, stringsAsFactors=F)
df.mod.VWC <- read.csv(path.mod.VWC, stringsAsFactors=F)

# make Date column
df.mod.temp$Date <- dmy_hm(df.mod.temp$Date12.DDMMYYYYhhmm.)
df.mod.VWC$Date <- dmy_hm(df.mod.VWC$Date12.DDMMYYYYhhmm.)

# subset
df.mod.temp <- subset(df.mod.temp, year(Date)>=yr.min & year(Date)<=yr.max)
df.mod.VWC <- subset(df.mod.VWC, year(Date)>=yr.min & year(Date)<=yr.max)

# format Date column
df.mod.temp$Date <- as.Date(df.mod.temp$Date)
df.mod.VWC$Date <- as.Date(df.mod.VWC$Date)

# for profile files, keep column for Date and anything beginning with "X" (these are the depths)
all.cols <- colnames(df.mod.temp)
cols.keep <- c("Date", all.cols[startsWith(all.cols, "X")])
df.mod.temp <- df.mod.temp[,cols.keep]
df.mod.VWC <- df.mod.VWC[,cols.keep]

# reformat as long-form
df.mod.temp.melt <- melt(df.mod.temp, id="Date", value.name="value.mod")
df.mod.temp.melt$depth.mm <- as.numeric(sub("X", "", df.mod.temp.melt$variable))

df.mod.VWC.melt <- melt(df.mod.VWC, id="Date", value.name="value.mod")
df.mod.VWC.melt$depth.mm <- as.numeric(sub("X", "", df.mod.VWC.melt$variable))

# add column names to match df.obs
df.mod.temp.melt$variable <- "temp"
df.mod.VWC.melt$variable <- "VWC"

# combine
df.mod.melt <- rbind(df.mod.temp.melt, df.mod.VWC.melt)

# merge with obs data frame
df.all <- merge(df.obs, df.mod.melt, by=c("Date", "variable", "depth.mm"), all.y=T)

# reorder
df.all <- df.all[order(df.all$variable, df.all$depth.mm, df.all$Date),]

# make calibration/validation column
df.all$period <- "cal"
df.all$period[year(df.all$Date) %in% val] <- "val"

# depths to plot for temp and VWC
depths.temp <- unique(df.all$depth.mm[is.finite(df.all$value.mean) & df.all$variable=="temp"])
depths.VWC <- unique(df.all$depth.mm[is.finite(df.all$value.mean) & df.all$variable=="VWC"])

# calculate fit metrics
df.fit <- summarize(group_by(subset(df.all, is.finite(value.mean)), variable, depth.mm, period),
                    RMSE = RMSE(value.mod, value.mean),
                    NRMSE = NRMSE(value.mod, value.mean),
                    R2 = R2(value.mod, value.mean),
                    NSE = NashSutcliffe(value.mod, value.mean))
df.fit.melt <- melt(df.fit, id=c("variable", "depth.mm", "period"), value.name="fit", variable.name="metric")

# plot
p.temp.all <-
  ggplot(subset(df.all, variable=="temp" & depth.mm %in% depths.temp), aes(x=Date)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(aes(y=value.mod), color="black") +
  geom_line(aes(y=value.mean), color="red") +
  facet_wrap(~depth.mm) +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())
ggsave(paste0(path.plot, "p.temp.all_", version, ".png"), p.temp.all, width=12, height=6, units="in")

p.VWC.all <-
  ggplot(subset(df.all, variable=="VWC" & depth.mm %in% depths.VWC), aes(x=Date)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(aes(y=value.mod), color="black") +
  geom_line(aes(y=value.mean), color="blue") +
  facet_wrap(~depth.mm) +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())
ggsave(paste0(path.plot, "p.VWC.all_", version, ".png"), p.VWC.all, width=8, height=6, units="in")

# plot of only common depths
p.temp.VWC <-
  ggplot(subset(df.all, depth.mm %in% depths.VWC), aes(x=Date)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_line(aes(y=value.mod), color="black") +
  geom_line(aes(y=value.mean), color="forestgreen") +
  facet_grid(variable~depth.mm, scales="free") +
  scale_x_date(expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())
ggsave(paste0(path.plot, "p.temp.VWC_", version, ".png"), p.temp.VWC, width=10, height=6, units="in")

# plot fit
p.fit <-
  ggplot(df.fit.melt, aes(x=factor(depth.mm), y=fit, fill=period)) +
  geom_hline(yintercept=0, color="gray65") +
  geom_bar(stat="identity", position="dodge") +
  facet_grid(metric~variable, scales="free") +
  scale_x_discrete(name="Depth [mm]") +
  scale_y_continuous(name="Fit") +
  scale_fill_discrete(name="Period") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")
ggsave(paste0(path.plot, "p.fit_", version, ".png"), p.fit, width=8, height=6, units="in")
