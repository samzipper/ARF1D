## ARFlux_ThawDepth_Plot.R
# This script is intended to plot thaw depth data for ARFlux sites.

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/ARF1D/"

require(lubridate)
require(ggplot2)
require(dplyr)
require(reshape2)
require(stringr)
require(gridExtra)
require(zoo)

# load data
df <- read.csv(paste0(git.dir, "data/ARFlux/Raw/2008-2016ARF_ThawDepths.csv"), stringsAsFactors=F)

# convert date
df$date <- dmy(df$Date)
df$Date <- NULL

# melt
df.melt <- melt(df, id=c("date"), value.name="depth.cm")
df.melt$depth.cm <- as.numeric(df.melt$depth.cm)

# make fire.sev column
df.melt$fire.sev <- str_split_fixed(df.melt$variable, pattern="[.]", n=3)[,3]

# daily average
df.melt.d <- dplyr::summarize(group_by(df.melt, fire.sev, date),
                              thaw.depth.cm = mean(depth.cm, na.rm=T))
# plot
ggplot(df.melt.d, aes(x=date, y=thaw.depth.cm, color=fire.sev)) + geom_point() + scale_y_reverse()
