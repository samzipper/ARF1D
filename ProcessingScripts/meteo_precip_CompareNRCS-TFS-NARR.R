## meteo_precip_CompareNRCS-NARR-TFS.R
#' This script is intended to compare precipitation data from 
#' the three available data sources: NRCS, NARR, and TFS.

rm(list=ls())

# git directory for relative paths
#git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"
git.dir <- "C:/Users/Sam/WorkGits/ARF1D/"

require(lubridate)
require(dplyr)
require(ggplot2)

# paths to daily files
path.NRCS <- paste0(git.dir, "data/NRCS/Toolik_AllYears.csv")
path.NARR <- paste0(git.dir, "data/meteo/Raw/NARR_Daily_1979-2015.csv")
path.TFS <- paste0(git.dir, "data/meteo/Daily_1989-2015_TFS_GapFill.csv")

# load files
df.NRCS <- read.csv(path.NRCS, stringsAsFactors=F, skip=3)
df.NARR <- read.csv(path.NARR, stringsAsFactors=F)
df.TFS <- read.csv(path.TFS, stringsAsFactors=F)

# trim to only date and precip
df.NRCS <- df.NRCS[,c("year", "DOY", "precip.in")]
df.NARR <- df.NARR[,c("year", "DOY", "precip.mm")]
df.TFS <- df.TFS[,c("Year", "DOY", "precip.mm")]

# conversions
df.NRCS$precip.mm <- df.NRCS$precip.in*25.4
df.NRCS$precip.in <- NULL

colnames(df.NRCS) <- c("year", "DOY", "precip.mm.NRCS")
colnames(df.NARR) <- c("year", "DOY", "precip.mm.NARR")
colnames(df.TFS) <- c("year", "DOY", "precip.mm.TFS")

# summarize to annual
df.NRCS.yr <- summarize(group_by(df.NRCS, year),
                        precip.mm.NRCS = sum(precip.mm.NRCS))
df.NARR.yr <- summarize(group_by(df.NARR, year),
                        precip.mm.NARR = sum(precip.mm.NARR))
df.TFS.yr <- summarize(group_by(df.TFS, year),
                        precip.mm.TFS = sum(precip.mm.TFS))

# get rid of 2009 TFS - weird data
df.TFS.yr <- subset(df.TFS.yr, year!=2009)

# bind NARR and TFS
df.yr <- merge(df.NARR.yr, df.TFS.yr, by="year")

# plot
p.NARR.TFS <-
  ggplot(df.yr, aes(x=precip.mm.NARR, y=precip.mm.TFS)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point() +
  stat_smooth(method="lm")