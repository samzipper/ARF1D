## ARFlux_MergeLTER+Ameriflux.R
#' This script is intended to take output from ARFlux data downloaded from 
#' the ARC LTER portal (ARFlux_2008-2012_Daily.csv) and from the Ameriflux
#' portal (ARFlux-Ameriflux_2008-2010_Daily.csv) and merge them into a single
#' CSV file (ARFlux-Merged_2008-2012_Daily.csv). While the two datasets are 
#' basically the same, apart from some small differences which seem ~rounding
#' errors, the Ameriflux dataset has soil moisture while the LTER dataset does
#' not. However, the Ameriflux dataset only goes through the year 2010, while 
#' the LTER dataset includes 2011 and 2012. For the merged product, Ameriflux 
#' VWC data are merged into the LTER data frame where they exist.

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

require(ggplot2)
require(gridExtra)
require(lubridate)

setwd(paste0(git.dir, "data/ARFlux"))

# open LTER and Ameriflux files, which were already gap-filled and aggregated to daily
df.LTER <- read.csv("ARFlux_2008-2012_Daily.csv", stringsAsFactors=F)
df.Amer <- read.csv("ARFlux-Ameriflux_2008-2010_Daily.csv", stringsAsFactors=F)

# mess with date columns
df.LTER$Date <- ymd(df.LTER$Date)
df.Amer$Date <- ymd(df.Amer$Date)

# figure out which colnames are not present in the other
cols.missing <- colnames(df.Amer)[which(!colnames(df.Amer) %in% colnames(df.LTER) & !colnames(df.Amer) %in% c("Year", "DOY"))]

# merge Amer into LTER
df.merge <- merge(df.LTER, df.Amer[, c("Date", "fire", cols.missing)], all.x=T, by=c("Date", "fire"))

# write output
write.csv(df.merge, "ARFlux-Merged_2008-2012_Daily.csv", row.names=F)
