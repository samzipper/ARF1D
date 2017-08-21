## meteo_TFS-EDC_CompareToTFS-LTER.R
#' This compares data downloaded from the LTER database (TFS-LTER) to 
#' data downloaded from the TFS-EDC website.

rm(list=ls())

# git directory for relative paths
#git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"
git.dir <- "C:/Users/Sam/WorkGits/ARF1D/"

require(lubridate)
require(ggplot2)

## load data
df.LTER <- read.csv(paste0(git.dir, "data/meteo/Daily_1989-2015_TFS_GapFill.csv"), stringsAsFactors=F)
df.EDC <- read.csv(paste0(git.dir, "data/meteo/TFS-EDC_1988-2016_Daily.csv"), stringsAsFactors=F)

# make date and DOY column
df.LTER$Date <- ymd(df.LTER$Date)
df.EDC$date <- ymd(df.EDC$date)
colnames(df.EDC) <- c("Date", "EDC_Tair.C.mean", "EDC_Tair.C.min", "EDC_Tair.C.max", "EDC_RH", "EDC_P.kPa", "EDC_wind.m_s", "EDC_precip.mm", "EDC_rad.W_m2")

# join
df <- left_join(df.LTER, df.EDC)

# plot
ggplot(df, aes(x=Tair.C.mean, y=EDC_Tair.C.mean)) + geom_point()
ggplot(df, aes(x=precip.mm, y=EDC_precip.mm)) + geom_point()
ggplot(subset(df, Year>=2012), aes(x=precip.mm, y=EDC_precip.mm)) + geom_point()
ggplot(subset(df, Year>=2012), aes(x=Tair.C.mean, y=EDC_Tair.C.mean)) + geom_point()

