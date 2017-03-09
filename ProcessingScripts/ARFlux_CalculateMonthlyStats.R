## ARFlux_CalculateMonthlyStats.R
#' This script is intended to calculate monthly statistics describing 
#' meteorological variables from subdaily data.

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

require(dplyr)
require(lubridate)

# read in data (use all sites)
df.in <- read.csv(paste0(git.dir, "data/ARFlux/ARFlux-Merged_2008-2012_Daily.csv"), stringsAsFactors=F)

# make month column
df.in$Date <- ymd(df.in$Date)
df.in$Year <- year(df.in$Date)
df.in$Month <- month(df.in$Date)

# calculate temperature range
df.in$Tair.C.range <- df.in$Tair.C.max - df.in$Tair.C.min

# trim df.in to relevant columns
df.in <- df.in[,c("fire", "Date", "Year", "Month", "Tair.C.range", "windSpeed.m.s", "windSpeed.m.s.std")]

# aggregate to daily values for each fire
df.mo <- summarize(group_by(df.in, Month),
                   Tair.C.range = mean(Tair.C.range, na.rm=T))#,
#                   windSpeed.m.s = mean(windSpeed.m.s, na.rm=T),
#                   windSpeed.m.s.std = mean(windSpeed.m.s.std, na.rm=T))

# save output
write.csv(df.mo, paste0(git.dir, "data/ARFlux/ARFlux-Merged_MonthlyStats.csv"), row.names=F)
