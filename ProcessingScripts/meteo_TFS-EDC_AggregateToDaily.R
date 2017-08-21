## meteo_TFS-EDC_AggregateToDaily.R
# This aggregates 1 hour Toolik Field Data Environmental Data Center
# meteorological data to daily means/sums.
#
# Raw data from here: http://toolik.alaska.edu/edc/abiotic_monitoring/data_query.php
# Conditions of use: http://toolik.alaska.edu/edc/about/conditions_of_use.php

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/ARF1D/"

require(lubridate)
require(ggplot2)
require(dplyr)
require(reshape2)
require(gridExtra)
source(paste0(git.dir, "ProcessingScripts/FindMissingDates.R"))

setwd(paste0(git.dir, "data/meteo/"))

# read in CSV
df.1hr <- read.csv("Raw/TFS-EDC_1-hour_data.csv")

# make date column
df.1hr$date <- ymd(df.1hr$date)

# merge datasets from different elevations, using 5 m as master elevation
df.1hr$Tair.C <- df.1hr$air_temp_5m
df.1hr$Tair.C[is.na(df.1hr$Tair.C)] <- df.1hr$air_temp_3m[is.na(df.1hr$Tair.C)]
df.1hr$Tair.C[is.na(df.1hr$Tair.C)] <- df.1hr$air_temp_1m[is.na(df.1hr$Tair.C)]

df.1hr$RH <- df.1hr$rh_5m
df.1hr$RH[is.na(df.1hr$RH)] <- df.1hr$rh_3m[is.na(df.1hr$RH)]
df.1hr$RH[is.na(df.1hr$RH)] <- df.1hr$rh_1m[is.na(df.1hr$RH)]

df.1hr$wind <- df.1hr$wind_sp_5m
df.1hr$wind[is.na(df.1hr$wind)] <- df.1hr$wind_sp_1m[is.na(df.1hr$wind)]

# summarize relevant data
df.1hr.d <- summarize(group_by(df.1hr, date),
                      Tair.C.mean = mean(Tair.C, na.rm=T),
                      Tair.C.min = min(Tair.C, na.rm=T),
                      Tair.C.max = max(Tair.C, na.rm=T),
                      RH = mean(RH),
                      P.kPa = mean(barometer_mbar, na.rm=T)/10,  # convert to kPa
                      wind.m_s = mean(wind, na.rm=T),
                      precip.mm = sum(rain, na.rm=T),
                      rad.W_m2 = mean(pyranometer*1000)) # conver to W/m2

# save as CSV
write.csv(df.1hr.d, "TFS-EDC_1988-2016_Daily.csv", row.names=F)
