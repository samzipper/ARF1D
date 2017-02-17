## meteo_TFS-EDC_AggregateToDaily.R
# This aggregates 1 hour and 3 hour Toolik Field Data Environmental Data Center
# meteorological data to daily means/sums.
#
# Raw data from here: http://toolik.alaska.edu/edc/abiotic_monitoring/data_query.php
# Conditions of use: http://toolik.alaska.edu/edc/about/conditions_of_use.php

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

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

# summarize relevant data
df.1hr.d <- summarize(group_by(df.1hr, date),
                      Tair.C.mean = mean(air_temp_3m, na.rm=T),
                      Tair.C.min = min(air_temp_3m, na.rm=T),
                      Tair.C.max = max(air_temp_3m, na.rm=T),
                      RH = mean(rh_3m),
                      P.kPa = mean(barometer_mbar, na.rm=T)/10,  # convert to kPa
                      precip.mm = sum(rain, na.rm=T))

# save as CSV
write.csv(df.1hr.d, "TFS-EDC_1988-2015_Daily.csv", row.names=F)
