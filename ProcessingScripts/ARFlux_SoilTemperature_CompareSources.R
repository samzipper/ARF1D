## ARFlux_SoilTemperature_CompareSources.R
#' This script is intended to compare soil temperature data downloaded from LTER website
#' with soil temperature data received via email from Adrian Rocha.

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
source(paste0(git.dir, "ProcessingScripts/FindMissingDates.R"))

## load data downloaded from LTER website
for (fire in c("Moderate", "Severe", "Unburned")){
  # load data
  df.f <- read.csv(paste0(git.dir, "data/ARFlux/Raw/ARFlux_2008-2012_", fire, ".csv"), stringsAsFactors=F)
  
  # retain only soil temp
  df.f <- subset(df.f, select=c("Year", "DOY", "Soil.Temperature"))
  
  # make fire column
  df.f$fire.sev <- fire
  
  if (exists("df.web")){
    df.web <- rbind(df.web, df.f)
  } else {
    df.web <- df.f
  }
}

# make date column
df.web$datetime <- ymd(paste0(df.web$Year, "-01-01")) + seconds(round(df.web$DOY*86400))
df.web$date <- as.Date(df.web$datetime)

# get rid of bad data
df.web$Soil.Temperature[df.web$Soil.Temperature < -25 | df.web$Soil.Temperature > 25] <- NaN

# summarize to daily
df.web.d <- dplyr::summarize(group_by(subset(df.web, select=c("date", "fire.sev", "Soil.Temperature")), fire.sev, date),
                             Tsoil.C = mean(Soil.Temperature, na.rm=T))

# make data source column
df.web.d$datasource <- "web"

# make depth column (just put as 5 cm for comparison)
df.web.d$depth.cm <- 5

# plot
ggplot(df.web.d, aes(x=date, y=Tsoil.C, color=fire.sev)) + geom_line()

## load data emailed from Adrian
df.email <- read.csv(paste0(git.dir, "data/ARFlux/Raw/ZipperSoilTempData08-16.csv"), stringsAsFactors=F, skip=1)

# fix column names
colnames(df.email) <- c("date", "Year.DOY", "Severe.5cm", "Moderate.5cm", "Unburned.5cm", 
                        "Severe.10cm", "Moderate.10cm", "Unburned.10cm", 
                        "Severe.20cm", "Moderate.20cm", "Unburned.20cm",
                        "Severe.30cm", "Moderate.30cm", "Unburned.30cm",
                        "Severe.40cm", "Moderate.40cm", "Unburned.40cm",
                        "Depth", "Severe.Air", "Moderate.Air", "Unburned.Air")

# make date column
df.email$date <- mdy(df.email$date)

# check for missing/duplicated dates
FindMissingDates(df.email$date)
df.email$date[duplicated(df.email$date)]

# melt to long-form
df.email.melt <- melt(df.email, id=c("date", "Year.DOY"), value.name="Tsoil.C", variable.name="fire.depth", stringsAsFactors=F)

# conver to numeric
df.email.melt$Tsoil.C <- as.numeric(df.email.melt$Tsoil.C)

# get fire and depth info
df.email.melt <- cbind(df.email.melt, as.data.frame(str_split(df.email.melt$fire.depth, pattern="[.]", n=2, simplify=T), stringsAsFactors=F))
colnames(df.email.melt) <- c("date", "Year.DOY", "fire.depth", "Tsoil.C", "fire.sev", "depth.cm")
df.email.melt$depth.cm <- as.numeric(gsub("cm", "", df.email.melt$depth.cm))

# make data source column
df.email.melt$datasource <- "email"

# get rid of useless columns
df.email.melt$Year.DOY <- NULL
df.email.melt$fire.depth <- NULL

## combine data frames
df <- rbind.data.frame(df.web.d, df.email.melt)

## plot
ggplot(subset(df, depth.cm==5 & year(date)<=2012), aes(x=date, y=Tsoil.C, color=datasource)) +
  geom_line() + facet_grid(fire.sev~.)

ggplot(subset(df, depth.cm==5), aes(x=date, y=Tsoil.C, color=datasource)) +
  geom_line() + facet_grid(fire.sev~.)
