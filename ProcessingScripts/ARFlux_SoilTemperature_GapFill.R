## ARFlux_SoilTemperature_GapFill.R
#' This script is intended to gap-fill soil temperature data received via email from Adrian Rocha.

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

## control scripts
max.gapfill.days <- 7  # maximum number of days to fill in via linear interpolation

## load data emailed from Adrian
df.email <- read.csv(paste0(git.dir, "data/ARFlux/Raw/ZipperSoilTempData08-16.csv"), stringsAsFactors=F, skip=1)

# fix column names
colnames(df.email) <- c("date", "Year.DOY", "Severe.5cm", "Moderate.5cm", "Unburned.5cm", 
                        "Severe.10cm", "Moderate.10cm", "Unburned.10cm", 
                        "Severe.20cm", "Moderate.20cm", "Unburned.20cm",
                        "Severe.30cm", "Moderate.30cm", "Unburned.30cm",
                        "Severe.40cm", "Moderate.40cm", "Unburned.40cm",
                        "Depth", "Severe.Air", "Moderate.Air", "Unburned.Air")
df.email$Depth <- NULL

# make date column
df.email$date <- mdy(df.email$date)

# get rid of 2017 data
df.email <- subset(df.email, year(date) != 2017)

# check for missing/duplicated dates
FindMissingDates(df.email$date)
df.email$date[duplicated(df.email$date)]

# gap-fill columns using linear interpolation
for (j in 3:20){
  df.email[,j] <- as.numeric(na.approx(as.numeric(df.email[,j]), maxgap=max.gapfill.days, na.rm=F))
}

# melt to long-form
df.email.melt <- melt(df.email, id=c("date", "Year.DOY"), value.name="Tsoil.C", variable.name="fire.depth", stringsAsFactors=F)

# convert to numeric
df.email.melt$Tsoil.C <- as.numeric(df.email.melt$Tsoil.C)

# get fire and depth info
df.email.melt <- cbind(df.email.melt, as.data.frame(str_split(df.email.melt$fire.depth, pattern="[.]", n=2, simplify=T), stringsAsFactors=F))
colnames(df.email.melt) <- c("date", "Year.DOY", "fire.depth", "Tsoil.C", "fire.sev", "depth.cm")
df.email.melt$depth.cm <- as.numeric(gsub("cm", "", df.email.melt$depth.cm))

# get rid of useless columns
df.email.melt$Year.DOY <- NULL
df.email.melt$fire.depth <- NULL

## plot everything
ggplot(df.email.melt, aes(x=date, y=Tsoil.C)) +
  geom_line() + facet_grid(depth.cm~fire.sev)

## now, gap-fill 5 cm data
# subset to data only
df <- subset(df.email.melt, depth.cm==5)
df$DOY <- yday(df$date)
df$datasource <- "observed"

# get list of dates with data for all sites
dates.sev <- subset(df, is.finite(Tsoil.C) & fire.sev=="Severe")$date
dates.mod <- subset(df, is.finite(Tsoil.C) & fire.sev=="Moderate")$date
dates.unb <- subset(df, is.finite(Tsoil.C) & fire.sev=="Unburned")$date
dates.all <- dates.sev
dates.all <- dates.all[dates.all %in% dates.mod]
dates.all <- dates.all[dates.all %in% dates.unb]

# summarize by DOY and fire
df.dates <- subset(df, date %in% dates.all)
df.d <- dplyr::summarize(group_by(df.dates, DOY, fire.sev),
                         Tsoil.mean = mean(Tsoil.C, na.rm=T))

# merge DOY averages with observed data
df <- left_join(df, df.d, by=c("DOY", "fire.sev"))

# gap-fill data
df$datasource[is.na(df$Tsoil.C)] <- "gapfill"
df$Tsoil.C[is.na(df$Tsoil.C)] <- df$Tsoil.mean[is.na(df$Tsoil.C)]
df$group <- "all"

# plot
ggplot(df, aes(x=date, y=Tsoil.C)) + geom_line(aes(color=datasource, group=group)) + facet_grid(.~fire.sev)

# save output data
write.csv(subset(df, select=c("date", "Tsoil.C", "fire.sev", "datasource")), paste0(git.dir, "data/ARFlux/ARFlux_SoilTemp_2008-2016.csv"), row.names=F)
