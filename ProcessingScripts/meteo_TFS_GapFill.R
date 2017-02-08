## meteo_TFS_GapFill.R
# This gap-fills Toolik Field Station meteorological data with the following rules:
#  1. If data from any of the flux towers are available, they are used in the order: Unburned, Moderate, Severe
#       -This never ends up being used
#  2. If there are <= 3 consecutive days missing, gap-fill via linear interpolation
#       -This never ends up being used
#  3. If there are > 3 consecutive days missing, use the long-term mean for that DOY.

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

require(lubridate)
require(ggplot2)
require(dplyr)
require(reshape2)
require(gridExtra)
source(paste0(git.dir, "ProcessingScripts/FindMissingDates.R"))

setwd(paste0(git.dir, "data"))

# filename of baseline met data
fname.met <- "meteo/Raw/Daily_1988_currentTFSDaily.csv"

# read in baseline met data
df <- read.csv(fname.met, stringsAsFactors=F)

# read in ARFlux data for gap-filling
df.ARFlux <- read.csv(paste0(git.dir, "data/ARFlux/ARFlux_2008-2012_Daily.csv"), stringsAsFactors=F)

# set up column classes
df.ARFlux$Date <- ymd(df.ARFlux$Date)

# data starts in June 1988 and precip data only exists until end of 2013; trim to start January 1 1989 and end December 31 2013
df <- subset(df, Year>=1989 & Year<=2013)

# convert date column
df$Date <- dmy(df$Date)

# collect columns of interest
df.met <- data.frame(Date = df$Date,
                     Year = df$Year,
                     DOY = yday(df$Date),
                     Station = df$Station,
                     Tair.C.mean = as.numeric(df$Daily_AirTemp_Mean_C),
                     Tair.C.min = as.numeric(df$Daily_AirTemp_AbsMin_C),
                     Tair.C.max = as.numeric(df$Daily_AirTemp_AbsMax_C),
                     precip.mm = as.numeric(df$Daily_Precip_Total_mm),
                     wind.m_s = as.numeric(df$Daily_windsp_mean_msec),
                     rad.W_m2 = as.numeric(df$Daily_globalrad_total_jcm2)*100*100/86400,  # convert J/cm2 to W/m2
                     Comments = df$Comments)

# get indices of missing variables
i.missing.Tair.mean <- which(is.na(df.met$Tair.C.mean))
i.missing.Tair.min <- which(is.na(df.met$Tair.C.min))
i.missing.Tair.max <- which(is.na(df.met$Tair.C.max))
i.missing.precip.mm <- which(is.na(df.met$precip.mm))
i.missing.wind.m_s <- which(is.na(df.met$wind.m_s))
i.missing.rad.W_m2 <- which(is.na(df.met$rad.W_m2))

# Gap-Fill Mean Air Temperature ------------------------------------------------

## Rule 1: Gapfilling from fire dataset

# scroll through missing dates
for (i.missing in i.missing.Tair.mean){
  # get missing date
  date.missing <- df.met$Date[i.missing]
  
  # check to see if it exists at any of the ARFlux sites
  i.ARFlux.unburned <- which(df.ARFlux$Date==date.missing & df.ARFlux$fire=="Unburned" & is.finite(df.ARFlux$Tair.C.mean))
  i.ARFlux.moderate <- which(df.ARFlux$Date==date.missing & df.ARFlux$fire=="Moderate" & is.finite(df.ARFlux$Tair.C.mean))
  i.ARFlux.severe <- which(df.ARFlux$Date==date.missing & df.ARFlux$fire=="Severe" & is.finite(df.ARFlux$Tair.C.mean))
  
  # if so, fill it in
  if (length(i.ARFlux.unburned)>0){
    df.met$Tair.C.mean[i.missing] <- df.ARFlux$Tair.C.mean[i.ARFlux.unburned]
    
    # update comments
    df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; Tair.C.mean gap-filled based on ARFlux unburned")
    
  } else if (length(i.ARFlux.moderate)>0){
    df.met$Tair.C.mean[i.missing] <- df.ARFlux$Tair.C.mean[i.ARFlux.moderate]
    
    # update comments
    df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; Tair.C.mean gap-filled based on ARFlux moderate")
    
  } else if (length(i.ARFlux.severe)>0){
    df.met$Tair.C.mean[i.missing] <- df.ARFlux$Tair.C.mean[i.ARFlux.severe]
    
    # update comments
    df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; Tair.C.mean gap-filled based on ARFlux severe")
    
  }
}

## Rule 2: Linear interpolation
# (not necessary, no gaps are <= 3 days)

## Rule 3: Use long-term mean

# get indices of temperature that is still missing
i.missing.Tair.mean <- which(is.na(df.met$Tair.C.mean))

# summarize to mean for each DOY
df.DOY <- summarize(group_by(df.met, DOY),
                    Tair.C.mean = mean(Tair.C.mean, na.rm=T))

# scroll through missing dates
for (i.missing in i.missing.Tair.mean){
  # get missing date
  date.missing <- df.met$Date[i.missing]
  
  # find DOY
  DOY.missing <- yday(date.missing)
  
  # fill in df.met
  df.met$Tair.C.mean[i.missing] <- df.DOY$Tair.C.mean[DOY.missing]
  
  # update comments
  df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; Tair.C.mean gap-filled based on long-term mean for this DOY")
  
}

# check for missing data
i.missing.Tair.mean <- which(is.na(df.met$Tair.C.mean))

# Gap-Fill Min Air Temperature ------------------------------------------------

## Rule 1: Gapfilling from fire dataset

# scroll through missing dates
for (i.missing in i.missing.Tair.min){
  # get missing date
  date.missing <- df.met$Date[i.missing]
  
  # check to see if it exists at any of the ARFlux sites
  i.ARFlux.unburned <- which(df.ARFlux$Date==date.missing & df.ARFlux$fire=="Unburned" & is.finite(df.ARFlux$Tair.C.min))
  i.ARFlux.moderate <- which(df.ARFlux$Date==date.missing & df.ARFlux$fire=="Moderate" & is.finite(df.ARFlux$Tair.C.min))
  i.ARFlux.severe <- which(df.ARFlux$Date==date.missing & df.ARFlux$fire=="Severe" & is.finite(df.ARFlux$Tair.C.min))
  
  # if so, fill it in
  if (length(i.ARFlux.unburned)>0){
    df.met$Tair.C.min[i.missing] <- df.ARFlux$Tair.C.min[i.ARFlux.unburned]
    
    # update comments
    df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; Tair.C.min gap-filled based on ARFlux unburned")
    
  } else if (length(i.ARFlux.moderate)>0){
    df.met$Tair.C.min[i.missing] <- df.ARFlux$Tair.C.min[i.ARFlux.moderate]
    
    # update comments
    df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; Tair.C.min gap-filled based on ARFlux moderate")
    
  } else if (length(i.ARFlux.severe)>0){
    df.met$Tair.C.min[i.missing] <- df.ARFlux$Tair.C.min[i.ARFlux.severe]
    
    # update comments
    df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; Tair.C.min gap-filled based on ARFlux severe")
    
  }
}

## Rule 2: Linear interpolation
# (not necessary, no gaps are <= 3 days)

## Rule 3: Use long-term mean

# get indices of temperature that is still missing
i.missing.Tair.min <- which(is.na(df.met$Tair.C.min))

# summarize to mean for each DOY
df.DOY <- summarize(group_by(df.met, DOY),
                    Tair.C.min = mean(Tair.C.min, na.rm=T))

# scroll through missing dates
for (i.missing in i.missing.Tair.min){
  # get missing date
  date.missing <- df.met$Date[i.missing]
  
  # find DOY
  DOY.missing <- yday(date.missing)
  
  # fill in df.met
  df.met$Tair.C.min[i.missing] <- df.DOY$Tair.C.min[DOY.missing]
  
  # update comments
  df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; Tair.C.min gap-filled based on long-term mean for this DOY")
  
}

# check for missing data
i.missing.Tair.min <- which(is.na(df.met$Tair.C.min))

# Gap-Fill Max Air Temperature ------------------------------------------------

## Rule 1: Gapfilling from fire dataset

# scroll through missing dates
for (i.missing in i.missing.Tair.max){
  # get missing date
  date.missing <- df.met$Date[i.missing]
  
  # check to see if it exists at any of the ARFlux sites
  i.ARFlux.unburned <- which(df.ARFlux$Date==date.missing & df.ARFlux$fire=="Unburned" & is.finite(df.ARFlux$Tair.C.max))
  i.ARFlux.moderate <- which(df.ARFlux$Date==date.missing & df.ARFlux$fire=="Moderate" & is.finite(df.ARFlux$Tair.C.max))
  i.ARFlux.severe <- which(df.ARFlux$Date==date.missing & df.ARFlux$fire=="Severe" & is.finite(df.ARFlux$Tair.C.max))
  
  # if so, fill it in
  if (length(i.ARFlux.unburned)>0){
    df.met$Tair.C.max[i.missing] <- df.ARFlux$Tair.C.max[i.ARFlux.unburned]
    
    # update comments
    df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; Tair.C.max gap-filled based on ARFlux unburned")
    
  } else if (length(i.ARFlux.moderate)>0){
    df.met$Tair.C.max[i.missing] <- df.ARFlux$Tair.C.max[i.ARFlux.moderate]
    
    # update comments
    df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; Tair.C.max gap-filled based on ARFlux moderate")
    
  } else if (length(i.ARFlux.severe)>0){
    df.met$Tair.C.max[i.missing] <- df.ARFlux$Tair.C.max[i.ARFlux.severe]
    
    # update comments
    df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; Tair.C.max gap-filled based on ARFlux severe")
    
  }
}

## Rule 2: Linear interpolation
# (not necessary, no gaps are <= 3 days)

## Rule 3: Use long-term mean

# get indices of temperature that is still missing
i.missing.Tair.max <- which(is.na(df.met$Tair.C.max))

# summarize to mean for each DOY
df.DOY <- summarize(group_by(df.met, DOY),
                    Tair.C.max = mean(Tair.C.max, na.rm=T))

# scroll through missing dates
for (i.missing in i.missing.Tair.max){
  # get missing date
  date.missing <- df.met$Date[i.missing]
  
  # find DOY
  DOY.missing <- yday(date.missing)
  
  # fill in df.met
  df.met$Tair.C.max[i.missing] <- df.DOY$Tair.C.max[DOY.missing]
  
  # update comments
  df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; Air temp gap-filled based on long-term mean for this DOY")
  
}

# check for missing data
i.missing.Tair.max <- which(is.na(df.met$Tair.C.max))

# Gap-Fill Precip ------------------------------------------------

## Rule 1: Gapfilling from fire dataset
# no precip data at flux towers

## Rule 3: Use long-term mean

# get indices of temperature that is still missing
i.missing.precip.mm <- which(is.na(df.met$precip.mm))

# summarize to mean for each DOY
df.DOY <- summarize(group_by(df.met, DOY),
                    precip.mm = mean(precip.mm, na.rm=T))

# scroll through missing dates
for (i.missing in i.missing.precip.mm){
  # get missing date
  date.missing <- df.met$Date[i.missing]
  
  # find DOY
  DOY.missing <- yday(date.missing)
  
  # fill in df.met
  df.met$precip.mm[i.missing] <- df.DOY$precip.mm[DOY.missing]
  
  # update comments
  df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; precip.mm gap-filled based on long-term mean for this DOY")
  
}

# check for missing data
i.missing.precip.mm <- which(is.na(df.met$precip.mm))

# Gap-Fill Wind Speed ------------------------------------------------

## Rule 1: Gapfilling from fire dataset

# scroll through missing dates
for (i.missing in i.missing.wind.m_s){
  # get missing date
  date.missing <- df.met$Date[i.missing]
  
  # check to see if it exists at any of the ARFlux sites
  i.ARFlux.unburned <- which(df.ARFlux$Date==date.missing & df.ARFlux$fire=="Unburned" & is.finite(df.ARFlux$windSpeed.m2.s))
  i.ARFlux.moderate <- which(df.ARFlux$Date==date.missing & df.ARFlux$fire=="Moderate" & is.finite(df.ARFlux$windSpeed.m2.s))
  i.ARFlux.severe <- which(df.ARFlux$Date==date.missing & df.ARFlux$fire=="Severe" & is.finite(df.ARFlux$windSpeed.m2.s))
  
  # if so, fill it in
  if (length(i.ARFlux.unburned)>0){
    df.met$wind.m_s[i.missing] <- df.ARFlux$windSpeed.m2.s[i.ARFlux.unburned]
    
    # update comments
    df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; wind.m_s gap-filled based on ARFlux unburned")
    
  } else if (length(i.ARFlux.moderate)>0){
    df.met$wind.m_s[i.missing] <- df.ARFlux$windSpeed.m2.s[i.ARFlux.moderate]
    
    # update comments
    df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; wind.m_s gap-filled based on ARFlux moderate")
    
  } else if (length(i.ARFlux.severe)>0){
    df.met$wind.m_s[i.missing] <- df.ARFlux$windSpeed.m2.s[i.ARFlux.severe]
    
    # update comments
    df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; wind.m_s gap-filled based on ARFlux severe")
    
  }
}

## Rule 3: Use long-term mean

# get indices of temperature that is still missing
i.missing.wind.m_s <- which(is.na(df.met$wind.m_s))

# summarize to mean for each DOY
df.DOY <- summarize(group_by(df.met, DOY),
                    wind.m_s = mean(wind.m_s, na.rm=T))

# scroll through missing dates
for (i.missing in i.missing.wind.m_s){
  # get missing date
  date.missing <- df.met$Date[i.missing]
  
  # find DOY
  DOY.missing <- yday(date.missing)
  
  # fill in df.met
  df.met$wind.m_s[i.missing] <- df.DOY$wind.m_s[DOY.missing]
  
  # update comments
  df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; wind.m_s gap-filled based on long-term mean for this DOY")
  
}

# check for missing data
i.missing.wind.m_s <- which(is.na(df.met$wind.m_s))

# Gap-Fill Global Radiation ------------------------------------------------

## Rule 1: Gapfilling from fire dataset

# scroll through missing dates
for (i.missing in i.missing.rad.W_m2){
  # get missing date
  date.missing <- df.met$Date[i.missing]
  
  # check to see if it exists at any of the ARFlux sites
  i.ARFlux.unburned <- which(df.ARFlux$Date==date.missing & df.ARFlux$fire=="Unburned" & is.finite(df.ARFlux$SWin.W.m2))
  i.ARFlux.moderate <- which(df.ARFlux$Date==date.missing & df.ARFlux$fire=="Moderate" & is.finite(df.ARFlux$SWin.W.m2))
  i.ARFlux.severe <- which(df.ARFlux$Date==date.missing & df.ARFlux$fire=="Severe" & is.finite(df.ARFlux$SWin.W.m2))
  
  # if so, fill it in
  if (length(i.ARFlux.unburned)>0){
    df.met$rad.W_m2[i.missing] <- df.ARFlux$SWin.W.m2[i.ARFlux.unburned]
    
    # update comments
    df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; rad.W_m2 gap-filled based on ARFlux unburned")
    
  } else if (length(i.ARFlux.moderate)>0){
    df.met$rad.W_m2[i.missing] <- df.ARFlux$SWin.W.m2[i.ARFlux.moderate]
    
    # update comments
    df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; rad.W_m2 gap-filled based on ARFlux moderate")
    
  } else if (length(i.ARFlux.severe)>0){
    df.met$rad.W_m2[i.missing] <- df.ARFlux$SWin.W.m2[i.ARFlux.severe]
    
    # update comments
    df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; rad.W_m2 gap-filled based on ARFlux severe")
    
  }
}

## Rule 3: Use long-term mean

# get indices of temperature that is still missing
i.missing.rad.W_m2 <- which(is.na(df.met$rad.W_m2))

# summarize to mean for each DOY
df.DOY <- summarize(group_by(df.met, DOY),
                    rad.W_m2 = mean(rad.W_m2, na.rm=T))

# scroll through missing dates
for (i.missing in i.missing.rad.W_m2){
  # get missing date
  date.missing <- df.met$Date[i.missing]
  
  # find DOY
  DOY.missing <- yday(date.missing)
  
  # fill in df.met
  df.met$rad.W_m2[i.missing] <- df.DOY$rad.W_m2[DOY.missing]
  
  # update comments
  df.met$Comments[i.missing] <- paste0(df.met$Comments[i.missing], "; rad.W_m2 gap-filled based on long-term mean for this DOY")
  
}

# check for missing data
i.missing.rad.W_m2 <- which(is.na(df.met$rad.W_m2))

# Save output data --------------------------------------------------------

# save gap-filled
write.csv(df.met, paste0(git.dir, "data/meteo/Daily_1989-2013_TFS_GapFill.csv"), row.names=F)

## make some plots
df.met.long <- melt(df.met, id.vars=c("Date", "Year", "DOY", "Station", "Comments"))
p.met <- 
  ggplot(df.met.long, aes(x=DOY, y=value)) +
  geom_point(aes(color=Year), alpha=0.1) +
  facet_wrap(~variable, scales="free_y") +
  scale_x_continuous(name="Day of Year", expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")

ggsave(paste0(git.dir, "data/meteo/Daily_1989-2013_TFS_GapFill.png"), 
       p.met, width=8, height=6, units="in")
