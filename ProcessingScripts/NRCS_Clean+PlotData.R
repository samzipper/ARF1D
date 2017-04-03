## NRCS_Clean+PlotData.R
#' This script cleans and plots data from the NRCS soil climate station
#' at Toolik Lake.

rm(list=ls())

# git directory for relative paths
#git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"
git.dir <- "C:/Users/Sam/WorkGits/ARF1D/"

require(ggplot2)
require(gridExtra)
require(lubridate)
require(stringr)

# path to raw data
path.in <- paste0(git.dir, "data/NRCS/Toolik_AllYears.csv")

# path to save cleaned output and plot
path.out <- paste0(git.dir, "data/NRCS/Toolik_AllYears_Clean.csv")
path.plot <- paste0(git.dir, "data/NRCS/Toolik_AllYears_Clean_Plot.png")

# read in raw data, skipping top 3 lines
df.in <- read.csv(path.in, skip=3, stringsAsFactors=F)

# trim to only relevant columns
df.in <- df.in[,c("year", "DOY", "date", "temp.air.C", "precip.in", 
                  "temp.soil.120mm.campbell.C", "temp.soil.280mm.campbell.C", "temp.soil.650mm.campbell.C", 
                  "temp.soil.6mm.C", "temp.soil.87mm.C", 
                  "temp.soil.160mm.C", "temp.soil.236mm.C", "temp.soil.312mm.C", 
                  "temp.soil.387mm.C", "temp.soil.463mm.C", "temp.soil.616mm.C", 
                  "temp.soil.768mm.C", "temp.soil.978mm.C", "temp.soil.120mm.C", "VWC.120mm", 
                  "temp.soil.90mm.C", "VWC.90mm", "temp.soil.280mm.C.r1", "VWC.280mm.r1", 
                  "temp.soil.280mm.C.r2", "VWC.280mm.r2", "temp.soil.390mm.C", "VWC.390mm", 
                  "temp.soil.680mm.C", "VWC.680mm")]
df.in$date <- mdy(df.in$date)

# 2004 VWC data needs to be divided by 100
df.in$VWC.90mm[df.in$year==2004] <- df.in$VWC.90mm[df.in$year==2004]/100
df.in$VWC.120mm[df.in$year==2004] <- df.in$VWC.120mm[df.in$year==2004]/100
df.in$VWC.280mm.r1[df.in$year==2004] <- df.in$VWC.280mm.r1[df.in$year==2004]/100
df.in$VWC.280mm.r2[df.in$year==2004] <- df.in$VWC.280mm.r2[df.in$year==2004]/100
df.in$VWC.390mm[df.in$year==2004] <- df.in$VWC.390mm[df.in$year==2004]/100
df.in$VWC.680mm[df.in$year==2004] <- df.in$VWC.680mm[df.in$year==2004]/100

# vitel temperature probe minimum is -12.5 C (see comparison between campbell & vitel for 120 mm) - set values less than 12 C to NaN
df.in$temp.soil.6mm.C[df.in$temp.soil.6mm.C< -12.0] <- NaN
df.in$temp.soil.87mm.C[df.in$temp.soil.87mm.C< -12.0] <- NaN
df.in$temp.soil.160mm.C[df.in$temp.soil.160mm.C< -12.0] <- NaN
df.in$temp.soil.236mm.C[df.in$temp.soil.236mm.C< -12.0] <- NaN
df.in$temp.soil.312mm.C[df.in$temp.soil.312mm.C< -12.0] <- NaN
df.in$temp.soil.387mm.C[df.in$temp.soil.387mm.C< -12.0] <- NaN
df.in$temp.soil.463mm.C[df.in$temp.soil.463mm.C< -12.0] <- NaN
df.in$temp.soil.616mm.C[df.in$temp.soil.616mm.C< -12.0] <- NaN
df.in$temp.soil.768mm.C[df.in$temp.soil.768mm.C< -12.0] <- NaN
df.in$temp.soil.978mm.C[df.in$temp.soil.978mm.C< -12.0] <- NaN
df.in$temp.soil.120mm.C[df.in$temp.soil.120mm.C< -12.0] <- NaN
df.in$temp.soil.90mm.C[df.in$temp.soil.90mm.C< -12.0] <- NaN
df.in$temp.soil.280mm.C.r1[df.in$temp.soil.280mm.C.r1< -12.0] <- NaN
df.in$temp.soil.280mm.C.r2[df.in$temp.soil.280mm.C.r2< -12.0] <- NaN
df.in$temp.soil.390mm.C[df.in$temp.soil.390mm.C< -12.0] <- NaN
df.in$temp.soil.680mm.C[df.in$temp.soil.680mm.C< -12.0] <- NaN

# plot comparing things at same depth
p.temp <-
  ggplot(df.in, aes(x=date)) +
  geom_line(aes(y=temp.soil.90mm.C), color="red")

p.air.120mm.temp <-
  ggplot(df.in, aes(x=date)) +
  geom_line(aes(y=temp.soil.120mm.campbell.C), color="red") +
  geom_line(aes(y=temp.air.70mm.C), color="blue")

p.120mm.temp <-
  ggplot(df.in, aes(x=date)) +
  geom_line(aes(y=temp.soil.120mm.campbell.C), color="red") +
  geom_line(aes(y=temp.soil.120mm.C), color="blue")

p.280mm.temp <-
  ggplot(df.in, aes(x=date)) +
  geom_line(aes(y=temp.soil.280mm.C.r1), color="red") +
  geom_line(aes(y=temp.soil.280mm.C.r2), color="blue")

p.280mm.VWC <-
  ggplot(df.in, aes(x=date)) +
  geom_line(aes(y=VWC.280mm.r1), color="red") +
  geom_line(aes(y=VWC.280mm.r2), color="blue")

## melt
df.melt <- melt(df.in, id=c("year", "DOY", "date", "temp.air.C", "precip.in"))

# split up variable name
var.name <- strsplit(as.character(df.melt$variable), "[.]")

df.melt$depth.mm <- as.numeric(str_replace(unlist(var.name), "mm", ""))[is.finite(as.numeric(str_replace(unlist(var.name), "mm", "")))]
df.melt$variable <- sapply(var.name, "[[", 1)

# for duplicates at same date, get min/mean/max
df.melt.d <- summarize(group_by(df.melt, variable, depth.mm, date),
                       value.min = min(value, na.rm=T),
                       value.mean = mean(value, na.rm=T),
                       value.max = max(value, na.rm=T))

# put in order
df.melt.d <- df.melt.d[order(df.melt.d$variable, df.melt.d$depth.mm, df.melt.d$date),]

# save output file
write.csv(df.melt.d, path.out, row.names=F)

# make plot
p.temp.facet <-
  ggplot(subset(df.melt.d, variable=="temp"), aes(x=date, y=value.mean)) +
  geom_line() +
  facet_wrap(~depth.mm)

p.VWC.facet <-
  ggplot(subset(df.melt.d, variable=="VWC"), aes(x=date, y=value.mean)) +
  geom_line() +
  facet_wrap(~depth.mm)
