## geotop_output_ThetaTempProfileTimeseries.R
#' This script makes plots of soil liquid and ice content and 
#' temperature through time.

rm(list=ls())

# git directory for relative paths
#git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"
git.dir <- "C:/Users/Sam/WorkGits/ARF1D/"

require(lubridate)
require(ggplot2)
require(reshape2)
require(dplyr)
require(gridExtra)

#version
version <- "20170323-1hr-HomoSoilVG"

# burn severity site
fire <- "Unburned"  # options are: Unburned, Moderate, Severe

# path to soil input file
path.soil <- paste0(git.dir, "geotop/soil/SoilARF0001.txt")

# modeled data path
path.mod.temp <- paste0(git.dir, "geotop/output-tabs/soiltemp0001.txt")
path.mod.VWC <- paste0(git.dir, "geotop/output-tabs/thetaliq0001.txt")
path.mod.ice <- paste0(git.dir, "geotop/output-tabs/thetaice0001.txt")
path.mod.point <- paste0(git.dir, "geotop/output-tabs/point0001.txt")
path.mod.discharge <- paste0(git.dir, "geotop/output-tabs/discharge.txt")

# read in data
df.soil <- read.csv(path.soil)
df.mod.temp <- read.csv(path.mod.temp, stringsAsFactors=F)
df.mod.VWC <- read.csv(path.mod.VWC, stringsAsFactors=F)
df.mod.ice <- read.csv(path.mod.ice, stringsAsFactors=F)
df.mod.point <- read.csv(path.mod.point, stringsAsFactors=F)
df.mod.discharge <- read.csv(path.mod.discharge, stringsAsFactors=F)

# make Date column
df.mod.temp$Date <- dmy_hm(df.mod.temp$Date12.DDMMYYYYhhmm.)
df.mod.VWC$Date <- dmy_hm(df.mod.VWC$Date12.DDMMYYYYhhmm.)
df.mod.ice$Date <- dmy_hm(df.mod.ice$Date12.DDMMYYYYhhmm.)
df.mod.point$Date <- dmy_hm(df.mod.point$Date12.DDMMYYYYhhmm.)
df.mod.discharge$Date <- dmy_hm(df.mod.discharge$DATE.day.month.year.hour.min.)

# for profile files, keep column for Date and anything beginning with "X" (these are the depths)
all.cols <- colnames(df.mod.temp)
cols.keep <- c("Date", all.cols[startsWith(all.cols, "X")])

# figure out depth for each column
cols.depth <- as.numeric(sub("X", "", cols.keep))

# add depth to soil column
df.soil$depth.mm <- cols.depth[is.finite(cols.depth)]

# trim
df.mod.temp <- df.mod.temp[,cols.keep]
df.mod.VWC <- df.mod.VWC[,cols.keep]
df.mod.ice <- df.mod.ice[,cols.keep]

# melt
df.mod.temp.melt <- melt(df.mod.temp, id=c("Date"))
df.mod.VWC.melt <- melt(df.mod.VWC, id=c("Date"))
df.mod.ice.melt <- melt(df.mod.ice, id=c("Date"))

# calculate depth
df.mod.temp.melt$depth.mm <- as.numeric(sub("X", "", df.mod.temp.melt$variable))
df.mod.VWC.melt$depth.mm <- as.numeric(sub("X", "", df.mod.VWC.melt$variable))
df.mod.ice.melt$depth.mm <- as.numeric(sub("X", "", df.mod.ice.melt$variable))

# add soil data
df.mod.temp.melt <- merge(df.mod.temp.melt, df.soil, by=c("depth.mm"), all.x=T)
df.mod.VWC.melt <- merge(df.mod.VWC.melt, df.soil, by=c("depth.mm"), all.x=T)
df.mod.ice.melt <- merge(df.mod.ice.melt, df.soil, by=c("depth.mm"), all.x=T)

# combine into a single data frame
df.all <- data.frame(Date = df.mod.temp.melt$Date,
                     depth.mm = df.mod.temp.melt$depth.mm,
                     Dz = df.mod.temp.melt$Dz,
                     porosity = df.mod.temp.melt$vwc_s,
                     Tsoil.C = df.mod.temp.melt$value,
                     VWC.liq = df.mod.VWC.melt$value,
                     VWC.ice = df.mod.ice.melt$value)

# calculate saturation
df.all$sat.liq <- df.all$VWC.liq/df.all$porosity
df.all$sat.ice <- df.all$VWC.ice/df.all$porosity
df.all$sat.air <- (1-df.all$sat.ice-df.all$sat.liq)
df.all$sat.air[df.all$sat.air<0] <- 0

# max annual thaw depth (used to test for stability)
df.mod.point$year <- year(df.mod.point$Date)
df.mod.point.yr <- summarize(group_by(df.mod.point, year),
                             thaw.depth.max = max(lowest_thawed_soil_depth.mm.))
p.thaw.yr <-
  ggplot(df.mod.point.yr, aes(x=year, y=thaw.depth.max)) +
  geom_line() +
  scale_y_reverse(name="Max Annual Thaw Depth [mm]") +
  scale_x_continuous(name="Year") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.temp <-
  ggplot() +
  geom_tile(data=df.all, aes(x=Date, y=depth.mm, fill=Tsoil.C, height=Dz)) +
  geom_line(data=df.mod.point, aes(x=Date, y=lowest_thawed_soil_depth.mm.), color="black") +
  scale_y_reverse(name="Depth [mm]", expand=c(0,0)) +
  scale_x_datetime(expand=c(0,0)) +
  scale_fill_gradient2(name="Temperature [C]", low="blue", high="red") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")

p.VWC.liq <-
  ggplot() +
  geom_tile(data=df.all, aes(x=Date, y=depth.mm, fill=VWC.liq, height=Dz)) +
  geom_line(data=df.mod.point, aes(x=Date, y=highest_water_table_depth.mm.), color="black") +
  geom_line(data=df.mod.point, aes(x=Date, y=lowest_water_table_depth.mm.), color="black") +
  scale_y_reverse(name="Depth [mm]", expand=c(0,0)) +
  scale_x_datetime(expand=c(0,0)) +
  scale_fill_gradient(name="VWC [m3/m3]", high="blue", low="red") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")

p.VWC.ice <-
  ggplot() +
  geom_tile(data=df.all, aes(x=Date, y=depth.mm, fill=VWC.ice, height=Dz)) +
  geom_line(data=df.mod.point, aes(x=Date, y=highest_water_table_depth.mm.), color="black") +
  geom_line(data=df.mod.point, aes(x=Date, y=lowest_water_table_depth.mm.), color="black") +
  scale_y_reverse(name="Depth [mm]", expand=c(0,0)) +
  scale_x_datetime(expand=c(0,0)) +
  scale_fill_gradient(name="Ice [m3/m3]", high="blue", low="red") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")

p.sat.liq <-
  ggplot() +
  geom_tile(data=df.all, aes(x=Date, y=depth.mm, fill=sat.liq, height=Dz)) +
  geom_line(data=df.mod.point, aes(x=Date, y=highest_water_table_depth.mm.), color="black") +
  geom_line(data=df.mod.point, aes(x=Date, y=lowest_water_table_depth.mm.), color="black") +
  scale_y_reverse(name="Depth [mm]", expand=c(0,0)) +
  scale_x_datetime(expand=c(0,0)) +
  scale_fill_gradient(name="Pore Liquid Fraction [-]", high="blue", low="red", limits=c(0,1)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")

p.sat.ice <-
  ggplot() +
  geom_tile(data=df.all, aes(x=Date, y=depth.mm, fill=sat.ice, height=Dz)) +
  geom_line(data=df.mod.point, aes(x=Date, y=highest_water_table_depth.mm.), color="black") +
  geom_line(data=df.mod.point, aes(x=Date, y=lowest_water_table_depth.mm.), color="black") +
  scale_y_reverse(name="Depth [mm]", expand=c(0,0)) +
  scale_x_datetime(expand=c(0,0)) +
  scale_fill_gradient(name="Pore Ice Fraction [-]", high="blue", low="red", limits=c(0,1)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")

p.sat.air <-
  ggplot() +
  geom_tile(data=df.all, aes(x=Date, y=depth.mm, fill=sat.air, height=Dz)) +
  geom_line(data=df.mod.point, aes(x=Date, y=highest_water_table_depth.mm.), color="black") +
  geom_line(data=df.mod.point, aes(x=Date, y=lowest_water_table_depth.mm.), color="black") +
  scale_y_reverse(name="Depth [mm]", expand=c(0,0)) +
  scale_x_datetime(expand=c(0,0)) +
  scale_fill_gradient(name="Pore Air Fraction [-]", high="red", low="blue", limits=c(0,1)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")