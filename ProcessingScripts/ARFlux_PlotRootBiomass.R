## ARFlux_PlotRootBiomass.R
#' This script is intended to read in a root biomass dataset 
#' (from http://arc-lter.ecosystems.mbl.edu/2011arfrootbiomasscn-byquad )
#' and make some plots.

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

require(dplyr)
require(reshape2)
require(ggplot2)
require(stringr)

# path to root biomass CSV file
path.in <- paste0(git.dir, "data/ARFlux/ARFlux_RootBiomass.csv")

# read in root biomass CSV file
df.in <- read.csv(path.in, stringsAsFactors=F)

# renames columns of interest
colnames(df.in)[colnames(df.in)=="Fire.disturbance"] <- "fire"
colnames(df.in)[colnames(df.in)=="Organic.layer.depth..cm."] <- "org.depth.cm"
colnames(df.in)[colnames(df.in)=="Organic.Total.root.biomass..g.m2."] <- "org.biomass.g_m2"
colnames(df.in)[colnames(df.in)=="Mineral.layer.depth..cm."] <- "min.depth.cm"
colnames(df.in)[colnames(df.in)=="Mineral.Total.root.biomass..g.m2."] <- "min.biomass.g_m2"

# save only columns of interest
df.in <- df.in[,c("fire","org.depth.cm","org.biomass.g_m2","min.depth.cm","min.biomass.g_m2")]

# reshape to plot
df.in.melt <- melt(df.in, id.vars=c("fire"), stringsAsFactors=F)

# convert columns
df.in.melt$variable <- as.character(levels(df.in.melt$variable)[as.numeric(df.in.melt$variable)])
df.in.melt$variable[str_detect(df.in.melt$variable, "biomass")] <- "biomass.g_m2"
df.in.melt$variable[str_detect(df.in.melt$variable, "depth")] <- "depth.cm"

# recast
df.in.cast <- dcast(df.in.melt, fire ~ variable, value.var="value")

# plot
p.root.org <- 
  ggplot(df.in, aes(x=org.biomass.g_m2, y=org.depth.cm)) +
  geom_point() +
  facet_wrap(~fire, ncol=3, scales="free") +
  stat_smooth(method="loess") +
  scale_y_reverse() +
  theme_bw()

p.root.min <- 
  ggplot(df.in, aes(x=min.biomass.g_m2, y=min.depth.cm)) +
  geom_point() +
  facet_wrap(~fire, ncol=3, scales="free") +
  stat_smooth(method="loess") +
  scale_y_reverse() +
  theme_bw()
  