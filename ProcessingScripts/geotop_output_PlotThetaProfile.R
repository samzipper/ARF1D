## geotop_output_PlotThetaProfile.R
#' This script is intended to read in model output from GEOtop with
#' soil moisture profiles and plot them.

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

require(ggplot2)
require(reshape2)
require(dplyr)
require(gridExtra)

# version name
version <- "20170210-1mo-NoFlow"

## paths
# modeled data path
path.liq <- paste0(git.dir, "geotop/output-tabs/thetaliq0001.txt")
path.ice <- paste0(git.dir, "geotop/output-tabs/thetaice0001.txt")

# save plot path
path.plot <- paste0(git.dir, "geotop/output-plots/ThetaProfiles_", version, ".png")

## read in data
df.liq <- read.csv(path.liq, stringsAsFactors=F)
df.ice <- read.csv(path.ice, stringsAsFactors=F)

## rearrange model data into long-form data frame
# keep column for TimeFromStart and anything beginning with "X" (these are the depths)
all.cols <- colnames(df.liq)
colnames(df.liq)[startsWith(all.cols, "TimeFromStart")] <- "TimeFromStart"
colnames(df.ice)[startsWith(all.cols, "TimeFromStart")] <- "TimeFromStart"
all.cols <- colnames(df.liq)
cols.keep <- all.cols[startsWith(all.cols, "TimeFromStart") | startsWith(all.cols, "X")]

# melt into long-form data frame
df.liq.long <- melt(df.liq[,cols.keep], id.vars="TimeFromStart", variable.name="depth.mm", value.name="liq")
df.ice.long <- melt(df.ice[,cols.keep], id.vars="TimeFromStart", variable.name="depth.mm", value.name="ice")

# get rid of leading X in column depth.mm
df.liq.long$depth.mm <- as.numeric(sub("X", "", df.liq.long$depth.mm))
df.ice.long$depth.mm <- as.numeric(sub("X", "", df.ice.long$depth.mm))

# merge liquid and ice
df.long <- merge(df.liq.long, df.ice.long)

# calculate liquid+ice saturation and liquid fraction
df.long$sat <- df.long$liq + df.long$ice
df.long$liq.fraction <- df.long$liq/df.long$sat

# get rid of time=0
df.long <- subset(df.long, TimeFromStart>0)

## plots for temperature profiles at different timesteps
p.sat.depth.year <- 
  ggplot(df.long, aes(y=sat, x=depth.mm/1000, color=factor(round(TimeFromStart)))) +
  geom_line(alpha=0.5) +
  scale_y_continuous(name="Liquid + Ice Content [m3/m3]") +
  scale_x_reverse(name="Depth [m]", expand=c(0,0)) +
  scale_color_manual(name="Year", values=colorRampPalette(c("#116611", "#162955"))(length(unique(df.liq.long$TimeFromStart)))) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid=element_blank())

ggsave(path.plot, p.sat.depth.year,
       width=6, height=6, units="in")
