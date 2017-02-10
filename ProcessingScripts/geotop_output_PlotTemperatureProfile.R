## geotop_output_PlotTemperatureProfile.R
#' This script is intended to read in model output from GEOtop with
#' temperature profiles and plot them.

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
path.mod <- paste0(git.dir, "geotop/output-tabs/soiltemp0001.txt")

# save plot path
path.plot <- paste0(git.dir, "geotop/output-plots/TemperatureProfiles_", version, ".png")

## read in data
df.mod <- read.csv(path.mod, stringsAsFactors=F)

## rearrange model data into long-form data frame
# keep column for TimeFromStart and anything beginning with "X" (these are the depths)
all.cols <- colnames(df.mod)
colnames(df.mod)[startsWith(all.cols, "TimeFromStart")] <- "TimeFromStart"
all.cols <- colnames(df.mod)
cols.keep <- all.cols[startsWith(all.cols, "TimeFromStart") | startsWith(all.cols, "X")]

# melt into long-form data frame
df.mod.long <- melt(df.mod[,cols.keep], id.vars="TimeFromStart", variable.name="depth.mm", value.name="Tsoil")

# get rid of leading X in column depth.mm
df.mod.long$depth.mm <- as.numeric(sub("X", "", df.mod.long$depth.mm))

# get rid of time=0
df.mod.long <- subset(df.mod.long, TimeFromStart>0)

## plots for temperature profiles at different timesteps
p.Tsoil.depth.year <- 
  ggplot(df.mod.long, aes(y=Tsoil, x=depth.mm/1000, color=factor(round(TimeFromStart)))) +
  geom_line(alpha=0.5) +
  scale_y_continuous(name="Soil Temperature [C]") +
  scale_x_reverse(name="Depth [m]", expand=c(0,0)) +
  scale_color_manual(name="Year", values=colorRampPalette(c("#116611", "#162955"))(length(unique(df.mod.long$TimeFromStart)))) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid=element_blank())

ggsave(path.plot, p.Tsoil.depth.year,
       width=6, height=6, units="in")
