## geotop_SoilPropertyCalibration_SelectBest+Plot_NRCS.R
#' This is intended to compare modeled and measured data:
#'   -Soil temperature
#'   -Soil moisture
#'   
#'   These data are from the NRCS soil climate site.

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/SoilPropertyCalibrationBasic/ARF1D/"
#git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

# which geotop version
geo.dir <- "geotop_NRCS/"
#geo.dir <- "geotop/"

require(lubridate)
require(ggplot2)
require(reshape2)
require(dplyr)
require(gridExtra)
source(paste0(git.dir, "ProcessingScripts/FitMetrics.R"))
which.closest <- function(x, vec){
  # function returns the index (i) of the element in the vector (vec) which is closest to the input value (x).
  # if x is equidistant from multiple elements in x, all are returned.
  # a tolerance parameter is built in to avoid rounding errors.
  # 
  eps <- .Machine$double.eps^0.5  # tolerance parameter
  
  i <- which(abs(vec-x)-min(abs(vec-x), na.rm=T) < eps)
  
  return(i)
}

# flags for what to do
process <- F   # process model output data (if false: load already-processed model output data)
plot <- T      # make plots (if false: no plots)

# path to save plot
path.plot <- paste0(git.dir, geo.dir, "SoilPropertyCalibration/")

## calibrate with pre-2004, validate with post-2004
cal <- seq(1998,2004,1)
val <- seq(2005,2011,1)

yr.min <- min(c(cal, val))
yr.max <- max(c(cal, val))

## load observation data and preprocess
# paths
path.obs <- paste0(git.dir, "data/NRCS/Toolik_AllYears_Clean.csv")

# read in data
df.obs <- read.csv(path.obs, stringsAsFactors=F)

colnames(df.obs)[colnames(df.obs)=="date"] <- "Date"

# add date column
df.obs$Date <- ymd(df.obs$Date)

# subset to only calibration/validation period
df.obs <- subset(df.obs, year(Date)>=yr.min & year(Date)<=yr.max)

# bad data from 2004: temp 90 cm, temp 390 cm, temp 680 cm
df.obs <- df.obs[!(year(df.obs$Date)==2004 & df.obs$variable=="temp" & df.obs$depth.mm %in% c(90,390,680)),]

# depths to plot for temp and VWC
depths.temp <- unique(df.obs$depth.mm[is.finite(df.obs$value.mean) & df.obs$variable=="temp"])
depths.VWC <- unique(df.obs$depth.mm[is.finite(df.obs$value.mean) & df.obs$variable=="VWC"])

# process data
if (process){
  
  # numbers to process
  numbers.all <- seq(1,100)

  ## read in input data with soil hydraulic properties
  df.in <- read.csv(paste0(git.dir, geo.dir, "SoilPropertyCalibration/SoilPropertyCalibration_Input.csv"))
  df.in$number <- sprintf("%04d", df.in$number)
  
  ## Process model output data
  for (num in numbers.all){
    # parameter set number
    number <- sprintf("%04d", num)
    
    # modeled data path
    path.mod.temp <- paste0(git.dir, geo.dir, "SoilPropertyCalibration/output_", number, "_soiltemp0001.txt")
    path.mod.VWC <- paste0(git.dir, geo.dir, "SoilPropertyCalibration/output_", number, "_thetaliq0001.txt")
    
    # read in data
    df.mod.temp <- read.csv(path.mod.temp, stringsAsFactors=F)
    df.mod.VWC <- read.csv(path.mod.VWC, stringsAsFactors=F)
    
    # make Date column
    df.mod.temp$Date <- as.Date(dmy_hm(df.mod.temp$Date12.DDMMYYYYhhmm.))
    df.mod.VWC$Date <- as.Date(dmy_hm(df.mod.VWC$Date12.DDMMYYYYhhmm.))
    
    # check if model finished or crashed
    if (ymd(paste0(year(max(df.mod.temp$Date)), "-", month(max(df.mod.temp$Date)), "-", day(max(df.mod.temp$Date))))==as.Date(ymd("2011-12-31"))){
      mod.finish <- T
    } else {
      mod.finish <- F
    }
    
    # subset
    df.mod.temp <- subset(df.mod.temp, year(Date)>=yr.min & year(Date)<=yr.max)
    df.mod.VWC <- subset(df.mod.VWC, year(Date)>=yr.min & year(Date)<=yr.max)
    
    if (mod.finish){
      
      # keep column for Date and anything beginning with "X" (these are the depths)
      all.cols <- colnames(df.mod.temp)
      cols.keep <- c("Date", all.cols[startsWith(all.cols, "X")])
      df.mod.temp <- df.mod.temp[,cols.keep]
      df.mod.VWC <- df.mod.VWC[,cols.keep]
      
      # reformat as long-form
      df.mod.temp.melt <- melt(df.mod.temp, id="Date", value.name="value.mod")
      df.mod.temp.melt$depth.mm <- as.numeric(sub("X", "", df.mod.temp.melt$variable))
      
      df.mod.VWC.melt <- melt(df.mod.VWC, id="Date", value.name="value.mod")
      df.mod.VWC.melt$depth.mm <- as.numeric(sub("X", "", df.mod.VWC.melt$variable))
      
      # add column names to match df.obs
      df.mod.temp.melt$variable <- "temp"
      df.mod.VWC.melt$variable <- "VWC"
      
      # combine
      df.mod.melt <- rbind(df.mod.temp.melt, df.mod.VWC.melt)
      
      # merge with obs data frame
      df.all <- merge(df.obs, df.mod.melt, by=c("Date", "variable", "depth.mm"), all.y=T)
      
      # get rid of entries without observations results
      df.all <- subset(df.all, is.finite(value.mean))
      
      # make calibration/validation column
      df.all$period <- "cal"
      df.all$period[year(df.all$Date) %in% val] <- "val"
      
      # calculate fit metrics
      df.fit <- summarize(group_by(subset(df.all, is.finite(value.mean)), variable, depth.mm, period),
                          RMSE = RMSE(value.mod, value.mean),
                          NRMSE = NRMSE(value.mod, value.mean),
                          R2 = R2(value.mod, value.mean),
                          NSE = NashSutcliffe(value.mod, value.mean))
      
      # add some extra data
      df.fit$number <- number
      
    }
    
    # add to overall output data frame
    if (exists("df.fit.all")){
      df.fit.all <- rbind(df.fit.all, df.fit)
    } else {
      df.fit.all <- df.fit
    }
    
    # status update
    print(paste0(number, " complete"))
  }
  
  # summarize by number, variable, and period
  df.fit.all.out <- 
    summarize(group_by(df.fit.all, number, period, variable),
              RMSE.mean = mean(RMSE),
              NRMSE.mean = mean(NRMSE),
              R2.mean = mean(R2),
              NSE.mean = mean(NSE))
  
  # merge with input
  df.out <- merge(df.in, df.fit.all.out, by="number", all=T)
  
  # save output files
  write.csv(df.fit, paste0(git.dir, geo.dir, "SoilPropertyCalibration/SoilPropertyCalibration_FitMetrics_All.csv"), row.names=F)
  write.csv(df.out, paste0(git.dir, geo.dir, "SoilPropertyCalibration/SoilPropertyCalibration_FitMetrics_Summary.csv"), row.names=F)
  
}  # end of process loop

if (plot){
  # read in output data
  df.out <- read.csv(paste0(git.dir, geo.dir, "SoilPropertyCalibration/SoilPropertyCalibration_FitMetrics_Summary.csv"), stringsAsFactors=F)
  
  # determine best as that which has the lowest combined NSE or calibration period
  df.out.total <- summarize(group_by(df.out, number, period),
                            NSE = sum(NSE.mean),
                            NRMSE = mean(NRMSE.mean))
  
  ## several different ways of choosing best overall...
  # max combined NSE
  #n.best.overall <- df.out.total$number[df.out.total$period=="cal"][which.max(df.out.total$NSE[df.out.total$period=="cal"])]
  
  # min mean NRMSE
  n.best.overall <- df.out.total$number[df.out.total$period=="cal"][which.min(df.out.total$NRMSE[df.out.total$period=="cal"])]
  
  # lowest VWC RMSE
  #n.best.overall <- df.out$number[df.out$period=="cal" & df.out$variable=="VWC"][which.min(df.out$RMSE[df.out$period=="cal" & df.out$variable=="VWC"])]
  
  ## Process model output data
  # parameter set number
  number <- sprintf("%04d", n.best.overall)
  
  # modeled data path
  path.mod.temp <- paste0(git.dir, geo.dir, "SoilPropertyCalibration/output_", number, "_soiltemp0001.txt")
  path.mod.VWC <- paste0(git.dir, geo.dir, "SoilPropertyCalibration/output_", number, "_thetaliq0001.txt")
  
  # read in data
  df.mod.temp <- read.csv(path.mod.temp, stringsAsFactors=F)
  df.mod.VWC <- read.csv(path.mod.VWC, stringsAsFactors=F)
  
  # make Date column
  df.mod.temp$Date <- as.Date(dmy_hm(df.mod.temp$Date12.DDMMYYYYhhmm.))
  df.mod.VWC$Date <- as.Date(dmy_hm(df.mod.VWC$Date12.DDMMYYYYhhmm.))
  
  # subset
  df.mod.temp <- subset(df.mod.temp, year(Date)>=yr.min & year(Date)<=yr.max)
  df.mod.VWC <- subset(df.mod.VWC, year(Date)>=yr.min & year(Date)<=yr.max)
  
  # keep column for Date and anything beginning with "X" (these are the depths)
  all.cols <- colnames(df.mod.temp)
  cols.keep <- c("Date", all.cols[startsWith(all.cols, "X")])
  df.mod.temp <- df.mod.temp[,cols.keep]
  df.mod.VWC <- df.mod.VWC[,cols.keep]
  
  # reformat as long-form
  df.mod.temp.melt <- melt(df.mod.temp, id="Date", value.name="value.mod")
  df.mod.temp.melt$depth.mm <- as.numeric(sub("X", "", df.mod.temp.melt$variable))
  
  df.mod.VWC.melt <- melt(df.mod.VWC, id="Date", value.name="value.mod")
  df.mod.VWC.melt$depth.mm <- as.numeric(sub("X", "", df.mod.VWC.melt$variable))
  
  # add column names to match df.obs
  df.mod.temp.melt$variable <- "temp"
  df.mod.VWC.melt$variable <- "VWC"
  
  # combine
  df.mod.melt <- rbind(df.mod.temp.melt, df.mod.VWC.melt)
  
  # merge with obs data frame
  df.all <- merge(df.obs, df.mod.melt, by=c("Date", "variable", "depth.mm"), all.y=T)
  
  # make calibration/validation column
  df.all$period <- "cal"
  df.all$period[year(df.all$Date) %in% val] <- "val"
  
  # calculate fit metrics
  df.fit <- summarize(group_by(subset(df.all, is.finite(value.mean)), variable, depth.mm, period),
                      RMSE = RMSE(value.mod, value.mean),
                      NRMSE = NRMSE(value.mod, value.mean),
                      R2 = R2(value.mod, value.mean),
                      NSE = NashSutcliffe(value.mod, value.mean))
  df.fit.melt <- melt(df.fit, id=c("variable", "depth.mm", "period"), value.name="fit", variable.name="metric")
  
  # plot
  p.temp.all <-
    ggplot(subset(df.all, variable=="temp" & depth.mm %in% depths.temp), aes(x=Date)) +
    geom_hline(yintercept=0, color="gray65") +
    geom_line(aes(y=value.mod), color="black") +
    geom_line(aes(y=value.mean), color="red") +
    facet_wrap(~depth.mm) +
    scale_x_date(expand=c(0,0)) +
    theme_bw() +
    theme(panel.grid=element_blank())
  ggsave(paste0(path.plot, "SoilPropertyCalibration_BestFit_p.temp.all_", sprintf("%04d", n.best.overall), ".png"), p.temp.all, width=12, height=6, units="in")
  
  p.VWC.all <-
    ggplot(subset(df.all, variable=="VWC" & depth.mm %in% depths.VWC), aes(x=Date)) +
    geom_hline(yintercept=0, color="gray65") +
    geom_line(aes(y=value.mod), color="black") +
    geom_line(aes(y=value.mean), color="blue") +
    facet_wrap(~depth.mm) +
    scale_x_date(expand=c(0,0)) +
    theme_bw() +
    theme(panel.grid=element_blank())
  ggsave(paste0(path.plot, "SoilPropertyCalibration_BestFit_p.VWC.all_", sprintf("%04d", n.best.overall), ".png"), p.VWC.all, width=8, height=6, units="in")
  
  # plot of only common depths
  p.temp.VWC <-
    ggplot(subset(df.all, depth.mm %in% depths.VWC), aes(x=Date)) +
    geom_hline(yintercept=0, color="gray65") +
    geom_line(aes(y=value.mod), color="black") +
    geom_line(aes(y=value.mean), color="forestgreen") +
    facet_grid(variable~depth.mm, scales="free") +
    scale_x_date(expand=c(0,0)) +
    theme_bw() +
    theme(panel.grid=element_blank())
  ggsave(paste0(path.plot, "SoilPropertyCalibration_BestFit_p.temp.VWC_", sprintf("%04d", n.best.overall), ".png"), p.temp.VWC, width=10, height=6, units="in")
  
  # plot fit
  p.fit <-
    ggplot(df.fit.melt, aes(x=factor(depth.mm), y=fit, fill=period)) +
    geom_hline(yintercept=0, color="gray65") +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(metric~variable, scales="free") +
    scale_x_discrete(name="Depth [mm]") +
    scale_y_continuous(name="Fit") +
    scale_fill_discrete(name="Period") +
    theme_bw() +
    theme(panel.grid=element_blank(),
          legend.position="bottom")
  ggsave(paste0(path.plot, "SoilPropertyCalibration_BestFit_p.fit_", sprintf("%04d", n.best.overall), ".png"), p.fit, width=8, height=6, units="in")
  
}  # end of plot loop