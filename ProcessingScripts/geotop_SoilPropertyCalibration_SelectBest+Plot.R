## geotop_SoilPropertyCalibration_SelectBest+Plot.R
#' This is intended to compare modeled and measured data:
#'   -Soil temperature
#'   -Soil moisture
#'   -Thaw depth
#'   
#'   These data are from the ARFlux sites only.

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

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
process <- T   # process model output data (if false: load already-processed model output data)
plot <- T      # make plots (if false: no plots)

## calibrate with odd years, validate with even
# this is done because the only winter with observational data
# is 2008-2009, and I want to have winter in both calibration
# and validation
cal <- c(2008, 2010, 2012)
val <- c(2009, 2011, 2013)

yr.min <- min(c(cal, val))
yr.max <- max(c(cal, val))

# process data
if (process){
  
  # numbers to process
  numbers.all <- seq(1,26)
  
  # read in input data with soil hydraulic properties
  df.in <- read.csv(paste0(git.dir, "geotop/SoilPropertyCalibration/SoilPropertyCalibration_Input.csv"))
  
  # loop through sites
  for (fire in c("Unburned", "Moderate", "Severe")){
    
    ## read in observed data and preprocess
    # paths
    path.obs.ARFlux <- paste0(git.dir, "data/ARFlux/ARFlux-Merged_2008-2012_Daily.csv")
    path.thaw <- paste0(git.dir, "data/ARFlux/ARFlux_ThawDepths_2008-2014.csv")
    
    # read in data
    df.obs.ARFlux <- read.csv(path.obs.ARFlux, stringsAsFactors=F)
    df.thaw <- read.csv(path.thaw, stringsAsFactors=F)
    
    # add date column
    df.obs.ARFlux$Date <- ymd(df.obs.ARFlux$Date)
    df.thaw$Date <- dmy(df.thaw$Date)
    
    # subset to only calibration/validation period
    df.obs.ARFlux <- subset(df.obs.ARFlux, year(Date)>=yr.min & year(Date)<=yr.max)
    df.thaw <- subset(df.thaw, year(Date)>=yr.min & year(Date)<=yr.max)
    
    ## pre-processing of observed data
    # trim observations to only this fire
    df.obs.ARFlux <- df.obs.ARFlux[df.obs.ARFlux$fire==fire,]
    df.thaw <- data.frame(Date = df.thaw$Date,
                          ThawDepth.mm = 10*as.numeric(df.thaw[,which(colnames(df.thaw)==paste0("Thaw.depth.", fire))]))
    
    # ARFlux depths are listed on repository, e.g. http://arc-lter.ecosystems.mbl.edu/2012arfluxunburned
    depths.ARFlux <- c(20, 60)  # depth is average of 2 cm and 6 cm sensor
    depths.ARFlux.VWC <- c(25)  # VWC data is at 2.5 cm
    
    # calculate mean depth
    depth.ARFlux.min <- min(depths.ARFlux)
    depth.ARFlux.mean <- mean(depths.ARFlux)
    depth.ARFlux.max <- max(depths.ARFlux)
    
    # trim ARFlux data frame
    df.obs.ARFlux <- df.obs.ARFlux[,c("Date", "Tsoil.C", "Tsoil.C.min", "Tsoil.C.max", "VWC", "VWC.min", "VWC.max")]
    
    ## Process model output data
    for (num in numbers.all){
      # parameter set number
      number <- sprintf("%04d", num)
      
      # modeled data path
      path.mod.temp <- paste0(git.dir, "geotop/SoilPropertyCalibration/output_", number, "_", fire, "_soiltemp0001.txt")
      path.mod.VWC <- paste0(git.dir, "geotop/SoilPropertyCalibration/output_", number, "_", fire, "_thetaliq0001.txt")
      
      # read in data
      df.mod.temp <- read.csv(path.mod.temp, stringsAsFactors=F)
      df.mod.VWC <- read.csv(path.mod.VWC, stringsAsFactors=F)
      
      # make Date column
      df.mod.temp$Date <- dmy_hm(df.mod.temp$Date12.DDMMYYYYhhmm.)
      df.mod.VWC$Date <- dmy_hm(df.mod.VWC$Date12.DDMMYYYYhhmm.)
      
      # check if model finished or crashed
      if (ymd(paste0(year(max(df.mod.temp$Date)), "-", month(max(df.mod.temp$Date)), "-", day(max(df.mod.temp$Date))))==as.Date(ymd("2013-12-31"))){
        mod.finish <- T
      } else {
        mod.finish <- F
      }
      
      # subset
      df.mod.temp <- subset(df.mod.temp, year(Date)>=yr.min & year(Date)<=yr.max)
      df.mod.VWC <- subset(df.mod.VWC, year(Date)>=yr.min & year(Date)<=yr.max)
      
      # check if any data is left - if not, model crashed!
      if (mod.finish){
        
        # format Date column
        df.mod.temp$Date <- format(df.mod.temp$Date, "%m/%d/%Y")
        df.mod.VWC$Date <- format(df.mod.VWC$Date, "%m/%d/%Y")
        
        # keep column for Date and anything beginning with "X" (these are the depths)
        all.cols <- colnames(df.mod.temp)
        cols.keep <- c("Date", all.cols[startsWith(all.cols, "X")])
        df.mod.temp <- df.mod.temp[,cols.keep]
        df.mod.VWC <- df.mod.VWC[,cols.keep]
        
        # figure out depth for each column
        cols.depth <- as.numeric(sub("X", "", cols.keep))
        
        # determine which columns are within the depth range of the observations
        cols.compare.ARFlux <- unlist(lapply(depths.ARFlux, FUN=which.closest, vec=cols.depth))
        cols.compare.ARFlux.VWC <- unlist(lapply(depths.ARFlux.VWC, FUN=which.closest, vec=cols.depth))
        
        # summarize model for each day
        Dates.all <- unique(df.mod.temp$Date)
        df.mod.day <- data.frame(Date = Dates.all, 
                                 Temp.ARFlux.mean = NaN,
                                 Temp.ARFlux.std = NaN,
                                 Temp.ARFlux.min = NaN,
                                 Temp.ARFlux.max = NaN,
                                 VWC.ARFlux.mean = NaN,
                                 VWC.ARFlux.std = NaN,
                                 VWC.ARFlux.min = NaN,
                                 VWC.ARFlux.max = NaN,
                                 ThawDepth.mm = NaN)
        for (d in Dates.all){
          # figure out index for this date
          i.d <- which(df.mod.day$Date==d)
          
          # summarize VWC and temperature
          df.mod.day$Temp.ARFlux.mean[i.d] <- mean(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
          df.mod.day$Temp.ARFlux.std[i.d] <- sd(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
          df.mod.day$Temp.ARFlux.min[i.d] <- min(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
          df.mod.day$Temp.ARFlux.max[i.d] <- max(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
          
          df.mod.day$VWC.ARFlux.mean[i.d] <- mean(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)
          df.mod.day$VWC.ARFlux.std[i.d] <- sd(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)
          df.mod.day$VWC.ARFlux.min[i.d] <- min(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)
          df.mod.day$VWC.ARFlux.max[i.d] <- max(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)
          
          ## figure out thaw depth: shallowest depth at which temperature > 0
          # extract all temperatures for that date
          df.mod.temp.profile <- data.frame(depth = cols.depth[is.finite(cols.depth)], 
                                            Temp = unlist(df.mod.temp[df.mod.temp$Date==d, is.finite(cols.depth)]))
          
          # find the shallowest temperature that is frozen
          df.mod.day$ThawDepth.mm[i.d] <- min(df.mod.temp.profile$depth[which(df.mod.temp.profile$Temp<0)])
        }
        
        ## make plots
        # convert dates
        df.mod.day$Date <- mdy(df.mod.day$Date)
        df.obs.ARFlux$Date <- ymd(df.obs.ARFlux$Date)
        
        # for each thaw date, summarize to mean and range
        df.thaw.day <- summarize(group_by(df.thaw, Date),
                                 ThawDepth.mm.mean = mean(ThawDepth.mm, na.rm=T),
                                 ThawDepth.mm.std = sd(ThawDepth.mm, na.rm=T),
                                 ThawDepth.mm.min = min(ThawDepth.mm, na.rm=T),
                                 ThawDepth.mm.max = max(ThawDepth.mm, na.rm=T))
        
        
        ## calculate fit metrics - calibration period
        # ARFlux
        df.fit.temp.ARFlux <- df.obs.ARFlux[is.finite(df.obs.ARFlux$Tsoil.C),]
        df.fit.temp.ARFlux$Temp.mod <- df.mod.day$Temp.ARFlux.mean[match(df.fit.temp.ARFlux$Date, df.mod.day$Date)]
        
        df.fit.VWC.ARFlux <- df.obs.ARFlux[is.finite(df.obs.ARFlux$VWC),]
        df.fit.VWC.ARFlux$VWC.mod <- df.mod.day$VWC.ARFlux.mean[match(df.fit.VWC.ARFlux$Date, df.mod.day$Date)]
        
        df.thaw.day$thaw.mod <- df.mod.day$ThawDepth.mm[match(df.thaw.day$Date, df.mod.day$Date)]
        
        # make calibration/validation column
        df.fit.temp.ARFlux$period <- "cal"
        df.fit.temp.ARFlux$period[year(df.fit.temp.ARFlux$Date) %in% val] <- "val"
        
        df.fit.VWC.ARFlux$period <- "cal"
        df.fit.VWC.ARFlux$period[year(df.fit.VWC.ARFlux$Date) %in% val] <- "val"
        
        df.thaw.day$period <- "cal"
        df.thaw.day$period[year(df.thaw.day$Date) %in% val] <- "val"
        
        # calculate fit statistics
        
        # Reload and plot best ----------------------------------------------------
        
        # make a matrix with fit statistics
        df.fit.table.cal <- data.frame(RMSE=numeric(3), NRMSE=numeric(3), NSE=numeric(3), R2=numeric(3),
                                       row.names=c("Thaw Depth","Tower Soil Temp", "Tower VWC"))
        df.fit.table.cal["Tower Soil Temp",] <- c(RMSE(subset(df.fit.temp.ARFlux, period=="cal")$Temp.mod, subset(df.fit.temp.ARFlux, period=="cal")$Tsoil.C),
                                                  NRMSE(subset(df.fit.temp.ARFlux, period=="cal")$Temp.mod, subset(df.fit.temp.ARFlux, period=="cal")$Tsoil.C),
                                                  NashSutcliffe(subset(df.fit.temp.ARFlux, period=="cal")$Temp.mod, subset(df.fit.temp.ARFlux, period=="cal")$Tsoil.C),
                                                  R2(subset(df.fit.temp.ARFlux, period=="cal")$Temp.mod, subset(df.fit.temp.ARFlux, period=="cal")$Tsoil.C))
        df.fit.table.cal["Tower VWC",] <- c(RMSE(subset(df.fit.VWC.ARFlux, period=="cal")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="cal")$VWC),
                                            NRMSE(subset(df.fit.VWC.ARFlux, period=="cal")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="cal")$VWC),
                                            NashSutcliffe(subset(df.fit.VWC.ARFlux, period=="cal")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="cal")$VWC),
                                            R2(subset(df.fit.VWC.ARFlux, period=="cal")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="cal")$VWC))
        df.fit.table.cal["Thaw Depth",] <- c(RMSE(subset(df.thaw.day, period=="cal")$thaw.mod, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean),
                                             NRMSE(subset(df.thaw.day, period=="cal")$thaw.mod, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean),
                                             NashSutcliffe(subset(df.thaw.day, period=="cal")$thaw.mod, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean),
                                             R2(subset(df.thaw.day, period=="cal")$thaw.mod, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean))
        
        
        df.fit.table.val <- data.frame(RMSE=numeric(3), NRMSE=numeric(3), NSE=numeric(3), R2=numeric(3),
                                       row.names=c("Thaw Depth","Tower Soil Temp", "Tower VWC"))
        df.fit.table.val["Tower Soil Temp",] <- c(RMSE(subset(df.fit.temp.ARFlux, period=="val")$Temp.mod, subset(df.fit.temp.ARFlux, period=="val")$Tsoil.C),
                                                  NRMSE(subset(df.fit.temp.ARFlux, period=="val")$Temp.mod, subset(df.fit.temp.ARFlux, period=="val")$Tsoil.C),
                                                  NashSutcliffe(subset(df.fit.temp.ARFlux, period=="val")$Temp.mod, subset(df.fit.temp.ARFlux, period=="val")$Tsoil.C),
                                                  R2(subset(df.fit.temp.ARFlux, period=="val")$Temp.mod, subset(df.fit.temp.ARFlux, period=="val")$Tsoil.C))
        df.fit.table.val["Tower VWC",] <- c(RMSE(subset(df.fit.VWC.ARFlux, period=="val")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="val")$VWC),
                                            NRMSE(subset(df.fit.VWC.ARFlux, period=="val")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="val")$VWC),
                                            NashSutcliffe(subset(df.fit.VWC.ARFlux, period=="val")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="val")$VWC),
                                            R2(subset(df.fit.VWC.ARFlux, period=="val")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="val")$VWC))
        df.fit.table.val["Thaw Depth",] <- c(RMSE(subset(df.thaw.day, period=="val")$thaw.mod, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean),
                                             NRMSE(subset(df.thaw.day, period=="val")$thaw.mod, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean),
                                             NashSutcliffe(subset(df.thaw.day, period=="val")$thaw.mod, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean),
                                             R2(subset(df.thaw.day, period=="val")$thaw.mod, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean))
        
        # melt and merge
        df.fit.cal <- df.fit.table.cal
        df.fit.cal$period <- "cal"
        df.fit.cal$var <- row.names(df.fit.cal)
        df.fit.val <- df.fit.table.val
        df.fit.val$period <- "val"
        df.fit.val$var <- row.names(df.fit.val)
        df.fit.out <- rbind(df.fit.cal, df.fit.val)  # combine cal & val
        df.fit.out$fire <- fire
        df.fit.out$number <- number
        
      } else {
        # if model failed, just put some NaNs in there
        df.fit.out <- data.frame(fire = fire, 
                                 number = number,
                                 RMSE = NaN,
                                 NRMSE = NaN,
                                 R2 = NaN,
                                 NSE = NaN,
                                 period = "fail",
                                 var = "fail")
      }
      
      # add to overall output data frame
      if (exists("df.fit.out.all")){
        df.fit.out.all <- rbind(df.fit.out.all, df.fit.out)
      } else {
        df.fit.out.all <- df.fit.out
      }
      
      # status update
      print(paste0(fire, " ", number, " complete"))
    }
  }
  
  # merge with df.in
  df.in$number <- sprintf("%04d", df.in$number)
  df.out <- merge(df.in, df.fit.out.all, all.y=T)
  
  # save output file
  write.csv(df.out, paste0(git.dir, "geotop/SoilPropertyCalibration/SoilPropertyCalibration_Output.csv"), row.names=F)
  
}  # end of process loop

if (plot){
  # read in output data
  df.out <- read.csv(paste0(git.dir, "geotop/SoilPropertyCalibration/SoilPropertyCalibration_Output.csv"), stringsAsFactors=F)
  
  # determine best as that which has the lowest combined NSE
  df.out.fire <- summarize(group_by(df.out, number, fire),
                           NSE.sum = sum(NSE),
                           n.cal = sum(period=="cal"),
                           n.val = sum(period=="val"))
  df.out.fire <- subset(df.out.fire, n.cal+n.val==6)
  n.best.severe <- subset(df.out.fire, fire=="Severe")$number[which.max(subset(df.out.fire, fire=="Severe")$NSE.sum)]
  n.best.moderate <- subset(df.out.fire, fire=="Moderate")$number[which.max(subset(df.out.fire, fire=="Moderate")$NSE.sum)]
  n.best.unburned <- subset(df.out.fire, fire=="Unburned")$number[which.max(subset(df.out.fire, fire=="Unburned")$NSE.sum)]
  
  df.out.number <- summarize(group_by(df.out, number),
                             NSE.sum = sum(NSE),
                             n.cal = sum(period=="cal"),
                             n.val = sum(period=="val"))
  n.best.overall <- df.out.number$number[which.max(df.out.number$NSE.sum)]
  
  ## for best overall, load and plot
  # load observations
  # loop through sites
  for (fire in c("Unburned", "Moderate", "Severe")){
    
    ## read in observed data and preprocess
    # paths
    path.obs.ARFlux <- paste0(git.dir, "data/ARFlux/ARFlux-Merged_2008-2012_Daily.csv")
    path.thaw <- paste0(git.dir, "data/ARFlux/ARFlux_ThawDepths_2008-2014.csv")
    
    # read in data
    df.obs.ARFlux <- read.csv(path.obs.ARFlux, stringsAsFactors=F)
    df.thaw <- read.csv(path.thaw, stringsAsFactors=F)
    
    # add date column
    df.obs.ARFlux$Date <- ymd(df.obs.ARFlux$Date)
    df.thaw$Date <- dmy(df.thaw$Date)
    
    # subset to only calibration/validation period
    df.obs.ARFlux <- subset(df.obs.ARFlux, year(Date)>=yr.min & year(Date)<=yr.max)
    df.thaw <- subset(df.thaw, year(Date)>=yr.min & year(Date)<=yr.max)
    
    ## pre-processing of observed data
    # trim observations to only this fire
    df.obs.ARFlux <- df.obs.ARFlux[df.obs.ARFlux$fire==fire,]
    df.thaw <- data.frame(Date = df.thaw$Date,
                          ThawDepth.mm = 10*as.numeric(df.thaw[,which(colnames(df.thaw)==paste0("Thaw.depth.", fire))]))
    
    # ARFlux depths are listed on repository, e.g. http://arc-lter.ecosystems.mbl.edu/2012arfluxunburned
    depths.ARFlux <- c(20, 60)  # depth is average of 2 cm and 6 cm sensor
    depths.ARFlux.VWC <- c(25)  # VWC data is at 2.5 cm
    
    # calculate mean depth
    depth.ARFlux.min <- min(depths.ARFlux)
    depth.ARFlux.mean <- mean(depths.ARFlux)
    depth.ARFlux.max <- max(depths.ARFlux)
    
    # trim ARFlux data frame
    df.obs.ARFlux <- df.obs.ARFlux[,c("Date", "Tsoil.C", "Tsoil.C.min", "Tsoil.C.max", "VWC", "VWC.min", "VWC.max")]
    
    ## Process model output data
    # parameter set number
    number <- sprintf("%04d", n.best.overall)
    
    # modeled data path
    path.mod.temp <- paste0(git.dir, "geotop/SoilPropertyCalibration/output_", number, "_", fire, "_soiltemp0001.txt")
    path.mod.VWC <- paste0(git.dir, "geotop/SoilPropertyCalibration/output_", number, "_", fire, "_thetaliq0001.txt")
    
    # read in data
    df.mod.temp <- read.csv(path.mod.temp, stringsAsFactors=F)
    df.mod.VWC <- read.csv(path.mod.VWC, stringsAsFactors=F)
    
    # make Date column
    df.mod.temp$Date <- dmy_hm(df.mod.temp$Date12.DDMMYYYYhhmm.)
    df.mod.VWC$Date <- dmy_hm(df.mod.VWC$Date12.DDMMYYYYhhmm.)
    
    # subset
    df.mod.temp <- subset(df.mod.temp, year(Date)>=yr.min & year(Date)<=yr.max)
    df.mod.VWC <- subset(df.mod.VWC, year(Date)>=yr.min & year(Date)<=yr.max)
    
    # format Date column
    df.mod.temp$Date <- format(df.mod.temp$Date, "%m/%d/%Y")
    df.mod.VWC$Date <- format(df.mod.VWC$Date, "%m/%d/%Y")
    
    # keep column for Date and anything beginning with "X" (these are the depths)
    all.cols <- colnames(df.mod.temp)
    cols.keep <- c("Date", all.cols[startsWith(all.cols, "X")])
    df.mod.temp <- df.mod.temp[,cols.keep]
    df.mod.VWC <- df.mod.VWC[,cols.keep]
    
    # figure out depth for each column
    cols.depth <- as.numeric(sub("X", "", cols.keep))
    
    # determine which columns are within the depth range of the observations
    cols.compare.ARFlux <- unlist(lapply(depths.ARFlux, FUN=which.closest, vec=cols.depth))
    cols.compare.ARFlux.VWC <- unlist(lapply(depths.ARFlux.VWC, FUN=which.closest, vec=cols.depth))
    
    # summarize model for each day
    Dates.all <- unique(df.mod.temp$Date)
    df.mod.day <- data.frame(Date = Dates.all, 
                             Temp.ARFlux.mean = NaN,
                             Temp.ARFlux.std = NaN,
                             Temp.ARFlux.min = NaN,
                             Temp.ARFlux.max = NaN,
                             VWC.ARFlux.mean = NaN,
                             VWC.ARFlux.std = NaN,
                             VWC.ARFlux.min = NaN,
                             VWC.ARFlux.max = NaN,
                             ThawDepth.mm = NaN)
    for (d in Dates.all){
      # figure out index for this date
      i.d <- which(df.mod.day$Date==d)
      
      # summarize VWC and temperature
      df.mod.day$Temp.ARFlux.mean[i.d] <- mean(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
      df.mod.day$Temp.ARFlux.std[i.d] <- sd(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
      df.mod.day$Temp.ARFlux.min[i.d] <- min(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
      df.mod.day$Temp.ARFlux.max[i.d] <- max(unlist(df.mod.temp[df.mod.temp$Date==d, cols.compare.ARFlux]), na.rm=T)
      
      df.mod.day$VWC.ARFlux.mean[i.d] <- mean(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)
      df.mod.day$VWC.ARFlux.std[i.d] <- sd(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)
      df.mod.day$VWC.ARFlux.min[i.d] <- min(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)
      df.mod.day$VWC.ARFlux.max[i.d] <- max(unlist(df.mod.VWC[df.mod.VWC$Date==d, cols.compare.ARFlux.VWC]), na.rm=T)
      
      ## figure out thaw depth: shallowest depth at which temperature > 0
      # extract all temperatures for that date
      df.mod.temp.profile <- data.frame(depth = cols.depth[is.finite(cols.depth)], 
                                        Temp = unlist(df.mod.temp[df.mod.temp$Date==d, is.finite(cols.depth)]))
      
      # find the shallowest temperature that is frozen
      df.mod.day$ThawDepth.mm[i.d] <- min(df.mod.temp.profile$depth[which(df.mod.temp.profile$Temp<0)])
    }
    
    ## make plots
    # convert dates
    df.mod.day$Date <- mdy(df.mod.day$Date)
    df.obs.ARFlux$Date <- ymd(df.obs.ARFlux$Date)
    
    # for each thaw date, summarize to mean and range
    df.thaw.day <- summarize(group_by(df.thaw, Date),
                             ThawDepth.mm.mean = mean(ThawDepth.mm, na.rm=T),
                             ThawDepth.mm.std = sd(ThawDepth.mm, na.rm=T),
                             ThawDepth.mm.min = min(ThawDepth.mm, na.rm=T),
                             ThawDepth.mm.max = max(ThawDepth.mm, na.rm=T))
    
    
    ## calculate fit metrics - calibration period
    # ARFlux
    df.fit.temp.ARFlux <- df.obs.ARFlux[is.finite(df.obs.ARFlux$Tsoil.C),]
    df.fit.temp.ARFlux$Temp.mod <- df.mod.day$Temp.ARFlux.mean[match(df.fit.temp.ARFlux$Date, df.mod.day$Date)]
    
    df.fit.VWC.ARFlux <- df.obs.ARFlux[is.finite(df.obs.ARFlux$VWC),]
    df.fit.VWC.ARFlux$VWC.mod <- df.mod.day$VWC.ARFlux.mean[match(df.fit.VWC.ARFlux$Date, df.mod.day$Date)]
    
    df.thaw.day$thaw.mod <- df.mod.day$ThawDepth.mm[match(df.thaw.day$Date, df.mod.day$Date)]
    
    # make calibration/validation column
    df.fit.temp.ARFlux$period <- "cal"
    df.fit.temp.ARFlux$period[year(df.fit.temp.ARFlux$Date) %in% val] <- "val"
    
    df.fit.VWC.ARFlux$period <- "cal"
    df.fit.VWC.ARFlux$period[year(df.fit.VWC.ARFlux$Date) %in% val] <- "val"
    
    df.thaw.day$period <- "cal"
    df.thaw.day$period[year(df.thaw.day$Date) %in% val] <- "val"
    
    # calculate fit statistics
    
    # Reload and plot best ----------------------------------------------------
    
    # make a matrix with fit statistics
    df.fit.table.cal <- data.frame(RMSE=numeric(3), NRMSE=numeric(3), NSE=numeric(3), R2=numeric(3),
                                   row.names=c("Thaw Depth","Tower Soil Temp", "Tower VWC"))
    df.fit.table.cal["Tower Soil Temp",] <- c(RMSE(subset(df.fit.temp.ARFlux, period=="cal")$Temp.mod, subset(df.fit.temp.ARFlux, period=="cal")$Tsoil.C),
                                              NRMSE(subset(df.fit.temp.ARFlux, period=="cal")$Temp.mod, subset(df.fit.temp.ARFlux, period=="cal")$Tsoil.C),
                                              NashSutcliffe(subset(df.fit.temp.ARFlux, period=="cal")$Temp.mod, subset(df.fit.temp.ARFlux, period=="cal")$Tsoil.C),
                                              R2(subset(df.fit.temp.ARFlux, period=="cal")$Temp.mod, subset(df.fit.temp.ARFlux, period=="cal")$Tsoil.C))
    df.fit.table.cal["Tower VWC",] <- c(RMSE(subset(df.fit.VWC.ARFlux, period=="cal")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="cal")$VWC),
                                        NRMSE(subset(df.fit.VWC.ARFlux, period=="cal")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="cal")$VWC),
                                        NashSutcliffe(subset(df.fit.VWC.ARFlux, period=="cal")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="cal")$VWC),
                                        R2(subset(df.fit.VWC.ARFlux, period=="cal")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="cal")$VWC))
    df.fit.table.cal["Thaw Depth",] <- c(RMSE(subset(df.thaw.day, period=="cal")$thaw.mod, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean),
                                         NRMSE(subset(df.thaw.day, period=="cal")$thaw.mod, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean),
                                         NashSutcliffe(subset(df.thaw.day, period=="cal")$thaw.mod, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean),
                                         R2(subset(df.thaw.day, period=="cal")$thaw.mod, subset(df.thaw.day, period=="cal")$ThawDepth.mm.mean))
    
    
    df.fit.table.val <- data.frame(RMSE=numeric(3), NRMSE=numeric(3), NSE=numeric(3), R2=numeric(3),
                                   row.names=c("Thaw Depth","Tower Soil Temp", "Tower VWC"))
    df.fit.table.val["Tower Soil Temp",] <- c(RMSE(subset(df.fit.temp.ARFlux, period=="val")$Temp.mod, subset(df.fit.temp.ARFlux, period=="val")$Tsoil.C),
                                              NRMSE(subset(df.fit.temp.ARFlux, period=="val")$Temp.mod, subset(df.fit.temp.ARFlux, period=="val")$Tsoil.C),
                                              NashSutcliffe(subset(df.fit.temp.ARFlux, period=="val")$Temp.mod, subset(df.fit.temp.ARFlux, period=="val")$Tsoil.C),
                                              R2(subset(df.fit.temp.ARFlux, period=="val")$Temp.mod, subset(df.fit.temp.ARFlux, period=="val")$Tsoil.C))
    df.fit.table.val["Tower VWC",] <- c(RMSE(subset(df.fit.VWC.ARFlux, period=="val")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="val")$VWC),
                                        NRMSE(subset(df.fit.VWC.ARFlux, period=="val")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="val")$VWC),
                                        NashSutcliffe(subset(df.fit.VWC.ARFlux, period=="val")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="val")$VWC),
                                        R2(subset(df.fit.VWC.ARFlux, period=="val")$VWC.mod, subset(df.fit.VWC.ARFlux, period=="val")$VWC))
    df.fit.table.val["Thaw Depth",] <- c(RMSE(subset(df.thaw.day, period=="val")$thaw.mod, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean),
                                         NRMSE(subset(df.thaw.day, period=="val")$thaw.mod, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean),
                                         NashSutcliffe(subset(df.thaw.day, period=="val")$thaw.mod, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean),
                                         R2(subset(df.thaw.day, period=="val")$thaw.mod, subset(df.thaw.day, period=="val")$ThawDepth.mm.mean))
    
    ## make plots
    # path to save plots
    path.plot.val.sub <- paste0(git.dir, "geotop/SoilPropertyCalibration/Plots_CalVal_Subsurface_", number, "_", fire, ".png")
    
    # plot subsurface temperature and VWC
    p.thaw.compare.ARFlux <- 
      ggplot() +
      annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
      geom_line(data=subset(df.mod.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=ThawDepth.mm), color="brown") +
      geom_segment(data=subset(df.thaw.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, xend=Date, y=ThawDepth.mm.mean-ThawDepth.mm.std, yend=ThawDepth.mm.mean+ThawDepth.mm.std)) +
      geom_point(data=subset(df.thaw.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=ThawDepth.mm.mean)) +
      scale_x_date(expand=c(0,0)) +
      scale_y_reverse(name="Tower Thaw Depth [mm]") +
      theme_bw() +
      theme(panel.grid=element_blank())
    
    p.temp.compare.ARFlux <-
      ggplot() +
      annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
      geom_hline(yintercept=0, color="gray65") +
      geom_line(data=subset(df.mod.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=Temp.ARFlux.mean), color="red") +
      geom_point(data=df.obs.ARFlux, aes(x=Date, y=Tsoil.C)) +
      scale_x_date(expand=c(0,0)) +
      scale_y_continuous(name=paste0("Tower Soil Temp [C]: ", depth.ARFlux.min, "-", depth.ARFlux.max, " mm")) +
      theme_bw() +
      theme(panel.grid=element_blank())
    
    p.VWC.compare.ARFlux <-
      ggplot() +
      annotate(geom="rect", xmin=ymd(paste0(cal, "-01-01")), xmax=ymd(paste0(cal, "-12-31")), ymin=-Inf, ymax=Inf, fill="gray90") +
      geom_hline(yintercept=0, color="gray65") +
      geom_line(data=subset(df.mod.day, year(Date)>=2008 & year(Date)<=2013), aes(x=Date, y=VWC.ARFlux.mean), color="blue") +
      geom_point(data=df.obs.ARFlux, aes(x=Date, y=VWC)) +
      scale_x_date(expand=c(0,0)) +
      scale_y_continuous(name=paste0("Tower VWC [m3/m3]: ", depths.ARFlux.VWC, " mm")) +
      theme_bw() +
      theme(panel.grid=element_blank())
    
    ggsave(path.plot.val.sub, arrangeGrob(p.thaw.compare.ARFlux, p.temp.compare.ARFlux, p.VWC.compare.ARFlux, 
                                          arrangeGrob(tableGrob(round(df.fit.table.cal, 3)), tableGrob(round(df.fit.table.val, 3)), ncol=2), 
                                          ncol=1, heights=c(1,1,1,0.75)),
           width=12, height=12, units="in")
  }  # end of fire loop
}  # end of plot loop