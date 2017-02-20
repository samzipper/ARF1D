## geotop_SnowDepthCalibration_CreateGeotopInputFiles.R
#' This is intended to create a bunch of Geotop input files that will
#' be used for snow depth calibration.

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

# path to baseline geotop input file
path.in <- paste0(git.dir, "geotop/geotop.inpts")

# path to folder to save output files
path.out <- paste0(git.dir, "geotop/SnowDepthCalibration/")

# variables to edit:
#  -SnowCorrFactor = {0.5,1.5}
#  -ThresTempRain = {0,4.0}
#  -ThresTempSnow = {-4.0,0}
SnowCorrFactor.all <- seq(0.5,1.5,0.1)
ThresTempRain.all <- seq(0,4,0.25)
ThresTempSnow.all <- seq(-4,0,0.25)
n.options <- length(SnowCorrFactor.all)*length(ThresTempRain.all)*length(ThresTempSnow.all)

# read in template geotop.inpt file
lines.in <- readLines(path.in)

# for each variable, figure out what line it is on, and then replace
line.SnowCorrFactor <- which(startsWith(lines.in, "SnowCorrFactor"))
line.ThresTempRain <- which(startsWith(lines.in, "ThresTempRain"))
line.ThresTempSnow <- which(startsWith(lines.in, "ThresTempSnow"))

# make data frame to hold output
df.all <- data.frame(number = sprintf("%04d", seq(1,n.options)),
                     SnowCorrFactor = numeric(n.options),
                     ThresTempRain = numeric(n.options),
                     ThresTempSnow = numeric(n.options))

# make new geotop input file for each
counter <- 0
for (SnowCorrFactor in SnowCorrFactor.all){
  for (ThresTempRain in ThresTempRain.all){
    for (ThresTempSnow in ThresTempSnow.all){
      # iterate counter
      counter <- counter+1
      
      # fill in data frame
      df.all$SnowCorrFactor[counter] <- SnowCorrFactor
      df.all$ThresTempRain[counter] <- ThresTempRain
      df.all$ThresTempSnow[counter] <- ThresTempSnow
      
      # make new geotop.inpts file
      lines.out <- lines.in
      lines.out[line.SnowCorrFactor] <- paste0("SnowCorrFactor = ", SnowCorrFactor)
      lines.out[line.ThresTempRain] <- paste0("ThresTempRain = ", ThresTempRain)
      lines.out[line.ThresTempSnow] <- paste0("ThresTempSnow = ", ThresTempSnow)
        
      # save geotop.inpts file
      writeLines(lines.out, paste0(path.out, df.all$number[counter], "_geotop.inpts"))
      
      # status update
      print(paste0(counter, " complete"))
    }
  }
}

# save data frame with all output
write.csv(df.all, paste0(path.out, "SnowDepthCalibration_Input.csv"), row.names=F)
