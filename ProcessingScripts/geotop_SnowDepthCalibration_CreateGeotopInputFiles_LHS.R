## geotop_SnowDepthCalibration_CreateGeotopInputFiles_LHS.R
#' This is intended to create a bunch of Geotop input files that will
#' be used for snow depth calibration. Parameters varied for calibration 
#' are generated using Latin Hypercube Sampling (package lhs).

rm(list=ls())

require(lhs)

# set seed
set.seed(1)

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

# path to baseline geotop input file
path.in <- paste0(git.dir, "geotop/geotop.inpts")

# path to folder to save output files
path.out <- paste0(git.dir, "geotop/SnowDepthCalibration/")

# define number of variables and number of samples
#   variables to edit:
#    -SnowCorrFactor = {0.25,1.5}
#    -ThresTempRain = {0,4.0}
#    -ThresTempSnow = {-4.0,0}
#    -IrriducibleWatSatSnow = 0.04
#    -SnowAgingCoeffVis = 0.25
#    -SnowAgingCoeffNIR = 0.45
n.var <- 3
n.sample <- 1000

# create latin hypercube sample
LHS.sample <- improvedLHS(n.sample, n.var)   # sample parameter space
df.all <- data.frame(number = sprintf("%04d", seq(1,n.sample)),
                     SnowCorrFactor = qunif(LHS.sample[,1], min=0.5, max=1.5),
                     ThresTempRain = qunif(LHS.sample[,2], min=0.0, max=4.0),
                     ThresTempSnow = qunif(LHS.sample[,3], min=-4.0, max=0.0),
                     IrriducibleWatSatSnow = 0.04,
                     SnowAgingCoeffVis = 0.25,
                     SnowAgingCoeffNIR = 0.45)

# read in template geotop.inpt file
lines.in <- readLines(path.in)

# for each variable, figure out what line it is on, and then replace
line.SnowCorrFactor <- which(startsWith(lines.in, "SnowCorrFactor"))
line.ThresTempRain <- which(startsWith(lines.in, "ThresTempRain"))
line.ThresTempSnow <- which(startsWith(lines.in, "ThresTempSnow"))
line.IrriducibleWatSatSnow <- which(startsWith(lines.in, "IrriducibleWatSatSnow"))
line.SnowAgingCoeffVis <- which(startsWith(lines.in, "SnowAgingCoeffVis"))
line.SnowAgingCoeffNIR <- which(startsWith(lines.in, "SnowAgingCoeffNIR"))

# make new geotop input file for each
for (counter in 1:n.sample){
      # make new geotop.inpts file
      lines.out <- lines.in
      lines.out[line.SnowCorrFactor] <- paste0("SnowCorrFactor = ", df.all$SnowCorrFactor[counter])
      lines.out[line.ThresTempRain] <- paste0("ThresTempRain = ", df.all$ThresTempRain[counter])
      lines.out[line.ThresTempSnow] <- paste0("ThresTempSnow = ", df.all$ThresTempSnow[counter])
      lines.out[line.IrriducibleWatSatSnow] <- paste0("IrriducibleWatSatSnow = ", df.all$IrriducibleWatSatSnow[counter])
      lines.out[line.SnowAgingCoeffVis] <- paste0("SnowAgingCoeffVis = ", df.all$SnowAgingCoeffVis[counter])
      lines.out[line.SnowAgingCoeffNIR] <- paste0("SnowAgingCoeffNIR = ", df.all$SnowAgingCoeffNIR[counter])
        
      # save geotop.inpts file
      writeLines(lines.out, paste0(path.out, df.all$number[counter], "_geotop.inpts"))
      
      # status update
      print(paste0(counter, " complete"))
}

# save data frame with all output
write.csv(df.all, paste0(path.out, "SnowDepthCalibration_Input.csv"), row.names=F)
