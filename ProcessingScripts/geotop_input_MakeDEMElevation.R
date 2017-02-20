## geotop_input_MakeDEMElevation.R
#' This script is intended to read in a DEM, convert all the data points to
#' a set number, and save it as a new DEM.
#' 

rm(list=ls())

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

require(raster)

# path to DEM (will be saved as same name)
DEM.path <- paste0(git.dir, "geotop/dem_mybasin_2.asc")

# elevation of all point with data
elev <- 760     # [m] - this is the Toolik Field Station elevation

# read in DEM
DEM <- raster(DEM.path)

# replace data
DEM[is.finite(DEM)] <- elev

# save output raster
writeRaster(DEM, filename=DEM.path, format="ascii", overwrite=T)
