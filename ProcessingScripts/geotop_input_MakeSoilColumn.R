## geotop_input_MakeSoilColumn.R
#' This script is intended to  create a soil input file for simulations
#' in GEOtop. The soil column will have two soil types: organic and mineral 
#' soil. Adjustable parameters will include the thickness and hydraulic
#' retention properties of each soil type.
#' 
#' The soil columns is structured as a series of uniform thickness layers for  
#' organic soil interval, and then linearly increasing thickness for the
#' mineral soil interval.

rm(list=ls())

# git directory for relative paths
#git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"
git.dir <- "C:/Users/Sam/WorkGits/ARF1D/"

# path to save output soil file
out.path <- paste0(git.dir, "geotop/soil/soilARF0001.txt")

# define soil layer properties
min.Dz <- 10       # [mm] - thickness of organic soil layers
total.Dz <- 10000  # [mm] - total soil thickness
nsoilay <- 50     # number of soil layers

## tunable soil parameters
# organic soil - using values from Jiang et al. (2015) SI for 1st layer
#   thermal conductivity & capacity from Kurylyk et al. (2016) WRR Table A1
org.z <- 174      # [mm] - thickness of organic soil (Jiang et al. 2015, Table 1, thickness at unburned site = 17.4 cm +/- 2.1 cm)
org.Ks <- 0.17    # [mm/s] - saturated hydraulic condutivity 
org.vwc_s <- 0.84  # [m3/m3] - saturated water content
org.vwc_r <- 0.04  # [m3/m3] - residual water content
org.VG_alpha <- 1.0*(1/1000) # [mm-1] - Van Genuchten alpha (Jiang et al. = 12.7 m-1)
org.VG_n <- 1.30       # [-] - Van Genuchten n
org.thermcond <- 0.25   # [W/m/K] - thermal conductivity of soil solids
org.thermcap <- 2.6E+6   # [J/m3/K] - thermal capacity of soil solids

# mineral soil - using values from Jiang et al. (2015) SI for 3rd layer
#   thermal conductivity & capacity from Kurylyk et al. (2016) WRR Table A1
min.Ks <- 0.021   # [mm/s] - saturated hydraulic condutivity 
min.vwc_s <- 0.56 # [m3/m3] - saturated water content
min.vwc_r <- 0.07 # [m3/m3] - residual water content
min.VG_alpha <- 2.41*(1/1000) # [mm-1] - Van Genuchten alpha (convert from 12.7 m-1)
min.VG_n <- 1.33       # [-] - Van Genuchten n
min.thermcond <- 1.62   # [W/m/K] - thermal conductivity of soil solids
min.thermcap <- 2.0E+6   # [J/m3/K] - thermal capacity of soil solids

## figure out number of organic and mineral layers based on thicknesses
nsoilay.org <- round(org.z/min.Dz)
nsoilay.min <- nsoilay - nsoilay.org

## figure out increment to increase mineral soil layer thickness with depth
# coefficient for incrementing
incrcoeff <- 0.0
for (j in 1:(nsoilay.min-1)){
  # figure out total increment coefficient
  incrcoeff <- incrcoeff + j
}

# incrementing constant
incconst <- ((total.Dz-nsoilay.org*min.Dz) - (min.Dz*nsoilay.min))/incrcoeff

## build soil layers
df.out <- data.frame(Dz = numeric(length=nsoilay),
                     z.tot = NaN,
                     Kh = NaN,
                     Kv = NaN,
                     vwc_r = NaN,
                     vwc_s = NaN,
                     VG_alpha = NaN,
                     VG_n = NaN,
                     SS = NaN,
                     thermcond = NaN,
                     thermcap = NaN)
df.out$Dz[1] <- min.Dz
df.out$z.tot[1] <- min.Dz
df.out$Kh[1] <- org.Ks
df.out$Kv[1] <- org.Ks
df.out$vwc_r[1] <- org.vwc_r
df.out$vwc_s[1] <- org.vwc_s
df.out$VG_alpha[1] <- org.VG_alpha
df.out$VG_n[1] <- org.VG_n
df.out$thermcond[1] <- org.thermcond
df.out$thermcap[1] <- org.thermcap
for (j in 1:(nsoilay-1)){
  if ((j-1)<nsoilay.org){
    # organic layers
    df.out$Dz[j+1] <- min.Dz
    df.out$z.tot[j+1] <- df.out$z.tot[j] + df.out$Dz[j+1]
    df.out$Kh[j+1] <- org.Ks
    df.out$Kv[j+1] <- org.Ks
    df.out$vwc_r[j+1] <- org.vwc_r
    df.out$vwc_s[j+1] <- org.vwc_s
    df.out$VG_alpha[j+1] <- org.VG_alpha
    df.out$VG_n[j+1] <- org.VG_n
    df.out$thermcond[j+1] <- org.thermcond
    df.out$thermcap[j+1] <- org.thermcap
  } else {
    # mineral soil
    df.out$Dz[j+1] <- min.Dz + (j-nsoilay.org)*incconst
    df.out$z.tot[j+1] <- df.out$z.tot[j] + df.out$Dz[j+1]
    df.out$Kh[j+1] <- min.Ks
    df.out$Kv[j+1] <- min.Ks
    df.out$vwc_r[j+1] <- min.vwc_r
    df.out$vwc_s[j+1] <- min.vwc_s
    df.out$VG_alpha[j+1] <- min.VG_alpha
    df.out$VG_n[j+1] <- min.VG_n
    df.out$thermcond[j+1] <- min.thermcond
    df.out$thermcap[j+1] <- min.thermcap
  }
}

# anything that does not vary based on depth
df.out$SS <- 1.00E-07       # Specific storativity - use default value

# save output file
write.csv(df.out, out.path, quote=F, row.names=F)
