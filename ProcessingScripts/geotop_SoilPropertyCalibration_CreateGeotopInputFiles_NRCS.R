## geotop_SoilPropertyCalibration_CreateGeotopInputFiles_NRCS.R
#' This is intended to create a bunch of Geotop input files that will
#' be used for snow depth calibration. Parameters varied for calibration 
#' are generated using Latin Hypercube Sampling (package lhs).

rm(list=ls())

# set seed
set.seed(1)

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/ARF1D/"

# which geotop version
geo.dir <- "geotop_NRCS/"
#geo.dir <- "geotop/"

require(lhs)

# function to calculate theta (used to define wilting point and field capacity)
VG_ThetaFromHead <- function(head, ThetaR, ThetaS, alpha, n){
  # This is intended to take ThetaS [m3 m-3], ThetaR [m3 m-3], alpha [m-1], and
  # n [-] values and produce theta estimates for a vector of head [m].
  if (head >= 0){
    theta <- ThetaS
  } else {
    theta <- ThetaR + (ThetaS - ThetaR)/((1+(alpha*abs(head))^n)^(1-1/n))
  }
  
  return(theta)
}

# path to save output soil file
out.path <- paste0(git.dir, geo.dir, "SoilPropertyCalibration/")

# define soil layer properties that are not being tuned
min.Dz   <- 10     # [mm] - thickness of organic soil layers
total.Dz <- 8000   # [mm] - total soil thickness
nsoilay  <- 70     # number of soil layers

## tunable soil parameters and 'default' values upon which distributions will be centered
# organic soil - using values from Jiang et al. (2015) SI for 1st layer
#   thermal conductivity & capacity from Kurylyk et al. (2016) WRR Table A1
org.z <- 130      # [mm] - thickness of organic soil (120 cm soil moisture sensor is clearly in peat, NRCS soil pedon description max peat depth is 13 cm)
#org.Ks <- 0.17    # [mm/s] - saturated hydraulic condutivity 
org.vwc_s <- 0.55  # [m3/m3] - saturated water content
org.vwc_r <- 0.01  # [m3/m3] - residual water content
#org.VG_alpha <- 12.7*(1/1000) # [mm-1] - Van Genuchten alpha (Jiang et al. = 12.7 m-1)
#org.VG_n <- 2.00       # [-] - Van Genuchten n
#org.thermcond <- 0.25   # [W/m/K] - thermal conductivity of soil solids
org.thermcap <- 2.6E+6   # [J/m3/K] - thermal capacity of soil solids

# calculate field capacity and wilting point
org.vwc_fc <- VG_ThetaFromHead(-3.3, org.vwc_r, org.vwc_s, org.VG_alpha*1000, org.VG_n)
org.vwc_wp <- VG_ThetaFromHead(-150, org.vwc_r, org.vwc_s, org.VG_alpha*1000, org.VG_n)

# mineral soil - using values from Carsel & Parrish (1988) for silt loam based on NRCS soil pedons
#   thermal conductivity & capacity from Kurylyk et al. (2016) WRR Table A1
#min.Ks <- 0.00125  # [mm/s] - saturated hydraulic condutivity 
min.vwc_s <- 0.41  # [m3/m3] - saturated water content
min.vwc_r <- 0.01  # [m3/m3] - residual water content
#min.VG_alpha <- 2.0*(1/1000)  # [mm-1] - Van Genuchten alpha (convert from 12.7 m-1)
#min.VG_n <- 1.41              # [-] - Van Genuchten n
#min.thermcond <- 1.62         # [W/m/K] - thermal conductivity of soil solids
min.thermcap <- 2.0E+6        # [J/m3/K] - thermal capacity of soil solids

## make LHS sample
# define number of variables and number of samples
n.var <- 8
n.sample <- 350

# are you adding to an existing sample, or starting from scratch?
new.sample <- F  # T = make new sample; F = augment old sample
if (new.sample){
  # create latin hypercube sample
  LHS.sample <- improvedLHS(n.sample, n.var)   # sample parameter space
} else {
  # how many samples were in your old sample you are adding to?
  n.sample.old <- 100
  
  # recreate the old sample
  # MAKE SURE YOU ARE USING THE SAME SET.SEED AT TOP OF SCRIPT
  LHS.sample.old <- improvedLHS(n.sample.old, n.var)
  
  # augment sample
  LHS.sample <- augmentLHS(LHS.sample.old, (n.sample-n.sample.old))
}
df.all <- data.frame(number = sprintf("%04d", seq(1,n.sample)),
                     org.Ks = qunif(LHS.sample[,1], min=1.0E-3, max=1.0E+0),        # lit review
                     org.vwc_s = org.vwc_s,                                         # NRCS observations
                     org.vwc_r = org.vwc_r,                                         # NRCS observations
                     org.vwc_fc = NaN,                                              # calculate based on SWCC
                     org.vwc_wp = NaN,                                              # calculate based on SWCC
                     org.VG_alpha = qunif(LHS.sample[,2], min=1E-3, max=1E-02),     # lit review
                     org.VG_n = qunif(LHS.sample[,3], min=1.2, max=2.1),            # lit review
                     org.thermcond = qunif(LHS.sample[,4], min=2.5E-1, max=2.5E+1), # lit review
                     org.thermcap = org.thermcap,                                   # Kurylyk et al. (2016) WRR
                     min.Ks = qunif(LHS.sample[,5], min=1E-07, max=1E-04),          # +/- 1 order of magnitude
                     min.vwc_s = min.vwc_s,                                         # NRCS observations
                     min.vwc_r = min.vwc_r,                                         # NRCS observations
                     min.vwc_fc = NaN,                                              # calculate based on SWCC
                     min.vwc_wp = NaN,                                              # calculate based on SWCC
                     min.VG_alpha = qunif(LHS.sample[,6], min=1E-03, max=1E-02),    # Carsell & Parrish
                     min.VG_n = qunif(LHS.sample[,7], min=1.3, max=2.0),            # Carsell & Parrish
                     min.thermcond = qunif(LHS.sample[,8], min=1, max=2),           # lit review
                     min.thermcap = min.thermcap)                                   # Kurylyk et al. (2016) WRR

# calculate vwc_fc and vwc_wp based on soil hydraulic properties
df.all$org.vwc_fc <- VG_ThetaFromHead(-3.3, df.all$org.vwc_r, df.all$org.vwc_s, df.all$org.VG_alpha*1000, df.all$org.VG_n)
df.all$org.vwc_wp <- VG_ThetaFromHead(-150, df.all$org.vwc_r, df.all$org.vwc_s, df.all$org.VG_alpha*1000, df.all$org.VG_n)

df.all$min.vwc_fc <- VG_ThetaFromHead(-3.3, df.all$min.vwc_r, df.all$min.vwc_s, df.all$min.VG_alpha*1000, df.all$min.VG_n)
df.all$min.vwc_wp <- VG_ThetaFromHead(-150, df.all$min.vwc_r, df.all$min.vwc_s, df.all$min.VG_alpha*1000, df.all$min.VG_n)

# save output file
write.csv(df.all, paste0(out.path, "SoilPropertyCalibration_Input.csv"), quote=F, row.names=F)

## make soil property files

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
                     vwc_fc = NaN,
                     vwc_wp = NaN,
                     VG_alpha = NaN,
                     VG_n = NaN,
                     SS = 1.00E-07,       # Specific storativity - use default value,
                     thermcond = NaN,
                     thermcap = NaN)
df.out$Dz[1] <- min.Dz
df.out$z.tot[1] <- min.Dz
for (j in 1:(nsoilay-1)){
  if ((j-1)<nsoilay.org){
    # organic layers - constant thickness
    df.out$Dz[j+1] <- min.Dz
    df.out$z.tot[j+1] <- df.out$z.tot[j] + df.out$Dz[j+1]
  } else {
    # mineral soil - increasing thickness with depth
    df.out$Dz[j+1] <- min.Dz + (j-nsoilay.org)*incconst
    df.out$z.tot[j+1] <- df.out$z.tot[j] + df.out$Dz[j+1]
  }
}

# now, make a separate soil parameter file for each set of properties
for (counter in 1:n.sample){
  # copy empty df.out file
  df.counter <- df.out
  
  # fill in organic and mineral properties
  df.counter$Kh[1:nsoilay.org] <- df.all$org.Ks[counter]
  df.counter$Kv[1:nsoilay.org] <- df.all$org.Ks[counter]
  df.counter$vwc_r[1:nsoilay.org] <- df.all$org.vwc_r[counter]
  df.counter$vwc_s[1:nsoilay.org] <- df.all$org.vwc_s[counter]
  df.counter$vwc_fc[1:nsoilay.org] <- df.all$org.vwc_fc[counter]
  df.counter$vwc_wp[1:nsoilay.org] <- df.all$org.vwc_wp[counter]
  df.counter$VG_alpha[1:nsoilay.org] <- df.all$org.VG_alpha[counter]
  df.counter$VG_n[1:nsoilay.org] <- df.all$org.VG_n[counter]
  df.counter$thermcond[1:nsoilay.org] <- df.all$org.thermcond[counter]
  df.counter$thermcap[1:nsoilay.org] <- df.all$org.thermcap[counter]
  
  df.counter$Kh[(nsoilay.org+1):nsoilay] <- df.all$min.Ks[counter]
  df.counter$Kv[(nsoilay.org+1):nsoilay] <- df.all$min.Ks[counter]
  df.counter$vwc_r[(nsoilay.org+1):nsoilay] <- df.all$min.vwc_r[counter]
  df.counter$vwc_s[(nsoilay.org+1):nsoilay] <- df.all$min.vwc_s[counter]
  df.counter$vwc_fc[(nsoilay.org+1):nsoilay] <- df.all$min.vwc_fc[counter]
  df.counter$vwc_wp[(nsoilay.org+1):nsoilay] <- df.all$min.vwc_wp[counter]
  df.counter$VG_alpha[(nsoilay.org+1):nsoilay] <- df.all$min.VG_alpha[counter]
  df.counter$VG_n[(nsoilay.org+1):nsoilay] <- df.all$min.VG_n[counter]
  df.counter$thermcond[(nsoilay.org+1):nsoilay] <- df.all$min.thermcond[counter]
  df.counter$thermcap[(nsoilay.org+1):nsoilay] <- df.all$min.thermcap[counter]
  
  # save output file
  write.csv(df.counter, paste0(out.path, df.all$number[counter], "_SoilNRCS0001.txt"), quote=F, row.names=F)
  
  # status update
  print(paste0(counter, " complete"))
}
