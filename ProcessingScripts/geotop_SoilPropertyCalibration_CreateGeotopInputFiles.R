## geotop_SoilPropertyCalibration_CreateGeotopInputFiles.R
#' This is intended to create a bunch of Geotop input files that will
#' be used for snow depth calibration. Parameters varied for calibration 
#' are generated using Latin Hypercube Sampling (package lhs).

rm(list=ls())

# set seed
set.seed(1)

# git directory for relative paths
git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"

require(lhs)

# path to save output soil file
out.path <- paste0(git.dir, "geotop/SoilPropertyCalibration/")

# define soil layer properties that are not being tuned
min.Dz <- 5        # [mm] - thickness of organic soil layers
total.Dz <- 10000  # [mm] - total soil thickness
nsoilay <- 100     # number of soil layers

## tunable soil parameters and 'default' values upon which distributions will be centered
## (14 tunable parameters total)
# organic soil - using values from Jiang et al. (2015) SI for 1st layer
#   thermal conductivity & capacity from Kurylyk et al. (2016) WRR Table A1
org.Ks <- 0.17    # [mm/s] - saturated hydraulic condutivity 
org.vwc_s <- 0.5  # [m3/m3] - saturated water content
org.vwc_r <- 0.1  # [m3/m3] - residual water content
org.VG_alpha <- 12.7*(1/1000) # [mm-1] - Van Genuchten alpha (convert from 12.7 m-1)
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

## make LHS sample
# define number of variables and number of samples
n.var <- 14
n.sample <- 1000

# create latin hypercube sample
LHS.sample <- improvedLHS(n.sample, n.var)   # sample parameter space
df.all <- data.frame(number = sprintf("%04d", seq(1,n.sample)),
                     org.Ks = qunif(LHS.sample[,1], min=org.Ks/10, max=org.Ks*10),                             # +/- 1 order of magnitude
                     org.vwc_s = qunif(LHS.sample[,2], min=org.vwc_s-0.1, max=org.vwc_s+0.1),                  # +/- 0.1
                     org.vwc_r = qunif(LHS.sample[,3], min=max(c(0.01, org.vwc_r-0.05)), max=org.vwc_r+0.05),  # +/- 0.05, but can't be <0.01
                     org.VG_alpha = qunif(LHS.sample[,4], min=org.VG_alpha*0.9, max=org.VG_alpha*1.1),         # +/- 10%
                     org.VG_n = qunif(LHS.sample[,5], min=org.VG_n*0.9, max=org.VG_n*1.1),                     # +/- 10%
                     org.thermcond = qunif(LHS.sample[,6], min=org.thermcond*0.75, max=org.thermcond*1.25),    # +/- 25%
                     org.thermcap = qunif(LHS.sample[,7], min=org.thermcap*0.75, max=org.thermcap*1.25),       # +/- 25%
                     min.Ks = qunif(LHS.sample[,8], min=min.Ks/10, max=min.Ks*10),                             # +/- 1 order of magnitude
                     min.vwc_s = qunif(LHS.sample[,9], min=min.vwc_s-0.1, max=min.vwc_s+0.1),                  # +/- 0.1
                     min.vwc_r = qunif(LHS.sample[,10], min=max(c(0.01, min.vwc_r-0.05)), max=min.vwc_r+0.05), # +/- 0.05, but can't be <0.01
                     min.VG_alpha = qunif(LHS.sample[,11], min=min.VG_alpha*0.9, max=min.VG_alpha*1.1),        # +/- 10%
                     min.VG_n = qunif(LHS.sample[,12], min=min.VG_n*0.9, max=min.VG_n*1.1),                    # +/- 10%
                     min.thermcond = qunif(LHS.sample[,13], min=min.thermcond*0.75, max=min.thermcond*1.25),   # +/- 25%
                     min.thermcap = qunif(LHS.sample[,14], min=min.thermcap*0.75, max=min.thermcap*1.25))      # +/- 25%

# save output file
write.csv(df.all, paste0(out.path, "SoilPropertyCalibration_Input.csv"), quote=F, row.names=F)

## cycle through burn/unburned sites and make soil property files for each
for (fire in c("Unburned", "Moderate", "Severe")){
  # Define organic layer thickness [mm], based on Jiang et al. 2015, Table 1
  if (fire=="Unburned"){
    org.z <- 174
  } else if (fire=="Moderate"){
    org.z <- 154
  } else {
    org.z <- 116
  }
  
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
    df.counter$VG_alpha[1:nsoilay.org] <- df.all$org.VG_alpha[counter]
    df.counter$VG_n[1:nsoilay.org] <- df.all$org.VG_n[counter]
    df.counter$thermcond[1:nsoilay.org] <- df.all$org.thermcond[counter]
    df.counter$thermcap[1:nsoilay.org] <- df.all$org.thermcap[counter]
    
    df.counter$Kh[(nsoilay.org+1):nsoilay] <- df.all$min.Ks[counter]
    df.counter$Kv[(nsoilay.org+1):nsoilay] <- df.all$min.Ks[counter]
    df.counter$vwc_r[(nsoilay.org+1):nsoilay] <- df.all$min.vwc_r[counter]
    df.counter$vwc_s[(nsoilay.org+1):nsoilay] <- df.all$min.vwc_s[counter]
    df.counter$VG_alpha[(nsoilay.org+1):nsoilay] <- df.all$min.VG_alpha[counter]
    df.counter$VG_n[(nsoilay.org+1):nsoilay] <- df.all$min.VG_n[counter]
    df.counter$thermcond[(nsoilay.org+1):nsoilay] <- df.all$min.thermcond[counter]
    df.counter$thermcap[(nsoilay.org+1):nsoilay] <- df.all$min.thermcap[counter]
    
    # save output file
    write.csv(df.counter, paste0(out.path, df.all$number[counter], "_", fire, "_SoilARF0001.txt"), quote=F, row.names=F)
    
    # status update
    print(paste0(fire, " ", counter, " complete"))
  }

}
