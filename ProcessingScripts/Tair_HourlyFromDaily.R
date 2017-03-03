Tair_HourlyFromDaily <- function(solar.time, Tmax, Tmin, Tmax.yesterday, Tmin.tomorrow){
  # this function applies equations 2.2 & 2.3 from Campbell & Norman "Environmental Biophysics" to 
  # calculate hourly temperature based on a diurnal temperature function, solar time, and Tmax/Tmin 
  # from the current day, the previous day, and the following day.
  #
  # solar.time should be in hours, where 12=solar noon
  # Tmax and Tmin should be in degrees C
  #
  # all inputs should be vectors
  
  # make empty vector for output
  Tair <- numeric(length=length(solar.time))
  
  # diurnal dimensionless 0-1 parameter for scaling temperature amplitude
  diurnal <- 0.44 - 0.46*sin((pi/12)*solar.time+0.9) + 0.11*sin(2*(pi/12)*solar.time+0.9)  # Campbell & Norman eq. 2.2
  
  # calculate temperature based on time (Eq 2.3)
  Tair[solar.time<=5] <- Tmax.yesterday[solar.time<=5]*diurnal[solar.time<=5] + Tmin[solar.time<=5]*(1-diurnal[solar.time<=5])
  Tair[solar.time>5 & solar.time<=14] <- 
    Tmax[solar.time>5 & solar.time<=14]*diurnal[solar.time>5 & solar.time<=14] + 
    Tmin[solar.time>5 & solar.time<=14]*(1-diurnal[solar.time>5 & solar.time<=14])
  Tair[solar.time>14] <- Tmax[solar.time>14]*diurnal[solar.time>14] + Tmin.tomorrow[solar.time>14]*(1-diurnal[solar.time>14])
  
  return(Tair)
}