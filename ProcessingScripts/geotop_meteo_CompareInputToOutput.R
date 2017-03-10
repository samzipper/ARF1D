## geotop_meteo_CompareInputToOutput.R
#' This script is intended to compare input meteo data and GEOtop
#' output to make sure they are the same.

rm(list=ls())

# git directory for relative paths
#git.dir <- "C:/Users/Sam/WorkGits/Permafrost/ARF1D/"
git.dir <- "C:/Users/Sam/WorkGits/ARF1D/"

require(ggplot2)
require(gridExtra)
require(lubridate)

# path to input and output, and where to save plots
#in.path <- paste0(git.dir, "geotop/meteo/meteoNARRhourlyDynamicWithSpinUp0001.txt")
#in.path <- paste0(git.dir, "geotop/meteo/meteoNARR3hourlyWithSpinUp0001.txt")
in.path <- paste0(git.dir, "geotop/meteo/meteoNARRdailyWithSpinUp0001.txt")
out.path <- paste0(git.dir, "geotop/output-tabs/point0001.txt")
plot.path <- paste0(git.dir, "geotop/output-plots/meteo_CompareInputToOutput_HourlyIn_4hrRun_DailyOut.png")

#in.path <- paste0("C:/Users/Sam/src/geotop/tests/1D/InfiltrationTrench/meteo/meteotrenchhour0001.txt")
#out.path <- paste0("C:/Users/Sam/src/geotop/tests/1D/InfiltrationTrench/output-tabs/point0001.txt")
#plot.path <- paste0("C:/Users/Sam/src/geotop/tests/1D/InfiltrationTrench/output-tabs/meteo_CompareInputToOutput.png")

# year to compare
yr <- 2010

# read input and output
df.in <- read.csv(in.path)
df.out <- read.csv(out.path)

# set NaNs
df.in[df.in==-9999] <- NaN
df.out[df.out==-9999] <- NaN

# make date column
df.in$Date <- dmy_hm(df.in$POSIX)
df.out$Date <- dmy_hm(df.out$Date12.DDMMYYYYhhmm.)

df.in$Year <- year(df.in$Date)
df.in$DOY <- yday(df.in$Date)+1

df.out$Year <- year(df.out$Date)
df.out$DOY <- yday(df.out$Date)

# trim to comparison year only
df.in <- subset(df.in, Year==yr)
df.out <- subset(df.out, Year==yr)

# subset to matching dates
df.in <- subset(df.in, DOY %in% df.out$DOY)
df.out <- subset(df.out, DOY %in% df.in$DOY)

# combine
df <- data.frame(Year = df.in$Year,
                 DOY = df.in$DOY,
                 Date = df.in$Date,
                 prec.in = df.in$Iprec*24,
                 #prec.in = df.in$Iprec,
                 prec.out = df.out$Psnow_over_canopy.mm.+df.out$Prain_over_canopy.mm.,
                 wind.in = df.in$WindSp,
                 wind.out = df.out$Wind_speed.m.s.,
                 Tair.in = df.in$AirT,
                 Tair.out = df.out$Tair.C.,
                 RH.in = df.in$RH,
                 RH.out = df.out$Relative_Humidity...*100,
                 SWin.in = df.in$Swglob,
                 SWin.out = df.out$SWin.W.m2.)

df.long <- rbind(data.frame(Year = df$Year,
                            DOY = df$DOY,
                            Date = df$Date,
                            prec = df$prec.in,
                            wind = df$wind.in,
                            Tair = df$Tair.in,
                            RH = df$RH.in,
                            SWin = df$SWin.in,
                            source = "input"),
                 data.frame(Year = df$Year,
                            DOY = df$DOY,
                            Date = df$Date,
                            prec = df$prec.out,
                            wind = df$wind.out,
                            Tair = df$Tair.out,
                            RH = df$RH.out,
                            SWin = df$SWin.out,
                            source = "output"))

# make plots
p.prec.scatter <-
  ggplot(df, aes(x=prec.in, y=prec.out)) +
  geom_point(shape=21) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  labs(title="Precipitation [mm/day]") +
  scale_x_continuous(name="Meteo Input File") +
  scale_y_continuous(name="Point Output File") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.prec.time <- 
  ggplot(df.long, aes(x=Date, y=prec, color=source)) +
  geom_line() +
  geom_hline(yintercept=0, color="gray65") +
  labs(title="Precipitation [mm/day]") +
#  scale_x_continuous(name="Date") +
  scale_y_continuous(name="Quantity") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")

p.Tair.scatter <-
  ggplot(df, aes(x=Tair.in, y=Tair.out)) +
  geom_point(shape=21) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  labs(title="Air Temperature [C]") +
  scale_x_continuous(name="Meteo Input File") +
  scale_y_continuous(name="Point Output File") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.Tair.time <- 
  ggplot(df.long, aes(x=Date, y=Tair, color=source)) +
  geom_line() +
  geom_hline(yintercept=0, color="gray65") +
  labs(title="Air Temperature [C]") +
#  scale_x_continuous(name="DOY") +
  scale_y_continuous(name="Quantity") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")

p.wind.scatter <-
  ggplot(df, aes(x=wind.in, y=wind.out)) +
  geom_point(shape=21) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  labs(title="Wind Speed [m/s]") +
  scale_x_continuous(name="Meteo Input File") +
  scale_y_continuous(name="Point Output File") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.wind.time <- 
  ggplot(df.long, aes(x=Date, y=wind, color=source)) +
  geom_line() +
  geom_hline(yintercept=0, color="gray65") +
  labs(title="Wind Speed [m/s]") +
#  scale_x_continuous(name="DOY") +
  scale_y_continuous(name="Quantity") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")

p.RH.scatter <-
  ggplot(df, aes(x=RH.in, y=RH.out)) +
  geom_point(shape=21) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  labs(title="Relative Humidity [%]") +
  scale_x_continuous(name="Meteo Input File") +
  scale_y_continuous(name="Point Output File") +
  theme_bw() +
  theme(panel.grid=element_blank())

p.RH.time <- 
  ggplot(df.long, aes(x=Date, y=RH, color=source)) +
  geom_line() +
  labs(title="Relative Humidity [%]") +
#  scale_x_continuous(name="DOY") +
  scale_y_continuous(name="Quantity") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")

p.SWin.scatter <-
  ggplot(df, aes(x=SWin.in, y=SWin.out)) +
  geom_point(shape=21) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  labs(title="Incoming SW Radiation [W/m2]") +
  scale_x_continuous(name="Meteo Input File") +
  scale_y_continuous(name="Point Output File") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")

p.SWin.time <- 
  ggplot(df.long, aes(x=Date, y=SWin, color=source)) +
  geom_line() +
  labs(title="Incoming SW Radiation [W/m2]") +
#  scale_x_continuous(name="DOY") +
  scale_y_continuous(name="Quantity") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position="bottom")

# save plot
ggsave(plot.path, arrangeGrob(p.Tair.time, p.Tair.scatter,
                              p.SWin.time, p.SWin.scatter,
                              p.prec.time, p.prec.scatter,
                              p.wind.time, p.wind.scatter,
                              p.RH.time, p.RH.scatter, 
                              ncol=2),
       width=16, height=12, units="in")
