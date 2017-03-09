setwd("C:/Users/Sam/WorkGits/Permafrost/ARF1D/data/ARFlux/Raw-Ameriflux")

df <- read.csv("AMF_US-An1_BASE_HH_1-1.csv", skip=2)
df <- read.csv("AMF_US-An2_BASE_HH_1-1.csv", skip=2)
df <- read.csv("AMF_US-An3_BASE_HH_1-1.csv", skip=2)

df$Date <- ymd_hm(df$TIMESTAMP_END)

df[df==-9999] <- NaN

p.VWC <-
  ggplot(df, aes(x=Date, y=SWC_1)) +
  geom_point()