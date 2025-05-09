library(nlme)
library(tidyverse)
library(readr)

data_zone<-read_csv("sixres_zone.csv")

data_zone$latlon<-paste(data_zone$lat, data_zone$lon)
data_zone[-which(duplicated(data_zone$latlon) == TRUE),]

model<-lme(
  fixed = turb ~zone,
  random = ~1|system,
  data=data_zone,
  correlation = corExp(form= ~lon+lat),
  method="REML"
)


# model <- gls(data_zone ~ time, correlation = corARMA(p = 1), data = data.frame(time = time))

