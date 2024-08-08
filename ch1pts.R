# A script for with working code for chapter 1 data analysis
# This combines from scripts to include working with FLAMe data, CRASR samples and ee.points transect queries 

# Note that everything with GEE TIF output is in its own script "dwl6lakes.R"
# Processing ACOLITE output is in "acolite6lakes.R"

# clear workspace
rm(list=ls())

# 1. load essential packages

library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(OpenStreetMap)
library(colorscience)
library(sf)
library(raster)
require(rasterVis)
require(RColorBrewer)
require(terrainr)
require(rstac)
require(terra)
library(nlme)
library(anytime)

# 2. Reading in data

# Reading the combined FLAMe and Sentinel-2 data for pts along boat path with essential data variables:
# lat, lon, DWL, predicted SD, NDTI, turbidity
# Criteria for this data was speed >0, unique lat/lon, and "no flag" (as specified by AH paper)
# Filters dropped obs from 48033 to 28847; see "s2flm_pts.R" for making of sausage
s2flame_znpts<-read_csv("six4m_zone.csv")


# Reading in data for distance from dam transects queried in ee.points Shiny
# inlcudes all output data plus NDTI and normalized DFD
dfd_transect<- read_csv("dfd_transect.csv")

# Reading in sample to sensor validation data (using method from AH paper) 
# refer to "station_YSI.R" for making of sausage
sample_ysi<-read_csv("sensor_sample_valmeans.csv")


# 4. Plotting data

# plotting NDTI x norm DFD
line_ndti<- dfd_transect %>% ggplot(aes(norm_dist, ndti, color = system)) + geom_line(size=1.25) + 
  xlab("Normalized Distance from Dam") + ylab("Normalized Difference Turbidity Index") +  
  theme_bw()+scale_color_manual(values = viridis::viridis(6, option = "C"),
                                labels=c("bonham" = "Bonham", "waco"="Waco","brownwood"="Brownwood",
                                         "ohivie"="O.H. Ivie","redbluff"="Red Bluff","arrowhead"="Arrowhead"),
                                name = "System")+
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

# plotting DWL x norm DFD
line_dwl_p<- dfd_transect %>% ggplot(aes(norm_dist, dwl, color = system)) + geom_line(size=1.25) + 
  xlab("Normalized Distance from Dam") + ylab("Dominant Wavelength") + theme_bw() +
  scale_color_manual(values = viridis::viridis(6, option = "C"),
                     labels=c("bonham" = "Bonham", "waco"="Waco","brownwood"="Brownwood",
                              "ohivie"="O.H. Ivie","redbluff"="Red Bluff","arrowhead"="Arrowhead"),
                     name = "System")+
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

# plotting relationship between grab sample data and mean sensor values
lmt_valmeans<-lm(turb_ysi_m~turb_lab,data=sample_ysi) 
summary(lmt_valmeans) # r2 = 0.9893

# plotting sensor~sample for all 6 lakes - same lm as above
all6t_r2 <-ggplot(sample_ysi,aes(y=turb_ysi_m,x=turb_lab)) + 
  geom_point(aes(color = system)) + 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+ theme_bw()
# code for plots of individual lakes and their stations + all combined are found in "station_ysi.R" 

# plotting SDD to predicted SDD

# NDTI~turbidity from boat path points
logEstimate <- lm(ndti~log(turb),data=s2flame_znpts)
xvec <- seq(min(s2flame_znpts$turb),max(s2flame_znpts$turb),length=1000)
logpred <- predict(logEstimate,newdata=data.frame(turb=xvec))
pred <- data.frame(x = xvec, y = logpred)

# simple plot showing R2 = 0.59 for all lakes
turbndti_p<-ggplot(s2flame_znpts,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

# plot done in the style of AH paper
turb_ndti_plot<-ggplot(s2flame_znpts,aes(x=turb,y=ndti)) +
  geom_point(size=0.1,col=gray(0.3))+
  geom_smooth(method="lm",lty=5,col="black",se=FALSE,lwd=1)+
  geom_line(data = pred, aes(x=x, y=y),lty=1,col="black",lwd=1)+
  #xlim(c(0,99))+
  xlab("In-lake Turbidity (NTU)")+
  ylab("Satellite Turbidity (NDTI)")+
  theme_bw()

# 4a. Summary statistics

# defining confidence level
confidence_level <- 0.95
alpha <- 1-confidence_level

# calculating means, CI, and range for variables by zone, can also specify group_by(system)

zonestats_turb<- s2flame_znpts %>% group_by(zone) %>% 
  summarise(n=n(),mean=mean(turb),sd=sd(turb),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(turb), max=max(turb),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))


zonestats_secchi<- s2flame_znpts %>% group_by(zone) %>% 
  summarise(n=n(),mean=mean(secchi),sd=sd(secchi),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(secchi), max=max(secchi),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

zonestats_ndti<- s2flame_znpts %>% group_by(zone) %>% 
  summarise(n=n(),mean=mean(ndti),sd=sd(ndti),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(dwl), max=max(dwl),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

zonestats_dwl<- s2flame_znpts %>% group_by(zone) %>% 
  summarise(n=n(),mean=mean(dwl),sd=sd(dwl),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(dwl), max=max(dwl),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

# note that results for DWL along boat path are not consistent with full system analysis
# same summary stats were calculated for the zone designated water pixel data from GEE output
## see "dwl5lakes.R" for processing of GEE output - histogram creation, zone designation, and stats ##

# 4b. AOV and Tukeys HSD post-hoc

# using the subset method (10% of data) applying SP code - full can be found in "s2flm_pts.R"

# concatenate lat/lon and remove duplicated lat/lon pairs
s2flame_znpts$latlon<-paste(s2flame_znpts$lat,s2flame_znpts$lon)
s2flame_znpts<-s2flame_znpts[-which(duplicated(s2flame_znpts$latlon)==TRUE),]
s2flame_znpts<-data.frame(s2flame_znpts)
s2flame_znpts$index<-1:length(s2flame_znpts[,1])

# use a random subset of data to avoid autocorrelation and speed up processing
# see Loken et al. 2019 https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JG005186
set.seed(1) # makes the subsetting reproducible
sampled <- sample(1:length(s2flame_znpts[,1]),size=2884) #2528) #200) # adjust subset size as needed #total is 28728
data_sampled<-s2flame_znpts[sampled,]

# fit new anova on random subset for turb
aov_turb <- aov(turb~zone+system,data_sampled)
posthoc <- TukeyHSD(aov_turb)
posthoc

# do diagnostics
plot(aov_turb,which=1)
plot(aov_turb,which=2)
plot(aov_turb,which=3)

# buiid data frame with residuals and fitted, to check model diagnostics
resid <- as.numeric(aov_turb$residuals)
resid_st <- as.numeric(rstudent(aov_turb)) #
fitted <-as.numeric(aov_turb$fitted.values)
value <-as.numeric(aov_turb$model$turb)
aov1_df<- data.frame(zone=aov_turb$model$zone,
                     system=aov_turb$model$system,
                     value,
                     resid,resid_st,fitted)

# this is a ggplot version of the above residual vs fitted, plot(aov1,which=1)
resid_vs_fitted<-aov1_df %>%
  ggplot(aes(x=fitted,y=resid))+
  geom_point()+
  theme_bw()
resid_vs_fitted

# do more diagnostics - acf and pacf of the fitted and residuals
png("aov1_acf_fitted.png")
acf_fitted<-acf(aov1_df$fitted)
dev.off()
png("aov1_pacf_fitted.png")
pacf_fitted<-pacf(aov1_df$fitted)
dev.off()
png("aov1_acf_resid.png")
acf_resid<-acf(aov1_df$resid)
dev.off()
png("aov1_pacf_resid.png")
pacf_resid<-pacf(aov1_df$resid)
dev.off()

# fit anova on random subset for predicted Secchi
aov_secchi <- aov(secchi~zone+system,data_sampled)
summary(aov_secchi)
secchi_posthoc <- TukeyHSD(aov_secchi)
secchi_posthoc
# buiid data frame with residuals and fitted, to check model diagnostics
resid <- as.numeric(aov_secchi$residuals)
resid_st <- as.numeric(rstudent(aov_secchi)) #
fitted <-as.numeric(aov_secchi$fitted.values)
value <-as.numeric(aov_secchi$model$secchi)
aov_secchi_df<- data.frame(zone=aov_secchi$model$zone,
                           system=aov_secchi$model$system,
                           value,
                           resid,resid_st,fitted)

# this is a ggplot version of the above residual vs fitted, plot(aov1,which=1)
resid_vs_fitted<-aov_secchi_df %>%
  ggplot(aes(x=fitted,y=resid))+
  geom_point()+
  theme_bw()
resid_vs_fitted

# do more diagnostics - acf and pacf of the fitted and residuals

acf_fitted<-acf(aov_secchi_df$fitted)
pacf_fitted<-pacf(aov_secchi_df$fitted)
acf_resid<-acf(aov_secchi_df$resid)
pacf_resid<-pacf(aov_secchi_df$resid)

# fit anova on random subset for predicted DWL
aov_dwl <- aov(dwl~zone+system,data_sampled)
summary(aov_dwl)
dwl_posthoc <- TukeyHSD(aov_dwl)
dwl_posthoc

# buiid data frame with residuals and fitted, to check model diagnostics
resid <- as.numeric(aov_dwl$residuals)
resid_st <- as.numeric(rstudent(aov_dwl)) #
fitted <-as.numeric(aov_dwl$fitted.values)
value <-as.numeric(aov_dwl$model$dwl)
aov_dwl_df<- data.frame(zone=aov_dwl$model$zone,
                        system=aov_dwl$model$system,
                        value,
                        resid,resid_st,fitted)

# this is a ggplot version of the above residual vs fitted, plot(aov1,which=1)
dwl_resid_vs_fitted<-aov_dwl_df %>%
  ggplot(aes(x=fitted,y=resid))+
  geom_point()+
  theme_bw()
dwl_resid_vs_fitted

# do more diagnostics - acf and pacf of the fitted and residuals

acf_fitted<-acf(aov_dwl_df$fitted)
pacf_fitted<-pacf(aov_dwl_df$fitted)
acf_resid<-acf(aov_dwl_df$resid)
pacf_resid<-pacf(aov_dwl_df$resid)


aov_ndti <- aov(ndti~zone+system,data_sampled)
summary(aov_ndti)
ndti_posthoc <- TukeyHSD(aov_ndti)
ndti_posthoc

# buiid data frame with residuals and fitted, to check model diagnostics
resid <- as.numeric(aov_ndti$residuals)
resid_st <- as.numeric(rstudent(aov_ndti)) #
fitted <-as.numeric(aov_ndti$fitted.values)
value <-as.numeric(aov_ndti$model$ndti )
aov_ndti_df<- data.frame(zone=aov_ndti$model$zone,
                         system=aov_ndti$model$system,
                         value,
                         resid,resid_st,fitted)

# this is a ggplot version of the above residual vs fitted, plot(aov1,which=1)
ndti_resid_vs_fitted<-aov_ndti_df %>%
  ggplot(aes(x=fitted,y=resid))+
  geom_point()+
  theme_bw()
ndti_resid_vs_fitted

# do more diagnostics - acf and pacf of the fitted and residuals

acf_fitted<-acf(aov_ndti_df$fitted)
pacf_fitted<-pacf(aov_ndti_df$fitted)
acf_resid<-acf(aov_ndti_df$resid)
pacf_resid<-pacf(aov_ndti_df$resid)


################################################################################
# Just to show how some things were determined

# 1. predicted Secchi
# adding NLA derived TSI to each flame point
# reading in data (all NLA reservoirs for 2012 + 2017 and CRASR sample data)
nla_crasr<-read_csv("nla_crasr.csv")
nla_crasr <- nla_crasr%>% dplyr::rename(secchi = `secchi (m)`, turb = `turb (ntu)`)
# lm for log(secchi)~log*turb
tsdlm<-lm(log(secchi)~log(turb), data = nla_crasr)
# predicting secchi values from EXO turbidity data
all_bind$secchi <- exp(predict(tsdlm, all_bind))
# calculating TSI from newly derived secchi
all_bind$tsi_sd <- 10 * (6-log(all_bind$secchi)/log(2))

# 2. Chromaticity calculation for Dominant Wavelength 
#Convert R,G, and B spectral reflectance to dwl based on CIE chromaticity color space, see Wang et al 2015.
fui.hue <- function(R, G, B) {
  
  require(colorscience)
  # chromaticity.diagram.color.fill()
  Xi <- 2.7689*R + 1.7517*G + 1.1302*B
  Yi <- 1.0000*R + 4.5907*G + 0.0601*B
  Zi <- 0.0565*G + 5.5943*B
  
  # calculate coordinates on chromaticity diagram
  x <-  Xi / (Xi + Yi +  Zi)
  y <-  Yi / (Xi + Yi +  Zi)
  z <-  Zi / (Xi + Yi +  Zi)
  
  # calculate hue angle
  alpha <- atan2( (x - 0.33), (y - 0.33)) * 180/pi
  
  # make look up table for hue angle to wavelength conversion
  cie <- cccie31 %>%
    mutate(a = atan2( (x - 0.33), (y - 0.33)) * 180/pi) %>%
    dplyr::filter(wlnm <= 700) %>%
    dplyr::filter(wlnm >=380)
  
  # find nearest dominant wavelength to hue angle
  wl <- cie[as.vector(sapply(alpha,function(x) which.min(abs(x - cie$a)))), 'wlnm']
  #out <- cbind(as.data.frame(alpha), as.data.frame(wl))
  
  return(wl)
}

# assigning bands to R,G,B
R <-df$B4
G <-df$B3
B <-df$B2

# 3. defining zones using bounding boxes, these include entire system
# separating lake data between river arm and main body  
s2_df<-read_csv("s2flm_6lakes_datause.csv")

#dm_sf <- st_as_sf(bnras, coords = c("lon", "lat"), crs = 4326)
#mapview(dm_sf)

s2_w<-s2_df %>% filter(system == "waco")
s2_a<-s2_df %>% filter(system == "arrowhead")
s2_bw<-s2_df %>% filter(system == "brownwood")
s2_iv<-s2_df %>% filter(system == "ivie")
s2_bn<-s2_df %>% filter(system == "bonham")
s2_rb<-s2_df %>% filter(system == "redbluff")

## defining bboxes for waco ##

wmb <- s2_w[s2_w$lat >= 31.530244 & s2_w$lat <= 31.6095146 & s2_w$lon >= -97.2520428 & s2_w$lon <= -97.1809749,]
wna<-s2_w[s2_w$lat >= 31.5821285 & s2_w$lat <= 31.6157569 & s2_w$lon >= -97.308004 & s2_w$lon <= -97.2523861,]
wsa<-s2_w[s2_w$lat >= 31.4886804 & s2_w$lat <= 31.5311218 & s2_w$lon >= -97.265776 & s2_w$lon <= -97.2266369,]

wmb$zone="body"
wna$zone="arm"
wsa$zone="arm"

waco_bind<-rbind(wmb,wna,wsa)
#write_csv(waco_bind, "waco_zone.csv")

## defining bboxes for ah ##
amb<-s2_a[s2_a$lat >= 33.6933513 & s2_a$lat <= 33.7744378 & s2_a$lon >= -98.41832729 & s2_a$lon <= -98.3015976,]
ara<-s2_a[s2_a$lat >= 33.6339153 & s2_a$lat <= 33.6936369 & s2_a$lon >= -98.4667358 & s2_a$lon <= -98.351036,]

amb$zone="body"
ara$zone="arm"

ah_bind<-rbind(amb,ara)

## defining bboxes for bw ##
bwmb<-s2_bw[s2_bw$lat >= 31.7955865 & s2_bw$lat <= 31.865884 & s2_bw$lon >= -99.0787206 & s2_bw$lon <= -98.9753804,]
bwna<-s2_bw[s2_bw$lat >= 31.8661754 & s2_bw$lat <= 31.928552 & s2_bw$lon >= -99.0583197 & s2_bw$lon <= -99.0109412,]
bwsa<-s2_bw[s2_bw$lat >= 31.8200945 & s2_bw$lat <= 31.857136 & s2_bw$lon >= -99.1343388 & s2_bw$lon <= -99.0787206,]

bwmb$zone="body"
bwna$zone="arm"
bwsa$zone="arm"

bw_bind<-rbind(bwmb,bwna,bwsa)

## defining bboxes for iv ##
imb<-s2_iv[s2_iv$lat >= 31.479319 & s2_iv$lat <= 31.5557074 & s2_iv$lon >= -99.7077124 & s2_iv$lon <= -99.574503,]
ira<-s2_iv[s2_iv$lat >= 31.5559999 & s2_iv$lat <= 31.6256021 & s2_iv$lon >= -99.8392049 & s2_iv$lon <= -99.659991,]

imb$zone="body"
ira$zone="arm"

iv_bind<-rbind(imb,ira)

## defining bboxes for bon ##
bnmb1<-s2_bn[s2_bn$lat >= 33.6512746 & s2_bn$lat <= 33.6556328 & s2_bn$lon >= -96.1523812 & s2_bn$lon <= -96.13573,]
bnmb2<-s2_bn[s2_bn$lat >= 33.6442593 & s2_bn$lat <= 33.6515473 & s2_bn$lon >= -96.1636449 & s2_bn$lon <= -96.1349774,]
bnmb1$zone="body"
bnmb2$zone="body"


bnrane<-s2_bn[s2_bn$lat >= 33.6556197 & s2_bn$lat <= 33.67240739 & s2_bn$lon >= -96.15240104 & s2_bn$lon <= -96.13789,]
bnranw<-s2_bn[s2_bn$lat >= 33.6516187 & s2_bn$lat <= 33.66576413 & s2_bn$lon >= -96.1694813 & s2_bn$lon <= -96.152401,]
bnraw<-s2_bn[s2_bn$lat >= 33.63775681 & s2_bn$lat <= 33.65047555 & s2_bn$lon >= -96.17806441 & s2_bn$lon <= -96.163731,] 
bnras<-s2_bn[s2_bn$lat >= 33.6397708 & s2_bn$lat <= 33.644201 & s2_bn$lon >= -96.1537165 & s2_bn$lon <= -96.13637867,]
bnrane$zone="arm"
bnranw$zone="arm"
bnraw$zone="arm"
bnras$zone="arm"

bn_bind<-rbind(bnmb1,bnmb2,bnrane,bnranw,bnraw,bnras)

## defining bboxes for rb ##
rmb<-s2_rb[s2_rb$lat >= 31.895208 & s2_rb$lat <= 31.9744585 & s2_rb$lon >= -103.9684497 & s2_rb$lon <= -103.8946353,]
rra<-s2_rb[s2_rb$lat >= 31.974749 & s2_rb$lat <= 32.0443283 & s2_rb$lon >= -104.01239497 & s2_rb$lon <= -103.9324008,]

rmb$zone="body"
rra$zone="arm"

rb_bind<-rbind(rmb,rra)


all_bind <- rbind(bn_bind,waco_bind,ah_bind,bw_bind,iv_bind,rb_bind)