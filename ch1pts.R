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
# re-applying all filters except SCL = 6, and applied the band2 threshold > 300, obs drop to 20013
# at this point the r2 is about the same (0.59 vs 0.6)
# applying SCL = 6 has r2 of 0.62 n = 19904 , removing Waco (which should be done) n=18317 and r2=0.71
# this is with log10 turb~ndti, turb~ndti has r2 = 0.73 but bonham and ivie show even larger ndti spread
s2flm_fullfilter<-read_csv("six4m_zone_b2thresh.csv") %>% filter(SCL=="6") %>% filter(!system=="waco")


# Reading in data for distance from dam transects queried in ee.points Shiny
# inlcudes all output data plus NDTI and normalized DFD
dfd_transect<- read_csv("dfd_transect.csv")

# Reading in sample to sensor validation data (using method from AH paper) 
# refer to "station_YSI.R" for making of sausage
sample_ysi<-read_csv("sensor_sample_valmeans.csv")

# Reading in sentinel-2 data for the station lat lons
# note that this calculates dwl/ndti for hi-qual images for AH (7/28) and Waco(7/25)
s2crasr<-read_csv("s2crasr_merge.csv")


# 4. Plotting data

ggplot(s2flame_znpts,aes(log10(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

ggplot(s2flm_fullfilter,aes(log10(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()


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

zonestats_tsi<- s2flame_znpts %>% group_by(zone) %>% 
  summarise(n=n(),mean=mean(tsi_sd),sd=sd(tsi_sd),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(tsi_sd), max=max(tsi_sd),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

s2flm_river<- s2flame_znpts %>% filter(zone=="arm")
s2flm_lacus<- s2flame_znpts %>% filter(zone=="body")

tsi_r<-s2flm_river$tsi_sd
tsi_l<-s2flm_lacus$tsi_sd
quantile(tsi_r, c(.05, 0.5, .95))
quantile(tsi_l, c(.05, 0.5, .95))

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


