rm(list = ls())

library(data.table)
library(dplyr)
library(rstac)
library(terra)
library(mapview)
library(httr)
library(Metrics)
library(geodrawr)
library(svDialogs)
library(rstac)
library(randomForest)
library(rasterVis)
library(RColorBrewer)
library(terrainr)
library(sf)
library(raster)
library(plotly)
library(tidyverse)
library(ncdf4)


# attempting to query FLAMe points from ACOLITE corrected L2W NETcdf
# Steps include 1) read in ncd file 2) extract bands of interest into rasters 3) stack rasters into single, multiband raster
# 4) read in FLAMe path as a df 5) convert query coords into a sf 6) ensure same crs 7) extract values from points and combine

# Step 1. read in NETcdf - trying for a single image (Brownwood 8/7) to start
bw_nc_fp<-"C:/Users/Malcolm_Macleod1/Documents/acolite_py_win/OUT/brownwood/L2W/806/S2B_MSI_2022_08_07_17_25_05_T14SMA_L2W.nc"
bw_nc<-nc_open(bw_nc_fp)

# Step 2a. extracting bands of interest - starting with B2, B3, B4 for rhos - later going to L2W products
# rhos = surface reflectance data (BOA)
rhos_blue<-ncvar_get(bw_nc, "rhos_492")
rhos_green<-ncvar_get(bw_nc, "rhos_559")
rhos_red<-ncvar_get(bw_nc, "rhos_665")

# close NetCDF file
nc_close(bw_nc)

# creating raster object for each band
rhos_blue_ras<-raster(bw_nc_fp, varname="rhos_492")
rhos_green_ras<-raster(bw_nc_fp, varname="rhos_559")
rhos_red_ras<-raster(bw_nc_fp, varname="rhos_665")


# Step 3a. stacking rhos rasters
rhos_stack<-stack(rhos_blue_ras,rhos_green_ras,rhos_red_ras)

# Step 4a. Read in FLAMe path
#s2flm_pts<-read_csv("s2flm_6lakes_datause.csv")
s2flm_pts<-read_csv("six4m_zone.csv")
s2flm_bw<-s2flm_pts %>% dplyr::filter(system=="brownwood")

# converting df to spatial object and setting same crs as stack
coordinates(s2flm_bw) <- ~lon + lat
# Set the CRS for the spatial points to WGS84
crs_points <- CRS("+proj=longlat +datum=WGS84")
projection(s2flm_bw) <- crs_points

# Reproject the spatial points to match the CRS of the raster stack
s2flm_bw_utm <- spTransform(s2flm_bw, CRS(projection(rhos_stack)))

# Print the extent of the raster stack and the reprojected points
print(extent(rhos_stack))
print(extent(s2flm_bw_utm))

# Extract values from raster stack at specified points
bw_extract <- raster::extract(rhos_stack, s2flm_bw)

# Combine the query points with the extracted values
bw_query_results <- data.frame(
  latitude = coordinates(s2flm_bw_utm)[, 2],
  longitude = coordinates(s2flm_bw_utm)[, 1],
  rhos_blue = bw_extract[, 1],
  rhos_green = bw_extract[, 2],
  rhos_red = bw_extract[, 3],
  system=s2flm_bw_utm[,2],
  turb = s2flm_bw_utm[,10]
) %>% dplyr::select(latitude,longitude,system=system.system,rhos_blue,rhos_green,rhos_red,turb=turb.turb)

bw_query_results<-bw_query_results %>% mutate(ndti = (rhos_red-rhos_green)/(rhos_red+rhos_green))

ggplot(bw_query_results,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

# doing this for the other NetCDF images
# starting with Bonham
bn_nc_fp<-"C:/Users/Malcolm_Macleod1/Documents/acolite_py_win/OUT/bonham/L2W/715/S2B_MSI_2022_07_15_17_14_29_T15STT_L2W.nc"
bn_nc<-nc_open(bn_nc_fp)

# Step 2a. extracting bands of interest - starting with B2, B3, B4 for rhos - later going to L2W products
# rhos = surface reflectance data (BOA)
bn_rhos_blue<-ncvar_get(bn_nc, "rhos_492")
bn_rhos_green<-ncvar_get(bn_nc, "rhos_559")
bn_rhos_red<-ncvar_get(bn_nc, "rhos_665")

# close NetCDF file
nc_close(bn_nc)

# creating raster object for each band
bn_rhos_blue_ras<-raster(bn_nc_fp, varname="rhos_492")
bn_rhos_green_ras<-raster(bn_nc_fp, varname="rhos_559")
bn_rhos_red_ras<-raster(bn_nc_fp, varname="rhos_665")


# Step 3a. stacking rhos rasters
bn_rhos_stack<-stack(bn_rhos_blue_ras,bn_rhos_green_ras,bn_rhos_red_ras)

# Step 4a. Read in FLAMe path
s2flm_bon<-s2flm_pts %>% dplyr::filter(system=="bonham")

# converting df to spatial object and setting same crs as stack
coordinates(s2flm_bon) <- ~lon + lat
# Set the CRS for the spatial points to WGS84
bn_crs_points <- CRS("+proj=longlat +datum=WGS84")
projection(s2flm_bon) <- bn_crs_points

# Reproject the spatial points to match the CRS of the raster stack
s2flm_bn_utm <- spTransform(s2flm_bon, CRS(projection(bn_rhos_stack)))

# Print the extent of the raster stack and the reprojected points
print(extent(bn_rhos_stack))
print(extent(s2flm_bn_utm))

# Extract values from raster stack at specified points
bon_extract <- raster::extract(bn_rhos_stack, s2flm_bon)

# Combine the query points with the extracted values
bn_query_results <- data.frame(
  latitude = coordinates(s2flm_bn_utm)[, 2],
  longitude = coordinates(s2flm_bn_utm)[, 1],
  rhos_blue = bon_extract[, 1],
  rhos_green = bon_extract[, 2],
  rhos_red = bon_extract[, 3],
  system=s2flm_bn_utm[,2],
  turb = s2flm_bn_utm[,10]
) %>% dplyr::select(latitude,longitude,system=system.system,rhos_blue,rhos_green,rhos_red,turb=turb.turb)

bn_query_results<-bn_query_results %>% mutate(ndti = (rhos_red-rhos_green)/(rhos_red+rhos_green))

ggplot(bn_query_results,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

## now with Waco
# waco 7/23
w723_nc_fp<-"C:/Users/Malcolm_Macleod1/Documents/acolite_py_win/OUT/waco/L2W/723/S2A_MSI_2022_07_23_17_25_19_T14RPV_L2W.nc"
w723_nc<-nc_open(w723_nc_fp)

# Step 2a. extracting bands of interest - starting with B2, B3, B4 for rhos - later going to L2W products
# rhos = surface reflectance data (BOA)
w723_rhos_blue<-ncvar_get(w723_nc, "rhos_492")
w723_rhos_green<-ncvar_get(w723_nc, "rhos_560")
w723_rhos_red<-ncvar_get(w723_nc, "rhos_665")

# close NetCDF file
nc_close(w723_nc)

# creating raster object for each band
w723_rhos_blue_ras<-raster(w723_nc_fp, varname="rhos_492")
w723_rhos_green_ras<-raster(w723_nc_fp, varname="rhos_560")
w723_rhos_red_ras<-raster(w723_nc_fp, varname="rhos_665")


# Step 3a. stacking rhos rasters
w723_rhos_stack<-stack(w723_rhos_blue_ras,w723_rhos_green_ras,w723_rhos_red_ras)

# Step 4a. Read in FLAMe path
s2flm_w<-s2flm_pts %>% dplyr::filter(system=="waco")

# converting df to spatial object and setting same crs as stack
coordinates(s2flm_w) <- ~lon + lat
# Set the CRS for the spatial points to WGS84
w723_crs_points <- CRS("+proj=longlat +datum=WGS84")
projection(s2flm_w) <- w723_crs_points

# Reproject the spatial points to match the CRS of the raster stack
s2flm_w723_utm <- spTransform(s2flm_w, CRS(projection(w723_rhos_stack)))

# Print the extent of the raster stack and the reprojected points
print(extent(w723_rhos_stack))
print(extent(s2flm_w723_utm))

# Extract values from raster stack at specified points
w723_extract <- raster::extract(w723_rhos_stack, s2flm_w)

# Combine the query points with the extracted values
w723_query_results <- data.frame(
  latitude = coordinates(s2flm_w723_utm)[, 2],
  longitude = coordinates(s2flm_w723_utm)[, 1],
  rhos_blue = w723_extract[, 1],
  rhos_green = w723_extract[, 2],
  rhos_red = w723_extract[, 3],
  system=s2flm_w723_utm[,2],
  turb = s2flm_w723_utm[,10]
) %>% dplyr::select(latitude,longitude,system=system.system,rhos_blue,rhos_green,rhos_red,turb=turb.turb)

w723_query_results<-w723_query_results %>% mutate(ndti = (rhos_red-rhos_green)/(rhos_red+rhos_green))

ggplot(w723_query_results,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()


## now with oh ivie
iv_nc_fp<-"C:/Users/Malcolm_Macleod1/Documents/acolite_py_win/OUT/ohivie/L2W/807/S2B_MSI_2022_08_07_17_25_20_T14RMV_L2W.nc"
iv_nc<-nc_open(iv_nc_fp)

# Step 2a. extracting bands of interest - starting with B2, B3, B4 for rhos - later going to L2W products
# rhos = surface reflectance data (BOA)
iv_rhos_blue<-ncvar_get(iv_nc, "rhos_492")
iv_rhos_green<-ncvar_get(iv_nc, "rhos_559")
iv_rhos_red<-ncvar_get(iv_nc, "rhos_665")

# close NetCDF file
nc_close(iv_nc)

# creating raster object for each band
iv_rhos_blue_ras<-raster(iv_nc_fp, varname="rhos_492")
iv_rhos_green_ras<-raster(iv_nc_fp, varname="rhos_559")
iv_rhos_red_ras<-raster(iv_nc_fp, varname="rhos_665")


# Step 3a. stacking rhos rasters
iv_rhos_stack<-stack(iv_rhos_blue_ras,iv_rhos_green_ras,iv_rhos_red_ras)

# Step 4a. Read in FLAMe path
s2flm_iv<-s2flm_pts %>% dplyr::filter(system=="ivie")

# converting df to spatial object and setting same crs as stack
coordinates(s2flm_iv) <- ~lon + lat
# Set the CRS for the spatial points to WGS84
iv_crs_points <- CRS("+proj=longlat +datum=WGS84")
projection(s2flm_iv) <- iv_crs_points

# Reproject the spatial points to match the CRS of the raster stack
s2flm_iv_utm <- spTransform(s2flm_iv, CRS(projection(iv_rhos_stack)))

# Print the extent of the raster stack and the reprojected points
print(extent(iv_rhos_stack))
print(extent(s2flm_iv_utm))

# Extract values from raster stack at specified points
iv_extract <- raster::extract(iv_rhos_stack, s2flm_iv)

# Combine the query points with the extracted values
iv_query_results <- data.frame(
  latitude = coordinates(s2flm_iv_utm)[, 2],
  longitude = coordinates(s2flm_iv_utm)[, 1],
  rhos_blue = iv_extract[, 1],
  rhos_green = iv_extract[, 2],
  rhos_red = iv_extract[, 3],
  system=s2flm_iv_utm[,2],
  turb = s2flm_iv_utm[,10]
) %>% dplyr::select(latitude,longitude,system=system.system,rhos_blue,rhos_green,rhos_red,turb=turb.turb)

iv_query_results<-iv_query_results %>% mutate(ndti = (rhos_red-rhos_green)/(rhos_red+rhos_green)) %>% filter(ndti<10 & ndti>-10)
# even with these filters the R2 is incredibly low - a product of negative red band SR values is my guess
ggplot(iv_query_results,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

## now with red bluff
rb_nc_fp<-"C:/Users/Malcolm_Macleod1/Documents/acolite_py_win/OUT/redbluff/L2W/808/S2A_MSI_2022_08_08_17_45_21_T13SER_L2W.nc"
rb_nc<-nc_open(rb_nc_fp)

# Step 2a. extracting bands of interest - starting with B2, B3, B4 for rhos - later going to L2W products
# rhos = surface reflectance data (BOA)
rb_rhos_blue<-ncvar_get(rb_nc, "rhos_492")
rb_rhos_green<-ncvar_get(rb_nc, "rhos_560")
rb_rhos_red<-ncvar_get(rb_nc, "rhos_665")

# close NetCDF file
nc_close(rb_nc)

# creating raster object for each band
rb_rhos_blue_ras<-raster(rb_nc_fp, varname="rhos_492")
rb_rhos_green_ras<-raster(rb_nc_fp, varname="rhos_560")
rb_rhos_red_ras<-raster(rb_nc_fp, varname="rhos_665")


# Step 3a. stacking rhos rasters
rb_rhos_stack<-stack(rb_rhos_blue_ras,rb_rhos_green_ras,rb_rhos_red_ras)

# Step 4a. Read in FLAMe path
s2flm_rb<-s2flm_pts %>% dplyr::filter(system=="redbluff")

# converting df to spatial object and setting same crs as stack
coordinates(s2flm_rb) <- ~lon + lat
# Set the CRS for the spatial points to WGS84
rb_crs_points <- CRS("+proj=longlat +datum=WGS84")
projection(s2flm_rb) <- rb_crs_points

# Reproject the spatial points to match the CRS of the raster stack
s2flm_rb_utm <- spTransform(s2flm_rb, CRS(projection(rb_rhos_stack)))

# Print the extent of the raster stack and the reprojected points
print(extent(rb_rhos_stack))
print(extent(s2flm_rb_utm))

# Extract values from raster stack at specified points
rb_extract <- raster::extract(rb_rhos_stack, s2flm_rb)

# Combine the query points with the extracted values
rb_query_results <- data.frame(
  latitude = coordinates(s2flm_rb_utm)[, 2],
  longitude = coordinates(s2flm_rb_utm)[, 1],
  rhos_blue = rb_extract[, 1],
  rhos_green = rb_extract[, 2],
  rhos_red = rb_extract[, 3],
  system=s2flm_rb_utm[,2],
  turb = s2flm_rb_utm[,10]
) %>% dplyr::select(latitude,longitude,system=system.system,rhos_blue,rhos_green,rhos_red,turb=turb.turb)

rb_query_results<-rb_query_results %>% mutate(ndti = (rhos_red-rhos_green)/(rhos_red+rhos_green)) 

ggplot(rb_query_results,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

## lastly with arrowhead
ah812_nc_fp<-"C:/Users/Malcolm_Macleod1/Documents/acolite_py_win/OUT/arrowhead/L2W/812/S2A_MSI_2022_08_12_17_24_41_T14SNC_L2W.nc"
ah812_nc<-nc_open(ah812_nc_fp)

# Step 2a. extracting bands of interest - starting with B2, B3, B4 for rhos - later going to L2W products
# rhos = surface reflectance data (BOA)
ah812_rhos_blue<-ncvar_get(ah812_nc, "rhos_492")
ah812_rhos_green<-ncvar_get(ah812_nc, "rhos_560")
ah812_rhos_red<-ncvar_get(ah812_nc, "rhos_665")

# close NetCDF file
nc_close(ah812_nc)

# creating raster object for each band
ah812_rhos_blue_ras<-raster(ah812_nc_fp, varname="rhos_492")
ah812_rhos_green_ras<-raster(ah812_nc_fp, varname="rhos_560")
ah812_rhos_red_ras<-raster(ah812_nc_fp, varname="rhos_665")


# Step 3a. stacking rhos rasters
ah812_rhos_stack<-stack(ah812_rhos_blue_ras,ah812_rhos_green_ras,ah812_rhos_red_ras)

# Step 4a. Read in FLAMe path
s2flm_ah<-s2flm_pts %>% dplyr::filter(system=="arrowhead")

# converting df to spatial object and setting same crs as stack
coordinates(s2flm_ah) <- ~lon + lat
# Set the CRS for the spatial points to WGS84
ah812_crs_points <- CRS("+proj=longlat +datum=WGS84")
projection(s2flm_ah) <- ah812_crs_points

# Reproject the spatial points to match the CRS of the raster stack
s2flm_ah812_utm <- spTransform(s2flm_ah, CRS(projection(ah812_rhos_stack)))

# Print the extent of the raster stack and the reprojected points
print(extent(ah812_rhos_stack))
print(extent(s2flm_ah812_utm))

# Extract values from raster stack at specified points
ah812_extract <- raster::extract(ah812_rhos_stack, s2flm_ah)

# Combine the query points with the extracted values
ah812_query_results <- data.frame(
  latitude = coordinates(s2flm_ah812_utm)[, 2],
  longitude = coordinates(s2flm_ah812_utm)[, 1],
  rhos_blue = ah812_extract[, 1],
  rhos_green = ah812_extract[, 2],
  rhos_red = ah812_extract[, 3],
  system=s2flm_ah812_utm[,2],
  turb = s2flm_ah812_utm[,10]
) %>% dplyr::select(latitude,longitude,system=system.system,rhos_blue,rhos_green,rhos_red,turb=turb.turb)

ah812_query_results<-ah812_query_results %>% mutate(ndti = (rhos_red-rhos_green)/(rhos_red+rhos_green)) 

ggplot(ah812_query_results,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

# merging all the query results into a single df
acolite_all6lakes<-rbind(iv_query_results,bn_query_results,w723_query_results,bw_query_results,rb_query_results,ah812_query_results) # 0.34 or 0.49 w/o iv
acolite_all4lakes<-rbind(w723_query_results,bw_query_results,rb_query_results,ah812_query_results) # 0.52
#iv_query_results
# bn_query_results,

# viz whole enchilada ndti_turb
ggplot(acolite_all6lakes,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()


#################################################################################################################
# DOING SAME STEPS AS ABOVE BUT WITH ONLY WATER PIXELS
# Steps include 1) read in ncd file 2) extract bands of interest into rasters 3) stack rasters into single, multiband raster
# 4) read in FLAMe path as a df 5) convert query coords into a sf 6) ensure same crs 7) extract values from points and combine

# Step 1b. read in NETcdf - trying for a single image (Brownwood 8/7) to start
bw_nc_fp<-"C:/Users/Malcolm_Macleod1/Documents/acolite_py_win/OUT/brownwood/L2W/806/S2B_MSI_2022_08_07_17_25_05_T14SMA_L2W.nc"
bw_nc<-nc_open(bw_nc_fp)

# Step 2b. extracting bands of interest - starting with B2, B3, B4 for rhos - later going to L2W products
# Rrs = remote sensing reflectance for water pixels
rrs_blue<-ncvar_get(bw_nc, "Rrs_492")
rrs_green<-ncvar_get(bw_nc, "Rrs_559")
rrs_red<-ncvar_get(bw_nc, "Rrs_665")

# close NetCDF file
nc_close(bw_nc)

# creating raster object for each band
rrs_blue_ras<-raster(bw_nc_fp, varname="Rrs_492")
rrs_green_ras<-raster(bw_nc_fp, varname="Rrs_559")
rrs_red_ras<-raster(bw_nc_fp, varname="Rrs_665")


# Step 3b. stacking rhos rasters
rrs_stack<-stack(rrs_blue_ras,rrs_green_ras,rrs_red_ras)

# Step 4b. Read in FLAMe path
#s2flm_pts<-read_csv("s2flm_6lakes_datause.csv")
s2flm_pts<-read_csv("six4m_zone.csv")
s2flm_bw<-s2flm_pts %>% dplyr::filter(system=="brownwood")

# converting df to spatial object and setting same crs as stack
coordinates(s2flm_bw) <- ~lon + lat
# Set the CRS for the spatial points to WGS84
crs_points <- CRS("+proj=longlat +datum=WGS84")
projection(s2flm_bw) <- crs_points

# Reproject the spatial points to match the CRS of the raster stack
s2flm_bw_utm <- spTransform(s2flm_bw, CRS(projection(rrs_stack)))

# Print the extent of the raster stack and the reprojected points
print(extent(rrs_stack))
print(extent(s2flm_bw_utm))

# Extract values from raster stack at specified points
bw_extract <- raster::extract(rrs_stack, s2flm_bw)

# Combine the query points with the extracted values
bw_query_results <- data.frame(
  latitude = coordinates(s2flm_bw_utm)[, 2],
  longitude = coordinates(s2flm_bw_utm)[, 1],
  rrs_blue = bw_extract[, 1],
  rrs_green = bw_extract[, 2],
  rrs_red = bw_extract[, 3],
  system=s2flm_bw_utm[,2],
  turb = s2flm_bw_utm[,10]
) %>% dplyr::select(latitude,longitude,system=system.system,rrs_blue,rrs_green,rrs_red,turb=turb.turb)

bw_query_results<-bw_query_results %>% mutate(ndti = (rrs_red-rrs_green)/(rrs_red+rrs_green))

ggplot(bw_query_results,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

# doing this for the other NetCDF images
# starting with Bonham
bn_nc_fp<-"C:/Users/Malcolm_Macleod1/Documents/acolite_py_win/OUT/bonham/L2W/715/S2B_MSI_2022_07_15_17_14_29_T15STT_L2W.nc"
bn_nc<-nc_open(bn_nc_fp)

# Step 2a. extracting bands of interest - starting with B2, B3, B4 for rrs - later going to L2W products
# rhos = surface reflectance data (BOA)
bn_rrs_blue<-ncvar_get(bn_nc, "Rrs_492")
bn_rrs_green<-ncvar_get(bn_nc, "Rrs_559")
bn_rrs_red<-ncvar_get(bn_nc, "Rrs_665")

# close NetCDF file
nc_close(bn_nc)

# creating raster object for each band
bn_rrs_blue_ras<-raster(bn_nc_fp, varname="Rrs_492")
bn_rrs_green_ras<-raster(bn_nc_fp, varname="Rrs_559")
bn_rrs_red_ras<-raster(bn_nc_fp, varname="Rrs_665")


# Step 3a. stacking rrs rasters
bn_rrs_stack<-stack(bn_rrs_blue_ras,bn_rrs_green_ras,bn_rrs_red_ras)

# Step 4a. Read in FLAMe path
s2flm_bon<-s2flm_pts %>% dplyr::filter(system=="bonham")

# converting df to spatial object and setting same crs as stack
coordinates(s2flm_bon) <- ~lon + lat
# Set the CRS for the spatial points to WGS84
bn_crs_points <- CRS("+proj=longlat +datum=WGS84")
projection(s2flm_bon) <- bn_crs_points

# Reproject the spatial points to match the CRS of the raster stack
s2flm_bn_utm <- spTransform(s2flm_bon, CRS(projection(bn_rrs_stack)))

# Print the extent of the raster stack and the reprojected points
print(extent(bn_rrs_stack))
print(extent(s2flm_bn_utm))

# Extract values from raster stack at specified points
bon_extract <- raster::extract(bn_rrs_stack, s2flm_bon)

# Combine the query points with the extracted values
bn_query_results <- data.frame(
  latitude = coordinates(s2flm_bn_utm)[, 2],
  longitude = coordinates(s2flm_bn_utm)[, 1],
  rrs_blue = bon_extract[, 1],
  rrs_green = bon_extract[, 2],
  rrs_red = bon_extract[, 3],
  system=s2flm_bn_utm[,2],
  turb = s2flm_bn_utm[,10]
) %>% dplyr::select(latitude,longitude,system=system.system,rrs_blue,rrs_green,rrs_red,turb=turb.turb)

bn_query_results<-bn_query_results %>% mutate(ndti = (rrs_red-rrs_green)/(rrs_red+rrs_green))

ggplot(bn_query_results,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

## now with Waco
# waco 7/23
w723_nc_fp<-"C:/Users/Malcolm_Macleod1/Documents/acolite_py_win/OUT/waco/L2W/723/S2A_MSI_2022_07_23_17_25_19_T14RPV_L2W.nc"
w723_nc<-nc_open(w723_nc_fp)

# Step 2a. extracting bands of interest - starting with B2, B3, B4 for rrs - later going to L2W products
# rrs = surface reflectance data (BOA)
w723_rrs_blue<-ncvar_get(w723_nc, "Rrs_492")
w723_rrs_green<-ncvar_get(w723_nc, "Rrs_560")
w723_rrs_red<-ncvar_get(w723_nc, "Rrs_665")

# close NetCDF file
nc_close(w723_nc)

# creating raster object for each band
w723_rrs_blue_ras<-raster(w723_nc_fp, varname="Rrs_492")
w723_rrs_green_ras<-raster(w723_nc_fp, varname="Rrs_560")
w723_rrs_red_ras<-raster(w723_nc_fp, varname="Rrs_665")


# Step 3a. stacking rrs rasters
w723_rrs_stack<-stack(w723_rrs_blue_ras,w723_rrs_green_ras,w723_rrs_red_ras)

# Step 4a. Read in FLAMe path
s2flm_w<-s2flm_pts %>% dplyr::filter(system=="waco")

# converting df to spatial object and setting same crs as stack
coordinates(s2flm_w) <- ~lon + lat
# Set the CRS for the spatial points to WGS84
w723_crs_points <- CRS("+proj=longlat +datum=WGS84")
projection(s2flm_w) <- w723_crs_points

# Reproject the spatial points to match the CRS of the raster stack
s2flm_w723_utm <- spTransform(s2flm_w, CRS(projection(w723_rrs_stack)))

# Print the extent of the raster stack and the reprojected points
print(extent(w723_rrs_stack))
print(extent(s2flm_w723_utm))

# Extract values from raster stack at specified points
w723_extract <- raster::extract(w723_rrs_stack, s2flm_w)

# Combine the query points with the extracted values
w723_query_results <- data.frame(
  latitude = coordinates(s2flm_w723_utm)[, 2],
  longitude = coordinates(s2flm_w723_utm)[, 1],
  rrs_blue = w723_extract[, 1],
  rrs_green = w723_extract[, 2],
  rrs_red = w723_extract[, 3],
  system=s2flm_w723_utm[,2],
  turb = s2flm_w723_utm[,10]
) %>% dplyr::select(latitude,longitude,system=system.system,rrs_blue,rrs_green,rrs_red,turb=turb.turb)

w723_query_results<-w723_query_results %>% mutate(ndti = (rrs_red-rrs_green)/(rrs_red+rrs_green))

ggplot(w723_query_results,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()


## now with oh ivie
iv_nc_fp<-"C:/Users/Malcolm_Macleod1/Documents/acolite_py_win/OUT/ohivie/L2W/807/S2B_MSI_2022_08_07_17_25_20_T14RMV_L2W.nc"
iv_nc<-nc_open(iv_nc_fp)

# Step 2a. extracting bands of interest - starting with B2, B3, B4 for rrs - later going to L2W products
# rrs = surface reflectance data (BOA)
iv_rrs_blue<-ncvar_get(iv_nc, "Rrs_492")
iv_rrs_green<-ncvar_get(iv_nc, "Rrs_559")
iv_rrs_red<-ncvar_get(iv_nc, "Rrs_665")

# close NetCDF file
nc_close(iv_nc)

# creating raster object for each band
iv_rrs_blue_ras<-raster(iv_nc_fp, varname="Rrs_492")
iv_rrs_green_ras<-raster(iv_nc_fp, varname="Rrs_559")
iv_rrs_red_ras<-raster(iv_nc_fp, varname="Rrs_665")


# Step 3a. stacking rrs rasters
iv_rrs_stack<-stack(iv_rrs_blue_ras,iv_rrs_green_ras,iv_rrs_red_ras)

# Step 4a. Read in FLAMe path
s2flm_iv<-s2flm_pts %>% dplyr::filter(system=="ivie")

# converting df to spatial object and setting same crs as stack
coordinates(s2flm_iv) <- ~lon + lat
# Set the CRS for the spatial points to WGS84
iv_crs_points <- CRS("+proj=longlat +datum=WGS84")
projection(s2flm_iv) <- iv_crs_points

# Reproject the spatial points to match the CRS of the raster stack
s2flm_iv_utm <- spTransform(s2flm_iv, CRS(projection(iv_rrs_stack)))

# Print the extent of the raster stack and the reprojected points
print(extent(iv_rrs_stack))
print(extent(s2flm_iv_utm))

# Extract values from raster stack at specified points
iv_extract <- raster::extract(iv_rrs_stack, s2flm_iv)

# Combine the query points with the extracted values
iv_query_results <- data.frame(
  latitude = coordinates(s2flm_iv_utm)[, 2],
  longitude = coordinates(s2flm_iv_utm)[, 1],
  rrs_blue = iv_extract[, 1],
  rrs_green = iv_extract[, 2],
  rrs_red = iv_extract[, 3],
  system=s2flm_iv_utm[,2],
  turb = s2flm_iv_utm[,10]
) %>% dplyr::select(latitude,longitude,system=system.system,rrs_blue,rrs_green,rrs_red,turb=turb.turb)

iv_query_results<-iv_query_results %>% mutate(ndti = (rrs_red-rrs_green)/(rrs_red+rrs_green)) %>% filter(ndti<10 & ndti>-10)
# even with these filters the R2 is incredibly low - a product of negative red band SR values is my guess
ggplot(iv_query_results,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

## now with red bluff
rb_nc_fp<-"C:/Users/Malcolm_Macleod1/Documents/acolite_py_win/OUT/redbluff/L2W/808/S2A_MSI_2022_08_08_17_45_21_T13SER_L2W.nc"
rb_nc<-nc_open(rb_nc_fp)

# Step 2a. extracting bands of interest - starting with B2, B3, B4 for rrs - later going to L2W products
# rrs = surface reflectance data (BOA)
rb_rrs_blue<-ncvar_get(rb_nc, "Rrs_492")
rb_rrs_green<-ncvar_get(rb_nc, "Rrs_560")
rb_rrs_red<-ncvar_get(rb_nc, "Rrs_665")

# close NetCDF file
nc_close(rb_nc)

# creating raster object for each band
rb_rrs_blue_ras<-raster(rb_nc_fp, varname="Rrs_492")
rb_rrs_green_ras<-raster(rb_nc_fp, varname="Rrs_560")
rb_rrs_red_ras<-raster(rb_nc_fp, varname="Rrs_665")


# Step 3a. stacking rrs rasters
rb_rrs_stack<-stack(rb_rrs_blue_ras,rb_rrs_green_ras,rb_rrs_red_ras)

# Step 4a. Read in FLAMe path
s2flm_rb<-s2flm_pts %>% dplyr::filter(system=="redbluff")

# converting df to spatial object and setting same crs as stack
coordinates(s2flm_rb) <- ~lon + lat
# Set the CRS for the spatial points to WGS84
rb_crs_points <- CRS("+proj=longlat +datum=WGS84")
projection(s2flm_rb) <- rb_crs_points

# Reproject the spatial points to match the CRS of the raster stack
s2flm_rb_utm <- spTransform(s2flm_rb, CRS(projection(rb_rrs_stack)))

# Print the extent of the raster stack and the reprojected points
print(extent(rb_rrs_stack))
print(extent(s2flm_rb_utm))

# Extract values from raster stack at specified points
rb_extract <- raster::extract(rb_rrs_stack, s2flm_rb)

# Combine the query points with the extracted values
rb_query_results <- data.frame(
  latitude = coordinates(s2flm_rb_utm)[, 2],
  longitude = coordinates(s2flm_rb_utm)[, 1],
  rrs_blue = rb_extract[, 1],
  rrs_green = rb_extract[, 2],
  rrs_red = rb_extract[, 3],
  system=s2flm_rb_utm[,2],
  turb = s2flm_rb_utm[,10]
) %>% dplyr::select(latitude,longitude,system=system.system,rrs_blue,rrs_green,rrs_red,turb=turb.turb)

rb_query_results<-rb_query_results %>% mutate(ndti = (rrs_red-rrs_green)/(rrs_red+rrs_green)) 

ggplot(rb_query_results,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

## lastly with arrowhead
ah812_nc_fp<-"C:/Users/Malcolm_Macleod1/Documents/acolite_py_win/OUT/arrowhead/L2W/812/S2A_MSI_2022_08_12_17_24_41_T14SNC_L2W.nc"
ah812_nc<-nc_open(ah812_nc_fp)

# Step 2a. extracting bands of interest - starting with B2, B3, B4 for rrs - later going to L2W products
# rrs = surface reflectance data (BOA)
ah812_rrs_blue<-ncvar_get(ah812_nc, "Rrs_492")
ah812_rrs_green<-ncvar_get(ah812_nc, "Rrs_560")
ah812_rrs_red<-ncvar_get(ah812_nc, "Rrs_665")

# close NetCDF file
nc_close(ah812_nc)

# creating raster object for each band
ah812_rrs_blue_ras<-raster(ah812_nc_fp, varname="Rrs_492")
ah812_rrs_green_ras<-raster(ah812_nc_fp, varname="Rrs_560")
ah812_rrs_red_ras<-raster(ah812_nc_fp, varname="Rrs_665")


# Step 3a. stacking rrs rasters
ah812_rrs_stack<-stack(ah812_rrs_blue_ras,ah812_rrs_green_ras,ah812_rrs_red_ras)

# Step 4a. Read in FLAMe path
s2flm_ah<-s2flm_pts %>% dplyr::filter(system=="arrowhead")

# converting df to spatial object and setting same crs as stack
coordinates(s2flm_ah) <- ~lon + lat
# Set the CRS for the spatial points to WGS84
ah812_crs_points <- CRS("+proj=longlat +datum=WGS84")
projection(s2flm_ah) <- ah812_crs_points

# Reproject the spatial points to match the CRS of the raster stack
s2flm_ah812_utm <- spTransform(s2flm_ah, CRS(projection(ah812_rrs_stack)))

# Print the extent of the raster stack and the reprojected points
print(extent(ah812_rrs_stack))
print(extent(s2flm_ah812_utm))

# Extract values from raster stack at specified points
ah812_extract <- raster::extract(ah812_rrs_stack, s2flm_ah)

# Combine the query points with the extracted values
ah812_query_results <- data.frame(
  latitude = coordinates(s2flm_ah812_utm)[, 2],
  longitude = coordinates(s2flm_ah812_utm)[, 1],
  rrs_blue = ah812_extract[, 1],
  rrs_green = ah812_extract[, 2],
  rrs_red = ah812_extract[, 3],
  system=s2flm_ah812_utm[,2],
  turb = s2flm_ah812_utm[,10]
) %>% dplyr::select(latitude,longitude,system=system.system,rrs_blue,rrs_green,rrs_red,turb=turb.turb)

ah812_query_results<-ah812_query_results %>% mutate(ndti = (rrs_red-rrs_green)/(rrs_red+rrs_green)) 

ggplot(ah812_query_results,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

# merging all the query results into a single df
acolite_all6lakes<-rbind(iv_query_results,bn_query_results,w723_query_results,bw_query_results,rb_query_results,ah812_query_results) # 0.34 or 0.49 w/o iv
acolite_all4lakes<-rbind(w723_query_results,bw_query_results,rb_query_results,ah812_query_results) # 0.52
#iv_query_results
# bn_query_results,

# viz whole enchilada ndti_turb
ggplot(acolite_all6lakes,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()
################################################################################################################

# L2W products - turb Nechad (2009, 2016), turb Dogliotti, turb Novoa, hue angle
turb_n09<-ncvar_get(bw_nc, "TUR_Nechad2009_665")
turb_n16<-ncvar_get(bw_nc, "TUR_Nechad2016_665")
turb_dog<-ncvar_get(bw_nc, "TUR_Dogliotti2015")
turb_nov<-ncvar_get(bw_nc, "TUR_Novoa2017")
hue_angle<-ncvar_get(bw_nc, "hue_angle")























#################################################################################################################
# reading in GeoTIFF transformed in QGIS
bw_tif <- rast("C:/Users/Malcolm_Macleod1/s2flame/geotiff_output/brownwood806convert.tif")
bw_df <-as.data.frame(bw_tif)
colnames(bw_df) <- c("b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","lat","lon","sza","vza","raa")
bw_bands <- bw_df %>% dplyr::select(b1,b2,b3,b4,b5,lat,lon)


bon_tif <- rast("C:/Users/Malcolm_Macleod1/s2flame/geotiff_output/bonham715convert.tif")
bon_df <-as.data.frame(bon_tif)
colnames(bon_df) <- c("b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","raa","lon","lat","sza","vza")
bon_bands <- bon_df %>% dplyr::select(b1,b2,b3,b4,b5,lat,lon)

# assigning bands to R,G,B
R <- bon_bands$b4
G <- bon_bands$b3
B <- bon_bands$b2

bon_bands <- bon_bands %>% mutate(dwl = fui.hue(R, G, B))

plot(bw_tif)

# reading in L2W
ah_tif<-rast("C:/Users/Malcolm_Macleod1/Documents/acolite_py_win/OUT/arrowhead/L2W/812/S2A_MSI_2022_08_12_17_24_41_T14SNC_L2W.nc")
ah_df <-as.data.frame(ah_tif)


# sinking ncdf info to txt
sink <- nc_open("C:/Users/Malcolm_Macleod1/Documents/acolite_py_win/OUT/arrowhead/L2W/812/S2A_MSI_2022_08_12_17_24_41_T14SNC_L2W.nc")

{
  sink('ah812l2wtest.txt')
  print(sink)
  sink()
}
#######################################################################################################################################