# a README for six_res - a repository with figures and data for MM chapter 1

### GITHUB REPOSITORY LINK<br> <https://github.com/malcolmacleod/six_res>

### FILES<br>

These are the main ones needed to run the chapter1 scripts:

*"six4m_zone_b2thresh.csv"* - Combined FLAMe and Sentinel-2 Level-2 SR data (atmospherically corrected by Sen2Cor and dowloaded from GEE) for points along boat path the Band 2 threshold determined from Lake Bonham has been applied to exclude cloud shadows. Lake Waco is still included and SCL == 6 filter has yet to be applied, this is done in script.

*"six4m_zone.csv"* - the same data without the cloud shadow mask, the one above is preferred and performs better

*"l2w6lakes.csv"* - Combined FLAMe and Sentinel-2 Level-2 water product from the ACOLITE atmospheric correction. TOA images were downloaded and ACOLITE was applied to them producing a more stringent product designed specifically for inland waters

*"dfd_transect.csv"* - Sentinel 2 Level 2 SR data for points queried along a longitudinal transect from dam to river arm for each reservoir. Unlike the boat path data, this uses the high quality image closest to sampling date for Waco and Arrowhead which had glint and clouds, respectively

*"sensor_sample_valmeans.csv"* - original validation data frame of discrete stations and YSI data aggregated for a buffer around the station. Only has turb and chl-a

*"valmeans_secchi.csv"* - validation data frame similar to the one above but includes all CRASR variables plus predicted Secchi disk depth, could be used for TSI too

*"s2crasr_merge.csv"* - Sentinel-2 Level 2 SR data for the station points with DWL and NDTI already calculated

In the "csv" folder there are some other files such as all the CRASR data for these stations and the NLA data with our stations for reference

### SCRIPTS<br>

"ch1pts.R" - This is the main script of functional code used to produce the plots and data analysis for this manuscript. It is a conglomeration of lines from more exploratory scripts which can be found in the "scripts" folder

"chapter1.rmd" This is a RMarkdown file that is currently structured to quickly test out plots and relationships but will be tidied up into an informative HTML with the essentials of this story

Within the "scripts" folder you will find "acolite6lakes.R" used for analysis of ACOLITE corrected L2W data, "dwl6lakes.R" is used for whole system calculations of dominant wavelength such as the histograms and DWL maps, "s2flm_pts.R" is the primary exploratory script for everything related to points along the boat path and distance from dam, "sdd_turb_model.R" shows how NLA reservoir data was used to predict Secchi along boat path and preliminary code for making 6panel Secchi map, "sp_secchimap_6panel.R" is SP's revision of the Secchi map plot to add scale bars, common scale, and the dam locations, "station_ysi.R" is the script used to create the validation data frames for CRASR to YSI data along with CRASR to S2

### PLOTS<br>

"ch1results" is the folder with all the plots with the "main_plots" subfolder containing all the plots that will be included in the manuscript, as figures are revised, this folder will be updated. "exploratory" folder has a lot of the additional plots that did not make it to "main_plots" and "practice_plots" has various iterations from the process
