## six_res - A repository with figures and data for MM chapter 1


# Criteria for this data was speed >0, unique lat/lon, and "no flag" (as specified by AH paper), and SCL =6
# Filters dropped obs from 48033 to 28847; see "s2flm_pts.R" for making of sausage
# removing the 1545 Wacos "obs" the r2=0.66
s2flame_znpts<-read_csv("six4m_zone.csv")%>% filter(!system=="waco")
