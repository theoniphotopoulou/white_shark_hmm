# In this script ARGOS satellite tracking data from 33 white sharks tagged 
# in South Africa are read in and processed for statistical analysis.
# Locations with error class Z have already been removed as invalid.

# These data and the associated results are presented in detail in 
# "Sex and size influence the spatiotemporal distribution of white sharks, 
# with implications for interactions with fisheries and spatial management 
# in the southwest Indian Ocean" by Alison Kock, Amanda T Lombard, Ryan Daly, 
# Victoria Goodall, Michael Meÿer, Ryan Johnson, Chris Fischer, Pieter Koen, 
# Dylan Irion, Enrico Gennari, Alison Towner, Oliver Jewell, Charlene da Silva, 
# Matt Dicken, Malcolm Smale, Theoni Photopoulou (2021). 

# The data custodian for these data is Alison Kock and 
# the code was written by Theoni Photopoulou (20211017)

#### Load required packages and set time zone ####
require(lubridate)
require(rgdal)
require(ggplot2)
require(ggmap)
require(pgirmess)
require(move)
require(Hmisc)
require(foreach)
require(tidyr)
require(moveVis)
require(trip)
require(purrr)
require(dplyr)
require(readr)
require(sf)
require(here)

here()
Sys.setenv(TZ='UTC')

# Graphics windows are generated using quartz(). If you are
# using a Windows machine please replace quartz() with dev.new()

#### Read in uncorrected Argos satellite white shark tracks #### 
# with associated metadata

# from local data folder
tracks <- read_csv(file=here("data","south_african_white_shark_ARGOS_tracks.csv"), 
                   col_type=cols(DeployID = col_character(),
                                 SPOT = col_double(),
                                 date = col_datetime(format = ""),
                                 type = col_character(),
                                 quality = col_character(),
                                 latitude = col_double(),
                                 longitude = col_double(),
                                 area_tagged = col_character(),
                                 TL = col_double(),
                                 maturity = col_character()
                                 )
                   )
head(tracks)
names(tracks)
dim(tracks)

# If loading data from published archive
# remotes::install_github("inbo/inborutils")
# inborutils:::download_zenodo(doi="10.5281/zenodo.3820359", 
#                              path = ".", 
#                              parallel = FALSE, 
#                              quiet = FALSE)

# ** This might be needed for the csv downloaded from Zenodo **
# Reformat Date (datetime)
# mydt <- dmy_hms(as.character(tracks$Date), tz="UTC", usetz=T); head(mydt)
# length(tracks$Date); length(mydt)
# tracks$Date <- mydt[-length(mydt)]
# length(unique(tracks$Date))

# Reorder the factor levels of the location quality 
# All points with Z quality have already been removed
table(tracks$quality)
class(tracks$quality)
qual <- factor(tracks$quality, levels=c("B","A","0","1","2","3")) 
head(qual)
tracks$quality <- qual
table(tracks$quality)
tracks <- droplevels(tr)
table(tracks$quality)

# Order the data by shark ID and Date and check
# for any duplicated time stamps. This should yield 
# an empty tibble
tracks <- tracks[order(tracks$DeployID, tracks$date),]
tracks %>% 
  group_by(DeployID,date) %>% 
  filter(n()>1)

#### Do simple location filtering ####
# Do some simple, sensible filtering before you use 
# a model to filter the locations

library(doParallel)
library(argosfilter)

# Split the data into a list object where each element
# belongs to a separate animal
split_data <- split(tracks, tracks$DeployID)
length(split_data)

# Apply filter using a speed filter of 3 metres per sec 
# (vmax=3) after Nakamura et al. 2011 
# This produces a character vector of location description
# The ang=-1 argument tells the function not to remove any
# "spikes" in the track
registerDoParallel(cores=2)
tracks$filtered <- foreach(i = 1:length(split_data), 
                           .combine = c) %dopar% {
  argosfilter::sdafilter(
    lat=split_data[[i]]$latitude, 
    lon=split_data[[i]]$longitude, 
    dtime=split_data[[i]]$date,
    lc=split_data[[i]]$quality, 
    ang=-1,
    vmax=3)
}
stopImplicitCluster()

# Check how many locations are flagged as erroneous ("removed")
table(tracks$filtered)

# Keep only the locations labelled "not" (ie not removed by filter)
# These data do not have any information on axis-specific error 
tracks <- tracks %>% 
  dplyr::filter(., filtered=="not") %>% 
  arrange(.,DeployID, date) %>%
  select(-filtered)

head(tracks)
dim(tracks)

# Split tracks by individual ID again, and check how many
# locations are retained for each individual
tracks <- tracks %>% 
  mutate(longitude1=longitude, 
         latitude1=latitude)
tl <- split(x=tracks, f=tracks$DeployID)
tracks %>% group_by(DeployID) %>% count() %>% print(n=Inf)

# We use a Lamberts Equal Area (LAEA) coordinate reference
# system for these data 
range(tracks$longitude); mean(tracks$longitude)
range(tracks$latitude); mean(tracks$latitude)

laea_crs <- "+proj=laea +lat_0=-33.5 +lon_0=27.8 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'"
wgs84 <- "+proj=longlat +units=km +datum=WGS84"

# Create a spatial points dataframe for each DeployID 
ff <- foreach(i=c(1:length(tl))) %do% {
  sub <- tl[[i]]
  sptr <- sf::st_as_sf(x = sub, 
                       coords = c("longitude1", "latitude1"),
                       crs=wgs84)
  sptr <- sf::st_transform(sptr, crs = st_crs(laea_crs))
}

# Check projection is LAEA
st_crs(ff[[1]])[[1]]
st_geometry(ff[[1]])
head(ff[[1]])

# This function works out Euclidean distance between two sets of points.
# It's an backup option for calculating distance between consecutive points
# within the LAEA projection as opposed to using Latitude and Longitude
my_dist <- function(x,y){sqrt((x - lag(x))^2 + (y - lag(y))^2)}

# ff is a list of spatial objects that are LAEA projected 
# Work out the length of each track segment (seglength) in km and 
# the time in days (timediff) elapsed between consecutive locations 
# using the spatial object
jj <- foreach(i=c(1:length(ff))) %do% {
  sub.df <- ff[[i]]
	sub.coords <- st_coordinates(ff[[i]])
	# seglen1 <- my_dist(x=sub.coords[,1], y=sub.coords[,2]) # this is in meters 
	
	# You can compare my distances with those computed using the my_dist()
	# function if you would prefer
	sub.ll <- sub.df %>% select(longitude, latitude) %>% as.matrix()
  seglen <- geodist::geodist(sub.ll[,1:2], 
                              measure = "geodesic", 
                              sequential = TRUE) # this is in meters 

	difftms <- difftime(time1=sub.df$date[2:nrow(sub.df)], 
	                    time2=sub.df$date[1:(nrow(sub.df)-1)], 
	                    tz="UTC", units="days")
	sub.df <- sub.coords %>%
	          cbind(as.data.frame(sub.df), 
	             timediff=c(NA,difftms), 
	             seglength=c(NA,seglen/1000)) %>% 
	  select(-geometry) 
	  
	sub.df
	}
	
glimpse(jj[[1]])

# Turn the jj list back into a dataframe with x and y coordinates
ttracks <- dplyr::bind_rows(jj) 
row.names(ttracks) <- NULL
class(ttracks)
dim(ttracks)
head(ttracks)

#### Simple plotting of the uncorrected tracks ####
ftdf_plot <- qmplot(longitude, latitude, data = ttracks,
                    geom="blank", maptype="watercolor", color="bw") +
  geom_point(data = ttracks, aes(x = longitude, y = latitude, colour = as.factor(DeployID)), size = 0.2) + 
  geom_path(data = ttracks, aes(x = longitude, y = latitude, colour = as.factor(DeployID)), size = 0.2) + 
  theme_bw(base_size=20) # + facet_wrap(~ ID)
quartz(); ftdf_plot 

# There are some big space and time gaps in the data
# Find where the difference between consecutive points 
# exceeds 3 days or 400km and create a new label to indicate unique 
# segments of tracks that should be treated separately
stracks <- ttracks
tgna <- ifelse(is.na(stracks$timediff), -6, stracks$timediff) 
dgna <- ifelse(is.na(stracks$seglength), -6, stracks$seglength)
# Only 400km or less allowed between points
distgap <- ifelse(dgna>400, 1, 0) 
# Only 3 days or less allowed between points
#   If the shark moves at 3m/sec it can cover
#   a maximum of 3*3600*24*4/1000 = 777.6 km in 3 days
timegap <- ifelse(tgna>3, 1, 0) 
cstg <- cumsum(timegap)
cstg <- ifelse(cstg < 10, paste("00", cstg, sep=""), 
               ifelse(cstg < 100, paste("0", cstg, sep=""), paste(cstg)))
# Create a trip tag that follows these time and distance rules
triptag <- paste(stracks$DeployID, cstg, sep="_")
stracks$triptag <- as.factor(triptag) 
table(stracks$triptag)

# See how many locations each trip has in it
triptally <- stracks %>%
  group_by(triptag) %>%
  dplyr::tally() 
print.data.frame(triptally)
  
# Attach the trip tag lable to the main dataset
fstracks <- left_join(stracks, triptally, by="triptag")

head(fstracks)
summary(fstracks$n)
length(unique(fstracks$DeployID))

# Simple rules to ensure model convergence below
#         ** Arrived at empirically **
# Only take segments more than 25 points long: 
fftracks <- filter(fstracks, n > 25) 

# No gaps longer than 3 days between points and 
# only segments longer than 25 locations works well 
# Shark Helen gets dropped at this stage and 
# we are left with 33 sharks
length(unique(fftracks$DeployID))
fftracks <- droplevels(fftracks)
dim(fftracks)

# We need to work out seglength and timediff (days) again 
# after removing points
# Create a spatial points dataframe 
ftdfsp <- sf::st_as_sf(x = fftracks, 
                       coords = c("X", "Y"),
                       crs=laea_crs)
# Check projection is LAEA
st_crs(ftdfsp)[[1]]

hh <- split(ftdfsp, ftdfsp$triptag)
st_crs(hh[[1]])[[1]]

# ff is a list of spatial objects that are LAEA projected 
# Work out the length of each track segment (seglength) in km and 
# the time in days (timediff) elapsed between consecutive locations 
# using the spatial object
kk <- foreach(i=c(1:length(hh))) %do% {
  sub.df <- hh[[i]]
  sub.coords <- st_coordinates(hh[[i]])
  #seglen <- my_dist(x=sub[,1], y=sub[,2]) # this is in meters 
  
  # You can compare my distances with those computed using the my_dist()
  # function if you would prefer
  sub.ll <- sub.df %>% select(longitude, latitude) %>% as.matrix()
  seglen <- geodist::geodist(sub.ll[,1:2], 
                             measure = "geodesic", 
                             sequential = TRUE) # this is in meters 
  
  difftms <- difftime(time1=sub.df$date[2:nrow(sub.df)], 
                      time2=sub.df$date[1:(nrow(sub.df)-1)], 
                      tz="UTC", units="days")
  sub.df <- sub.coords %>%
    cbind(as.data.frame(sub.df)) %>%
    mutate(timediff=c(NA,difftms), 
          seglength=c(NA,seglen/1000)) %>% 
    select(-geometry) 
  
  sub.df
}

glimpse(kk[[1]])

# Turn the kk list back into a dataframe with x and y coordinates
ftracks <- kk
glimpse(ftracks[[15]])
ftdf <- dplyr::bind_rows(ftracks) 
row.names(ftdf) <- NULL
class(ftdf)
dim(ftdf)
head(ftdf)

# TRACK STATS for all sharks except Helen 
#  (got dropped at the 25 location step)

# How many locations does each individual have?
nloc_df <- ftdf %>% 
  group_by(SPOT) %>% 
  tally() %>% print(n=Inf)

# What is the mean time in days between locations? 
ftdf %>% dplyr::select(timediff) %>% 
  summarise(min_timediff = min(., na.rm=TRUE),
            max_timediff = max(., na.rm=TRUE))
mean(ftdf$timediff, na.rm=TRUE)*24 # timediff is is days
sd(ftdf$timediff, na.rm=TRUE)*24

# What is the duration in days per track segment?
tripdurdf <- ftdf %>% group_by(triptag) %>% 
  summarise(tripdur = sum(timediff, na.rm=T)) # days
tripdurdf %>% print(n=Inf)
summary(tripdurdf$tripdur)

# What is the duration in days for each individual's cumulative track?
trackdurdf <- ftdf %>% group_by(SPOT) %>% 
  summarise(trackdur = sum(timediff, na.rm=T)) # days
trackdurdf %>% print(n=Inf)
summary(trackdurdf$trackdur)

# How many individual days were there locations on?
n_days <- ftdf %>% group_by(SPOT) %>% 
  dplyr::select(date) %>% 
  mutate(date_only = date(date)) %>%
  distinct(date_only) %>%
  tally() %>% print(n=Inf)

# Extract metadata from tracking data object
metadata <- ftdf %>% 
  group_by(DeployID) %>%
  summarise(SPOT=unique(SPOT),
            sex=unique(sex),
            TL=unique(TL),
            area_tagged=unique(area_tagged),
            maturity=unique(maturity))
head(metadata)
save(metadata, file=here("data","metadata.RData"))

# Create a table of number of days tracked and locations
locations_table <- left_join(n_days, nloc_df, by="SPOT") %>%
  rename(n_days = n.x, n_locs = n.y) %>%
  mutate(nloc_per_day = n_locs/n_days) %>% 
  left_join(., metadata, by="SPOT")
locations_table %>%
  print(n=Inf)

###
# there is no difference in: 
# the number of locations per day
# the number of days with locations
# the number of locations
# between Female and Male sharks
m_rates <- locations_table %>% filter(sex=="M") %>% 
  dplyr::select(nloc_per_day) %>% as.data.frame() 
f_rates <- locations_table %>% filter(sex=="F") %>% 
  dplyr::select(nloc_per_day) %>% as.data.frame() 
t.test(f_rates, m_rates)  
m_nlocs <- locations_table %>% filter(sex=="M") %>% 
  dplyr::select(n_locs) %>% as.data.frame() 
f_nlocs <- locations_table %>% filter(sex=="F") %>% 
  dplyr::select(n_locs) %>% as.data.frame() 
t.test(f_nlocs, m_nlocs)  
m_ndays <- locations_table %>% filter(sex=="M") %>% 
  dplyr::select(n_days) %>% as.data.frame() 
f_ndays <- locations_table %>% filter(sex=="F") %>% 
  dplyr::select(n_days) %>% as.data.frame() 
t.test(f_ndays, m_ndays)  
###


#### Model-based location filtering ####

# Correct and regularise locations using Ian Jonsen's foieGras package
require(TMB)
#devtools::install_github("ianjonsen/foieGras")
require(devtools)
require(foieGras)

# I fit a model using foieGras package version 0.2.2 and I load the 
# model object at the beginning of script 02momentuHMM_models.R 
# to use in the rest of the analysis. However, if you wanted to
# do this for your own data, this is how you would proceed:

# filter individual trips, not all trips from each individual
model_df <- dplyr::select(ftdf, "triptag", "date", "quality", 
                          "longitude", "latitude"); names(model_df)
names(model_df) <- c("id", "date", "lc", "lon", "lat")
head(model_df)
dim(model_df)

# fit a model with a 12hr time step
fgfits_12hrs <- fit_ssm(d=model_df, vmax = 3, min.dt = 60*6, pf = FALSE, 
                   time.step = 12, model = "crw", 
                   control=ssm_control(optim="optim", eval.max=1000))

fgfits_12hrs 
# "fitted" are the estimated locations at the time of observation, and 
# "predicted" are the estimated locations at the regularised time interval

#
# NOW GO TO SCRIPT: 02momentuHMM_models.R
#



