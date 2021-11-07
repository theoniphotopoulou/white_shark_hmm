# In this script ARGOS satellite tracking data from 33 white sharks tagged 
# in South Africa are statistically analysed using hidden Markov models
# implemented in the momentuHMM R package by McClintock & Michelot (2018).
# The locations used in the model have been modelled and regularised using
# the foieGras R package by Jonsen et al. (2019).

# momentuHMM modelling of shark movement behaviour
# using foieGras filtered locations

# These data and the associated results are presented in detail in 
# "sex and size influence the spatiotemporal distribution of white sharks, 
# with implications for interactions with fisheries and spatial management 
# in the southwest Indian Ocean" by Alison Kock, Amanda T Lombard, Ryan Daly, 
# Victoria Goodall, Michael Meÿer, Ryan Johnson, Chris Fischer, Pieter Koen, 
# Dylan Irion, Enrico Gennari, Alison Towner, Oliver Jewell, Charlene da Silva, 
# Matt Dicken, Malcolm Smale, Theoni Photopoulou (2021). 

# The data custodian for these data is Alison Kock and 
# the code was written by Theoni Photopoulou (20211107)

#### Load required packages and set time zone ####
require(momentuHMM)
require(circular)
require(tidyr)
require(lubridate)
require(dplyr)
require(ggplot2)
require(ggmap)
require(Hmisc)
require(foreach)
require(purrr)

require(TMB)
require(devtools)
#install_version("foieGras", version = "0.2.2", repos = "http://cran.r-project.org")
require(foieGras)
packageVersion("foieGras") 

require(rnaturalearth)
require(sf)
require(grid)
require(gridExtra)
require(cowplot)
require(stringr)
require(here)

here()
Sys.setenv(TZ='UTC')

# Graphics windows are generated using quartz(). If you are
# using a Windows machine please replace quartz() with quartz()

### Set the projection for the data
laea_crs <- "+proj=laea +lat_0=-33.5 +lon_0=27.8 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
wgs84 <- "+proj=longlat +units=km +datum=WGS84"

#### Load the model-filtered regularised tracking data
load(file=here("output","fgfits_12hrs.RData"))
ls()
 
# The predicted locations are a UTM unprojected object 
plocs <- fgpreds <- grab(fgfits_12hrs, "predicted")

# Get the x, y coordinates of the UTM unprojected object
st_crs(fgpreds)[[1]]
# It has a projection string of:
# +proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs
utmxy <- fgpreds %>% 
  sf::st_coordinates()

# Get the x, y coordinates of the LAEA projected object
laxy <- fgpreds %>% 
  sf::st_transform(laea_crs) %>% 
  sf::st_coordinates()

# Create a WGS84 lon, lat object
ll <- fgpreds %>% 
  sf::st_transform(4326) %>% 
  sf::st_coordinates()

# This strips the sf class from the tibble
sf::st_geometry(fgpreds) <- NULL

# Add lon, lat to the fgpreds dataframe
lllocs <- fgpreds %>% 
  mutate(lon = ll[,1], lat = ll[,2]) %>%
  dplyr::select(id, date, lon, lat)

# Add unprojected x,y to the fgpreds dataframe
utmlocs <- fgpreds %>% 
  mutate(x = utmxy[,1], y = utmxy[,2]) %>%
  dplyr::select(id, date, x, y, x.se, y.se)

# Add LAEA projected x,y to the fgpreds dataframe
laealocs <- fgpreds %>% 
  mutate(x = laxy[,1], y = laxy[,2]) %>%
  dplyr::select(id, date, x, y)

head(plocs) # original sf object
head(fgpreds) # original sf object without geometry column
head(utmlocs) # unprojected UTM object
head(laealocs) # LAEA projected object
head(lllocs) # lat lon projected object

# Plot the latitude-longitude data on the map
quartz()
ggplot(lllocs) +
  geom_point(aes(x = lon, y = lat, colour = as.factor(id)), size = 0.2) + 
  geom_path(aes(x = lon, y = lat, colour = as.factor(id)), size = 0.2) + 
  theme_bw(base_size=20) + 
  theme(legend.position="none")

# Plot the UTM data on the map with errors
quartz()
ggplot(utmlocs) +
  geom_point(aes(x = x, y = y, colour = as.factor(id)), size = 0.2) + 
  geom_segment(aes(x = x-(1.96*x.se), xend = x+(1.96*x.se), y = y, yend = y)) + 
  geom_segment(aes(x = x, xend = x, y = y-(1.96*y.se), yend = y+(1.96*y.se))) + 
  theme_bw(base_size=20) + 
  theme(legend.position="none")

# Use the LAEA projected data to fit an HMM to
predlocs <- laealocs
head(predlocs)

laeapreds <- plocs %>% 
  sf::st_transform(laea_crs) 
head(laeapreds)

lonlat12 <- lllocs
head(lonlat12)

# LAEA PROJECTED
# arrange by id and data and rename the id column so that momentuHMM knows what to look for
predlocs12 <- predlocs %>% 
  rename(datetime=date) %>%
  arrange(id, datetime) 
head(predlocs12) 

# Add lon and lat to the UTM object for plotting purposes only
utm_predlocs12 <- utmlocs %>% 
  rename(datetime=date) %>%
  arrange(id, datetime) %>%
  mutate(xmin = x-(1.96*x.se), xmax = x+(1.96*x.se),       # this gives the 95% CI in x
         ymin = y-(1.96*y.se), ymax = y+(1.96*y.se)) %>%   # this gives the 95% CI in y
  mutate(lon = lllocs$lon, lat = lllocs$lat, 
         laea.x = laealocs$x, laea.y = laealocs$y)

head(utm_predlocs12)

## Read in the metadata
load(file=here("data","metadata.RData"))
head(metadata)
names(metadata)

# attach the metadata to the main dataset fdftracks             
dim(metadata)
shark_id <- str_split(utm_predlocs12$id,"_", simplify=TRUE)[,1]
utm_predlocs12$DeployID <- shark_id
utm_predlocs <- left_join(utm_predlocs12, 
                          metadata, 
                          "DeployID") %>%
  rename(ID=id, 
         shark_id=DeployID)

head(utm_predlocs)
dim(utm_predlocs)

# Create a season table
season_tab <- data.frame(month = c(1:12), season = c(rep("summer",2),
                                                     rep("autumn",3),
                                                     rep("winter",3),
                                                     rep("spring",3),
                                                     "summer"))

# Create a dataframe including the projected locations
pl <- utm_predlocs %>% 
  mutate(yearday = yday(datetime), 
         month = month(datetime), 
         week = week(datetime)
        )
head(pl); 
pl <- left_join(pl, season_tab, "month")
pl$year <- year(pl$datetime)
pl$cosyd <- cos(2*pi*pl$yearday/365)
pl$sinyd <- sin(2*pi*pl$yearday/365)

head(pl)
table(pl$sex) # F=1, M=2
table(pl$month)
table(pl$week)
table(pl$season)
table(pl$year)
table(pl$area_tagged)
table(pl$maturity)
head(pl$cosyd)
head(pl$sinyd)

## Save dataset before you model it
save(pl, file=here::here("output","pl_utm.RData"))

# Prepare data for HMM
sharkdat <- prepData(data=pl, covNames = c("datetime","shark_id","area_tagged",
                                           "sex","maturity","TL",
                                           "yearday","month","week","season",
                                           "cosyd","sinyd"), 
                     coordNames = c("x", "y"), type="UTM")
head(sharkdat)
sharkdat %>% summarise(mean_step=mean(step, na.rm=TRUE), 
                       mean_angle=mean(angle, na.rm=TRUE))

sharkdat %>% group_by(ID) %>% tally() %>% print(n=Inf)
sharkdat %>% group_by(shark_id) %>% tally() %>% print(n=Inf)

dim(sharkdat)
summary(sharkdat$step)
summary(sharkdat$angle)

# Save prepData object for model fitting
save(sharkdat, file=here("output","sharkdat.RData"))

# Specify plotting parameters and make maps
alpha.trans <- 0.2
base_size <- 20
pointsize <- 1

# Create quick map - no projection
afr_countries <- ne_countries(scale="large", type="countries", continent="africa", 
                              returnclass = "sf")
safr_countries <- ne_countries(scale="large", continent="africa", 
                               country=c("south africa","mozambique","namibia","madagascar","botswana","zimbabwe"), 
                               returnclass = "sf")

# Make sure the order of the levels in the maturity variable
# is correct
class(sharkdat$maturity)
sharkdat$maturity <- factor(sharkdat$maturity, 
                             levels = c("Juvenile", "Sub_adult", "Adult"))
levels(sharkdat$maturity)

# Make sure the order of the levels in the sex variable
# is correct
class(sharkdat$sex)
sharkdat$sex <- as.factor(sharkdat$sex)
levels(sharkdat$sex)

# Look at the data from one of the sharks
# histogram of step length to get an idea of 
# what it looks like
sharkdat %>% 
  filter(shark_id=="Alisha") %>%
  ggplot() + 
  geom_histogram(aes(step/1000, stat(density))) + 
  xlab("Step (km)") 
# histogram of turning angle
sharkdat %>% filter(shark_id=="Alisha") %>%
  ggplot() + 
  geom_histogram(aes(angle, stat(density)), binwidth = 0.1) + 
  xlab("Turning angle (radians)")


#### Fit hidden Markon models #### 

# I have fitted these and load the model object
# but if you wanted to fit them to a different dataset this is wh

# Set state names and number of states
stateNames <- c("Resident", "Transient")
nbStates <- 2

# Set initial values
mu0 <- c(9000,45000)  # mean of steps
sigma0 <- c(8000,23000) # sd of steps
stepPar0 <- c(mu0,sigma0)

anglemu0 <- c(0,0) # mean of angles
anglekappa0 <- c(0.4,0.8) # concentration of angles, 
                          # higher number the more concentrated. needs to be positive. 
                          # 0.01 is a pretty much a flat line
anglePar0 <- c(anglekappa0)

# Fit a model without covariates
wsm0 <- fitHMM(sharkdat, nbStates=2,
                       dist = list(step = 'gamma', angle = 'wrpcauchy'),
                       Par0 = list(step=stepPar0, angle=anglePar0),
                       estAngleMean = list(angle = FALSE),
                       stateNames = stateNames)

# Fit a model with sex, maturity and yearday as covariates
formula <- ~ sex + maturity + cosinor(yearday, period=365)
wsm1 <- fitHMM(sharkdat, nbStates=2,
               dist = list(step = 'gamma', angle = 'wrpcauchy'),
               Par0 = list(step=stepPar0, angle=anglePar0),
               estAngleMean = list(angle = FALSE),
               formula = formula,
               stateNames = stateNames)

# Fit a model with yearday only as a covariate
formula <- ~ cosinor(yearday, period=365)
wsm2 <- fitHMM(sharkdat, nbStates=2,
               dist = list(step = 'gamma', angle = 'wrpcauchy'),
               Par0 = list(step=stepPar0, angle=anglePar0),
               estAngleMean = list(angle = FALSE),
               formula = formula,
               stateNames = stateNames)

# Fit a model with sex and yearday as covariates
formula <- ~ sex + cosinor(yearday, period=365)
wsm3 <- fitHMM(sharkdat, nbStates=2,
               dist = list(step = 'gamma', angle = 'wrpcauchy'),
               Par0 = list(step=stepPar0, angle=anglePar0),
               estAngleMean = list(angle = FALSE),
               formula = formula,
               stateNames = stateNames)

# Fit a model with maturity and yearday as covariates
formula <- ~ maturity + cosinor(yearday, period=365)
wsm4 <- fitHMM(sharkdat, nbStates=2,
               dist = list(step = 'gamma', angle = 'wrpcauchy'),
               Par0 = list(step=stepPar0, angle=anglePar0),
               estAngleMean = list(angle = FALSE),
               formula = formula,
               stateNames = stateNames)

# Fit a model with sex and maturity as covariates
formula <- ~ sex + maturity
wsm5 <- fitHMM(sharkdat, nbStates=2,
               dist = list(step = 'gamma', angle = 'wrpcauchy'),
               Par0 = list(step=stepPar0, angle=anglePar0),
               estAngleMean = list(angle = FALSE),
               formula = formula,
               stateNames = stateNames)

# Fit a model with sex, total length and yearday as covariates
formula <- ~ sex + TL + cosinor(yearday, period=365)
wsm6 <- fitHMM(sharkdat, nbStates=2, nlmPar=list(steptol=1e-8, iterlim=1000),
               dist = list(step = 'gamma', angle = 'wrpcauchy'),
               Par0 = list(step=stepPar0, angle=anglePar0),
               estAngleMean = list(angle = FALSE),
               formula = formula,
               stateNames = stateNames)

  # this is the same as the previous model except sin and cos of yearday
  # are specified manually, rather than with the cosinor() function,
  # to enable plotting
formula <- ~ sex + TL + cosyd + sinyd
wsm6a <- fitHMM(sharkdat, nbStates=2, nlmPar=list(steptol=1e-8, iterlim=1000),
               dist = list(step = 'gamma', angle = 'wrpcauchy'),
               Par0 = list(step=stepPar0, angle=anglePar0),
               estAngleMean = list(angle = FALSE),
               formula = formula,
               stateNames = stateNames)

# Fit a model with sex, total length acting on the step length mean
# and yearday as a covariate on the tpm
formula <- ~ TL + sex
DMestmean <- list(step = list(mean = ~ cosinor(yearday, period=365), sd = ~ 1),
                  angle = list(concentration = ~ 1))

# Additional initial values are required for step mean for 
# sex=M for each state.
# These appear in positions 2 and 4 of the vector below
step.mean <- c(9, 8, 7.9, 10, 7.9, 10, # mean1, mean1:cosinor, mean2, mean2:cosinor
               9, 8.8) # sd1, sd2

exp(step.mean)
angle.conc <- c(0.4, 0.8)
exp(angle.conc)

wsm7 <- fitHMM(sharkdat, nbStates=2, 
               dist = list(step = 'gamma', angle = 'wrpcauchy'),
               Par0 = list(step=step.mean, angle=angle.conc),
               DM=DMestmean,
               formula = formula,
               stateNames = stateNames)

###### LOAD FITTED MODELS (from above)
load(here::here("output","wsHMM_models.RData"))

# Look at the AIC values and AIC weights for each model
AIC(wsm0, wsm1, wsm2, wsm3, wsm4, wsm5, wsm6a, wsm7)
AICweights(wsm0, wsm1, wsm2, wsm3, wsm4, wsm5, wsm6a, wsm7)
# -> no evidence that the shape of the distributions changes with TL, sex or yearday

# Model wsm6 has the most supported: sex, total length and yearday as covariates
# I use the wsm6a version where sin and cosine are specified manually
# Models wsm6a and wsm1 are essentially the same model, except in wsm1 total length
# is categoried in to maturity levels, whereas in wsm6a total length is provided 
# as a numeric variable
print(wsm6a)
wsm6a$mle
wsm6a$mod$minimum
wsm6a$mod$code # nlm code 2: successive iterates within tolerance, current iterate is probably solution
plot(wsm6a, ask=TRUE, breaks=50)

# Model checking
ps_wsm6a <- pseudoRes(wsm6a)
qqnorm(ps_wsm6a$stepRes)
qqnorm(ps_wsm6a$angleRes)
# Both sets of residuals look nice and symmetrical
hist(ps_wsm6a$stepRes, breaks=50, xlim=c(-5,5))
hist(ps_wsm6a$angleRes, breaks=50, xlim=c(-5,5))
# There is some residual autocorrelation in step length that peters 
# out at a lag of about 60
acf(na.omit(ps_wsm6a$stepRes), lag=200) 
# This corresponds to about 30 days
60*12/24
acf(na.omit(ps_wsm6a$angleRes))

# Plot the pseudoresiduals for supplementary material
base_size <- 16

ps_wsm6a_stepdf <- data.frame(res = c(ps_wsm6a$stepRes,ps_wsm6a$angleRes),
                              var = c(rep("Step length",length(ps_wsm6a$stepRes)),
                                      rep("Turning angle",length(ps_wsm6a$angleRes)))
                              )

wsm6a_pseudores_plot <- ggplot(ps_wsm6a_stepdf) +
  geom_histogram(aes(x=res, y=..density..), colour="grey30", fill="white") +
  xlim(-5,5) +
  xlab("Pseudo-residuals") + ylab("Density") + 
  theme_bw(base_size = base_size) +
  theme(panel.grid=element_blank()) +
  facet_wrap(~ var)

wsm6a_pseudores_plot

ggsave(filename=here::here("figures","pseudoresiduals_wsm6a.jpg"), 
       plot=wsm6a_pseudores_plot,
       width=20, height=12, units="cm",dpi=700)

# Save best fitting model
save(wsm6a, file=here::here("output","wsm6a.RData"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Adult sizes
maxTL_adM <- sharkdat %>% filter(sex=="M") %>% 
  summarise(minTL=min(TL), maxTL=max(TL)) %>% as.numeric() %>% max()
maxTL_adF <- sharkdat %>% filter(sex=="F") %>% 
  summarise(minTL=min(TL), maxTL=max(TL)) %>% as.numeric() %>% max()

maxTL_adM; maxTL_adF
adultM_TL <- maxTL_adM; adultF_TL <- maxTL_adF

# Subadult sizes
meanTL_subadM <- sharkdat %>% filter(sex=="M", TL %in% c(300:360)) %>% 
  summarise(minTL=min(TL), maxTL=max(TL)) %>% as.numeric() %>% mean(.)
meanTL_subadF <- sharkdat %>% filter(sex=="F", TL %in% c(300:480)) %>% 
  summarise(minTL=min(TL), maxTL=max(TL)) %>% as.numeric() %>% mean(.)

meanTL_subadM; meanTL_subadF
subadultM_TL <- meanTL_subadM; subadultF_TL <- meanTL_subadF 

# Juvenile size
meanTL_juven <- sharkdat %>% filter(TL<300) %>% 
  summarise(minTL=min(TL), maxTL=max(TL)) %>% as.numeric() %>% mean(.)

meanTL_juven
juvenile_TL <- meanTL_juven;

# Look at the stationary distributions for different
# maturity levels at the beginning of the year
covs_adF <- data.frame(sex="F", TL=adultF_TL, yearday=1)
covs_adM <- data.frame(sex="M", TL=adultM_TL, yearday=1)

covs_subadF <- data.frame(sex="F", TL=meanTL_subadF, yearday=1)
covs_subadM <- data.frame(sex="M", TL=meanTL_subadM, yearday=1)

covs_juvF <- data.frame(sex="F", TL=meanTL_juven, yearday=1)
covs_juvM <- data.frame(sex="M", TL=meanTL_juven, yearday=1)

plotStationary(wsm6, covs=covs_adF, plotCI = TRUE)
plotStationary(wsm6, covs=covs_adM, plotCI = TRUE)

plotStationary(wsm6, covs=covs_subadF, plotCI = TRUE)
plotStationary(wsm6, covs=covs_subadM, plotCI = TRUE)

plotStationary(wsm6, covs=covs_juvF, plotCI = TRUE)
plotStationary(wsm6, covs=covs_juvM, plotCI = TRUE)

## Extract the Viterbi decoded states and state probabilities
# and attach them to the predictions dataframe
wsm6a_vit <- viterbi(wsm6a)
wsm6a_sprobs <- stateProbs(wsm6a)
wsm6a_preds <- sharkdat %>% 
  mutate(wsm6a_vit = factor(wsm6a_vit, labels=c("Resident", "Transient")),
         sex = factor(sex, labels=c("Female", "Male")),
         wsm6a_Ptrans = wsm6a_sprobs[,2])

head(wsm6a_preds)
table(wsm6a_vit)

# Check the levels of the maturity factor are in the right order
# Should be Juvenile -> Sub_adult -> Adult
levels(wsm6a_preds$maturity)

# Create a dataframe of predictions together with Viterbi decoded states
# and state probabilities
wshark_locs_states <- wsm6a_preds %>%
  dplyr::select(ID, shark_id, year, x, y, x.se, y.se, xmin, xmax, ymin, ymax, 
         lon, lat, laea.x, laea.y, datetime, TL, area_tagged, maturity, sex, 
         yearday, month, week, season, viterbi_state=wsm6a_vit, prob_transient=wsm6a_Ptrans) %>%
  mutate(lon = lllocs$lon, lat = lllocs$lat) 
head(wshark_locs_states)

# Make the predictions an sf object
wsm6a_preds_sf <- wsm6a_preds %>%
  st_as_sf(coords=c("laea.x","laea.y")) %>% 
  st_set_crs(laea_crs)

# Save for use in the plotting script 03ARGOSshark_plot_results.R
save(wshark_locs_states, wsm6a_preds_sf,
     file=here::here("output","wshark_pred_locs_states.RData"))

# Specify plotting parameters
alpha.trans <- 0.2
base_size <- 20
pointsize <- 1

# stationary state distribution - ALL
delta.avg <- as.vector(round(table(wsm6a_vit)/length(wsm6a_vit), digits=2)); delta.avg
# [1] 0.58 0.42
# MALES
wsm6a_predsM <- wsm6a_preds %>% filter(sex=="Male")
delta.avgM <- as.vector(round(table(wsm6a_predsM$wsm6a_vit)/length(wsm6a_predsM$wsm6a_vit), 
                              digits=2)); delta.avgM
# [1] 0.65 0.35
# FEMALES
wsm6a_predsF <- wsm6a_preds %>% filter(sex=="Female")
delta.avgF <- as.vector(round(table(wsm6a_predsF$wsm6a_vit)/length(wsm6a_predsF$wsm6a_vit), 
                              digits=2)); delta.avgF
# [1] 0.52 0.48

# Create plots 

###################################################################
# ~~~~~~~~~~~~ State-dependent density distributions ~~~~~~~~~~~~~~
###################################################################

require(scales); require(viridis)
state_cols <- viridis_pal(alpha = 1, begin=0.05, end=0.65, direction = 1, option = "D")(2)
marg_col <- viridis_pal(alpha = 1, begin=0.9, end=1, direction = 1, option = "D")(1)
show_col(state_cols)
show_col(marg_col)

################
### State-dependent densities of STEP LENGTH scaled by the equilibrium state densities

x <- seq(min(sharkdat$step, na.rm=TRUE), max(sharkdat$step, na.rm=TRUE), length=1000)
wsm6a$mle$step
mu1_step <- wsm6a$mle$step[1]
sd1_step <- wsm6a$mle$step[2]
mu2_step <- wsm6a$mle$step[3]
sd2_step <- wsm6a$mle$step[4]

# Work out state-dependent densities
# step
d1_step <- (dgamma(x, shape = mu1_step^2/sd1_step^2, scale = sd1_step^2/mu1_step))*delta.avg[1]
d2_step <- (dgamma(x, shape = mu2_step^2/sd2_step^2, scale = sd2_step^2/mu2_step))*delta.avg[2]

dmarg_step <- d1_step + d2_step

# Plot densities for maximum dive depth, using manual colours or a colour palette like `viridis` or `RColorBrewer`
step_dat <- sharkdat$step

quartz(); 
step_dens <- ggplot(data=data.frame(x=step_dat), aes(x,..density..)) + 
  geom_histogram(boundary=0, binwidth=5000, fill="grey90") + theme_bw() + 
  geom_line(data=data.frame(x=x, d1_step=d1_step), aes(x=x, y=d1_step, colour="Resident"), size=1.5) +
  geom_line(data=data.frame(x=x, d2_step=d2_step), aes(x, d2_step, colour="Transient"), size=1.5) +
  geom_line(data=data.frame(x=x, dmarg_step=dmarg_step), aes(x, dmarg_step, color="Marginal"), size=1.2, linetype="dashed") +
  scale_colour_manual(name="", values=c("Resident" = state_cols[1], "Transient" = state_cols[2],
                                                 "Marginal" =  marg_col), 
                      aesthetics=c("colour"), breaks=c("Resident","Transient","Marginal")) + 
  scale_x_continuous(breaks=c(0,10000,40000,100000,150000), labels=c("0","10","40","100","150"), limits=c(0,180000)) + 
  xlab("Step length (km)") + ylab("Density") + #ggtitle("State-Dependent Step Length Densities") +
  theme(legend.position=c(.7, .7)) +
  theme(text=element_text(size=20)); step_dens

################
### State-dependent densities of TURNING ANGLE scaled by the equilibrium state densities

x <- seq(min(sharkdat$angle, na.rm=TRUE), 
         max(sharkdat$angle, na.rm=TRUE), length=1000)
wsm6a$mle$angle
rho1_angle <- wsm6a$mle$angle[2]
rho2_angle <- wsm6a$mle$angle[4]

# Work out state-dependent densities
# step
d1_angle <- (dwrappedcauchy(x, mu = circular(0), rho = rho1_angle))*delta.avg[1]
d2_angle <- (dwrappedcauchy(x, mu = circular(0), rho = rho2_angle))*delta.avg[2]

dmarg_angle <- d1_angle + d2_angle

# Plot densities for maximum dive depth, using manual colours or a colour palette like `viridis` or `RColorBrewer`
angle_dat <- sharkdat$angle

quartz(); 
angle_dens <- ggplot(data=data.frame(x=angle_dat), aes(x,..density..)) + 
  geom_histogram(boundary=0, binwidth=0.2, fill="grey90") + theme_bw() + 
  geom_line(data=data.frame(x=x, d1_angle=d1_angle), aes(x, d1_angle, colour="Resident"), size=1.5) +
  geom_line(data=data.frame(x=x, d2_angle=d2_angle), aes(x, d2_angle, colour="Transient"), size=1.5) +
  geom_line(data=data.frame(x=x, dmarg_angle=dmarg_angle), aes(x, dmarg_angle, color="Marginal"), size=1.2, linetype="dashed") +
  scale_colour_manual(name="Densities", values=c("Resident" = state_cols[1], "Transient" = state_cols[2],
                                                 "Marginal" =  marg_col), 
                      aesthetics=c("colour"), breaks=c("Resident","Transient","Marginal")) + 
  xlab("Turning angle") + ylab("") + 
  theme(legend.position="none") + 
  theme(text=element_text(size=20)); angle_dens

# create sex labels
title_lab <- cowplot::ggdraw() + cowplot::draw_label("State-dependent densities", x = 0.35, y = 0.3, angle = 0,
                               vjust = 0, hjust = 0, size = 25); title_lab

# add labels to plot 
lay <- rbind(c(1,1,1,1,1,1,1,1),
             c(2,2,2,2,3,3,3,3),
             c(2,2,2,2,3,3,3,3),
             c(2,2,2,2,3,3,3,3),
             c(2,2,2,2,3,3,3,3),
             c(2,2,2,2,3,3,3,3),
             c(2,2,2,2,3,3,3,3))

quartz()
density_comp_plot <- grid.arrange(title_lab, step_dens, angle_dens,
                                     layout_matrix=lay)

ggsave(filename=here::here("figures","step_dep_densities.jpg"), 
       plot=density_comp_plot,
       width=30, height=16, units="cm",dpi=700)



###################################################################
# ~~~~~~~~~~~~~~~~ Tracks faceted by sex and state ~~~~~~~~~~~~~~~~
###################################################################

# Tracks faceted by sex (rows) and state (columns)
quartz()
wsm6a_state_map_bysexstate <- wsm6a_preds_sf %>% 
  ggplot() + 
  geom_sf(data=safr_countries) +
  geom_sf(aes(colour = wsm6a_vit), size=pointsize) + 
  theme_bw(base_size=base_size) +
  scale_colour_viridis_d(name="State", begin=0.05, end=0.65,
                         labels=stateNames, alpha=alpha.trans) +
  theme(strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(size = 15)) +
  facet_wrap(~ sex + wsm6a_vit) + 
  theme(legend.position = "none", 
        strip.text.x = element_blank(),
        strip.background = element_blank()); wsm6a_state_map_bysexstate

# create sex labels
F_lab <- cowplot::ggdraw() + cowplot::draw_label("Females", x = 0, y = 0.6, angle = 0,
                                 vjust = 0, hjust = 0, size = 25); F_lab
M_lab <- cowplot::ggdraw() + cowplot::draw_label("Males", x = 0, y = 0.6, angle = 0,
                                 vjust = 0, hjust = 0, size = 25); M_lab

# create common x and y labels
x.lab <- textGrob("Longitude", vjust = 0, hjust = 1,
                  gp=gpar(fontsize=18))

y.lab <- textGrob("Latitude", vjust = 1, hjust = 0,
                  gp=gpar(fontsize=18), rot=90)

# add labels to plot
lay <- rbind(c(1,1,1,1,1,1,1,2),
             c(1,1,1,1,1,1,1,3),
             c(1,1,1,1,1,1,1,3),
             c(1,1,1,1,1,1,1,3),
             c(1,1,1,1,1,1,1,4),
             c(1,1,1,1,1,1,1,4),
             c(1,1,1,1,1,1,1,4))

statebysex_grid_comp_plot <- grid.arrange(wsm6a_state_map_bysexstate, 
                                          SexLeg_plot,
                                          F_lab,M_lab,
                                          layout_matrix=lay)

statebysex_grid_comp_plotL <- grid.arrange(arrangeGrob(statebysex_grid_comp_plot, 
                                                       left = y.lab, bottom = x.lab))                         

ggsave(filename=here::here("figures","map_grid_state_by_sex.jpg"),
      plot=statebysex_grid_comp_plotL,
      width=30, height=20, units="cm",dpi=700)



###################################################################
# ~~~~~~ Tracks coloured by P(Transient) and faceted by sex ~~~~~~~
###################################################################

# Plot of P(Transient) faceted by sex
quartz()
wsm6a_stateprobs_bysex_map_12hr <- wsm6a_preds_sf %>% 
  ggplot() + 
  geom_sf(data=safr_countries) +
  #geom_sf(aes(colour=id), show.legend=FALSE)
  geom_sf(aes(colour = wsm6a_Ptrans), alpha=alpha.trans+0.1, size=pointsize) + 
  theme_bw(base_size=base_size) +
  scale_colour_viridis_c(name="P(Transient)", begin=0.05, end=0.95, limits=c(0,1)) +
  theme(strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(size = base_size)) +
  facet_wrap(~ sex); wsm6a_stateprobs_bysex_map_12hr
#ggsave(plot=wsm6_stateprobs_map_12hr, filename="wsm6_stateprobs_map_12hr.png",path=here("output/ws_HMM"), width=20, height=10, dpi=700)

# create common x and y labels
x.lab <- textGrob("Longitude", vjust = -6, hjust=1,
                  gp=gpar(fontsize=18))

y.lab <- textGrob("Latitude", vjust = 1,
                  gp=gpar(fontsize=18), rot=90)

pTransient_bysex_comp_plotL <- grid.arrange(arrangeGrob(wsm6a_stateprobs_bysex_map_12hr, left = y.lab, bottom = x.lab))                         

ggsave(filename=here::here("output/ws_HMM","map_pTransient_bysex.jpg"), 
       plot=pTransient_bysex_comp_plotL,
       width=30, height=17, units="cm",dpi=700)

##########################################################################
# ~~~ Tracks coloured by P(Transient) and faceted by sex and maturity  ~~~
##########################################################################

# Tracks faceted by sex (rows: F then M) and maturity (columns: Juvenile Sub_adult Adult) 
# and coloured by P(Transient) - with legend

pointsize <- 1.4
ls_lab_size <- 15

# Main plot
wsm6a_stateprobs_bymaturity_map <- wsm6a_preds_sf %>% 
  ggplot() + theme_bw() +
  xlab("Longitude") + ylab("Latitude") +
  geom_sf(data=safr_countries) +
  geom_sf(aes(colour = wsm6a_Ptrans), alpha=alpha.trans+0.1, size=pointsize) + 
  scale_colour_viridis_c(name="P(Transient)", limits=c(0,1),
                         begin = 0, end = 0.94) +
  facet_wrap(sex ~ maturity) +
  theme(strip.background = element_rect(colour=NA, fill=NA), 
        #strip.text.x = element_text(size = 15),
        strip.text.x = element_blank(),
        axis.text = element_text(size=10),
        legend.position = "none",
        plot.margin = unit(c(t = -3, r = -0.5, b = 1, l = 1), "cm"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0),
                                    size = ls_lab_size),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0),
                                    size = ls_lab_size)); wsm6a_stateprobs_bymaturity_map

# Legend plot
wsm6a_stateprobs_bymaturity_map_LEG <- wsm6a_preds_sf %>% 
  ggplot() + 
  geom_sf(data=safr_countries) +
  #geom_sf(aes(colour=id), show.legend=FALSE)
  geom_sf(aes(colour = wsm6a_Ptrans), alpha=alpha.trans*2, size=pointsize) + 
  theme_bw(base_size=15) +
  scale_colour_viridis_c(name="P(Transient)", limits=c(0,1),
                         begin = 0, end = 0.94) +
  theme(strip.background = element_rect(colour=NA, fill=NA), 
        strip.text.x = element_text(size = 15),
        legend.position = c(0.01,0.2)) +
  facet_wrap(sex ~ maturity); wsm6a_stateprobs_bymaturity_map_LEG 

#' Get point legend
pTransLeg <- cowplot::get_legend(wsm6a_stateprobs_bymaturity_map_LEG)
pTransLeg_plot <- ggdraw(pTransLeg) 

# create sex labels
F_lab <- ggdraw() + draw_label("Females", x = 0.15, y = 0.3, angle = 0,
                               vjust = 0, hjust = 0, size = 20); F_lab
M_lab <- ggdraw() + draw_label("Males", x = 0.15, y = 0.5, angle = 0,
                               vjust = 0, hjust = 0, size = 20); M_lab

#' create title plots
Title_plotJuv <- ggdraw() + draw_text("Juvenile", x = 0.95, y = 0.49, 
                                      vjust = 1, hjust = 1, size = 20); Title_plotJuv
Title_plotSubad <- ggdraw() + draw_text("Subadult", x = 0.88, y = 0.49,
                                      vjust = 1, hjust = 1, size = 20); Title_plotSubad
Title_plotAd <- ggdraw() + draw_text("Adult", x = 0.7, y = 0.49,
                                      vjust = 1, hjust = 1, size = 20); Title_plotAd

lay <- rbind(c(1,1,1,2,2,2,3,3,3,4,4,8),
             c(rep(5,times=9),4,4,8),
             c(rep(5,times=9),4,4,8),
             c(rep(5,times=9),6,6,8),
             c(rep(5,times=9),6,6,8),
             c(rep(5,times=9),6,6,8),
             c(rep(5,times=9),7,7,9),
             c(rep(5,times=9),7,7,9))

map_grid_pTrans_comp_plot <- grid.arrange(Title_plotJuv, Title_plotSubad, Title_plotAd, F_lab, 
                                          wsm6a_stateprobs_bymaturity_map, M_lab, nullGrob(),
                                          pTransLeg, nullGrob(), layout_matrix=lay)

ggsave(filename=here::here("figures","pTrans_by_sex_and_maturity.jpg"), 
       plot=map_grid_pTrans_comp_plot,
       width=30, height=17, units="cm",dpi=700)

#
# NOW GO TO SCRIPT: 03ARGOSshark_plot_results.R
#






