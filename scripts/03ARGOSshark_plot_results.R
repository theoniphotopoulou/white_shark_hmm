# In this script the results from a hidden Markov model, implemented in the 
# momentuHMM R package by McClintock & Michelot (2018), and fitted to 
# time-regular model-corrected ARGOS satellite tracking data from 33 white 
# are plotted. The locations used in the model have been modelled and regularised using
# the foieGras R package by Jonsen et al. (2019).

# These data and the associated results are presented in detail in 
# "sex and size influence the spatiotemporal distribution of white sharks, 
# with implications for interactions with fisheries and spatial management 
# in the southwest Indian Ocean" by Alison Kock, Amanda T Lombard, Ryan Daly, 
# Victoria Goodall, Michael Meÿer, Ryan Johnson, Chris Fischer, Pieter Koen, 
# Dylan Irion, Enrico Gennari, Alison Towner, Oliver Jewell, Charlene da Silva, 
# Matt Dicken, Malcolm Smale, Theoni Photopoulou (2021). 

# The data custodian for these data is Alison Kock and 
# the code was written by Theoni Photopoulou (20211017)

require(momentuHMM)
require(lubridate)
require(dplyr)
require(rgdal)
require(ggplot2)
require(viridis)
require(ggmap)
require(cowplot)
require(gridExtra)
require(grid)
require(pgirmess)
require(Hmisc)
require(MASS)
require(numDeriv)
require(foreach)
require(tidyr)
require(purrr)
require(foieGras)
devtools::install_github("ropenscilabs/rnaturalearthdata")
devtools::install_github("ropenscilabs/rnaturalearthhires")
require(rnaturalearth)
require(sf)
require(here)

here()
Sys.setenv(TZ='UTC')

# Graphics windows are generated using quartz(). If you are
# using a Windows machine please replace quartz() with quartz()

# Load predictions
load(file=here::here("output","wshark_pred_locs_states.RData"))
# Load model object for plotting
load(file=here::here("output","wsm6a.RData"))

# Work out TL summary statistics
wshark_locs_states %>% group_by(sex, shark_id) %>% 
  summarise(unique(sex)) %>%
  count() # 20 Female, 13 Male
ttt <- wshark_locs_states %>% group_by(sex) %>% 
  summarise(minTL=min(TL), maxTL=max(TL),
            meanTL=mean(TL), sdTL=sd(TL)) 
ttt %>%
  mutate(seTL=sdTL/sqrt(c(20,13))) %>%
  print(n=Inf)

# Look at numbers of animals in each maturity levels
wshark_locs_states %>% group_by(sex, maturity) %>%
  summarise(n_maturity=n_distinct(shark_id))

# Specify plotting parameters
alpha.trans <- 0.2
base_size <- 20
pointsize <- 1

# How many locations are there per year?
table(wshark_locs_states$year)
# How many active tags are there per year?
wshark_locs_states %>% 
  group_by(year) %>% 
  summarise(n=n_distinct(shark_id)) %>% 
  arrange(year) 

# maturity specifications based on size
## Juvenile (>1.75 - 3.0 m)
## Sub-adult (Male: >3 - 3.6 m; Female: >3 - 4.8 m)
## Adult (Male: >3.6 m; Female: >4.8 m)

# Create dataframes with maturity sizes for use in the plots below
  # Adult sizes
maxTL_adM <- wshark_locs_states %>% filter(sex=="Male", TL>360) %>% 
  summarise(minTL=min(TL), maxTL=max(TL)) %>% as.numeric() %>% max()
maxTL_adF <- wshark_locs_states %>% filter(sex=="Female", TL>480) %>% 
  summarise(minTL=min(TL), maxTL=max(TL)) %>% as.numeric() %>% max()

maxTL_adM; maxTL_adF
adultM_TL <- maxTL_adM; adultF_TL <- maxTL_adF
    
  # observed Adult sizes Female and Male
wshark_locs_states %>% 
  filter(sex=="Female", maturity=="Adult") %>% 
  distinct(TL)
wshark_locs_states %>% 
  filter(sex=="Male", maturity=="Adult") %>% 
  distinct(TL) %>%
  summarise(Min_observed_length=min(TL), 
            Max_observed_length=max(TL),
            Mean_observed_length=mean(TL),
            Median_observed_length=median(TL))

  # Sub_adult sizes
meanTL_subadM <- wshark_locs_states %>% filter(sex=="Male", TL %in% c(300:360)) %>% 
  summarise(minTL=min(TL), maxTL=max(TL)) %>% as.numeric() %>% mean(.)
meanTL_subadF <- wshark_locs_states %>% filter(sex=="Female", TL %in% c(300:480)) %>% 
  summarise(minTL=min(TL), maxTL=max(TL)) %>% as.numeric() %>% mean(.)

meanTL_subadM; meanTL_subadF
subadultM_TL <- meanTL_subadM; subadultF_TL <- meanTL_subadF 
      
  
# observed Sub_adult sizes Female and Male
wshark_locs_states %>% filter(sex=="Female", maturity=="Sub_adult") %>% 
  distinct(TL) %>%
  summarise(Min_observed_length=min(TL), 
            Max_observed_length=max(TL),
            Mean_observed_length=mean(TL),
            Median_observed_length=median(TL)) 
wshark_locs_states %>% filter(sex=="Male", maturity=="Sub_adult") %>% 
  distinct(TL) %>%
  summarise(Min_observed_length=min(TL), 
            Max_observed_length=max(TL),
            Mean_observed_length=mean(TL),
            Median_observed_length=median(TL)) 


  # Juvenile size
meanTL_juven <- wshark_locs_states %>% filter(TL<300) %>% 
  summarise(minTL=min(TL), maxTL=max(TL)) %>% as.numeric() %>% mean(.)

meanTL_juven
juvenile_TL <- meanTL_juven;
  
  # observed Juvenile sizes Female and Male
wshark_locs_states %>% filter(sex=="Female", maturity=="Juvenile") %>% 
  distinct(TL) %>%
  summarise(Min_observed_length=min(TL), 
            Max_observed_length=max(TL),
            Mean_observed_length=mean(TL),
            Median_observed_length=median(TL)) 
wshark_locs_states %>% filter(sex=="Male", maturity=="Juvenile") %>% 
  distinct(TL) %>%
  summarise(Min_observed_length=min(TL), 
            Max_observed_length=max(TL),
            Mean_observed_length=mean(TL),
            Median_observed_length=median(TL)) 


# Covariate dataframes for sex and size
covs_adF <- data.frame(sex="Female", TL=adultF_TL, yearday=1)
covs_adM <- data.frame(sex="Male", TL=adultM_TL, yearday=1)

covs_subadF <- data.frame(sex="Female", TL=meanTL_subadF, yearday=1)
covs_subadM <- data.frame(sex="Male", TL=meanTL_subadM, yearday=1)

covs_juvF <- data.frame(sex="Female", TL=meanTL_juven, yearday=1)
covs_juvM <- data.frame(sex="Male", TL=meanTL_juven, yearday=1)

# Define a range of values for your covariates. This is to visualise the 
# HMM results in two different ways:
#  1) By sex and size across the whole year: 
#     adult Female, adult Male, subadult Female, subadult Male, juvenile Female, juvenile Male
#  2) By time of year for each sex across size range: 
#     day 40 (Feb 9, for non-leap years) and day 220 (Aug 8, for non-leap years)
#     I chose these by looking at the high and low in seasonal cycle: e.g., which.max(probs_matF[,1])

# 1) ACROSS THE WHOLE YEAR
lengthout <- 365
# adult Female
covs_adF <- data.frame(sex=as.factor(rep("Female",lengthout)), 
                       TL=rep(adultF_TL,lengthout), 
                       cosyd=cos(2*pi*seq(1,365,1)/365), 
                       sinyd=sin(2*pi*seq(1,365,1)/365))
levels(covs_adF$sex) <- c("Female","Male"); table(covs_adF$sex)
# adult Male
covs_adM <- data.frame(sex=as.factor(rep("Male",lengthout)), 
                       TL=rep(adultM_TL,lengthout), 
                       cosyd=cos(2*pi*seq(1,365,1)/365), 
                       sinyd=sin(2*pi*seq(1,365,1)/365))
levels(covs_adM$sex) <- c("Male","Female"); table(covs_adM$sex)
covs_adM$sex <- relevel(covs_adM$sex, ref="Female"); table(covs_adM$sex)
# subadult Female
covs_subadF <- data.frame(sex=as.factor(rep("Female",lengthout)), 
                          TL=rep(meanTL_subadF,lengthout), 
                          cosyd=cos(2*pi*seq(1,365,1)/365), 
                          sinyd=sin(2*pi*seq(1,365,1)/365))
levels(covs_subadF$sex) <- c("Female","Male"); table(covs_subadF$sex)
# subadult Male
covs_subadM <- data.frame(sex=as.factor(rep("Male",lengthout)), 
                          TL=rep(meanTL_subadM,lengthout), 
                          cosyd=cos(2*pi*seq(1,365,1)/365), 
                          sinyd=sin(2*pi*seq(1,365,1)/365))
levels(covs_subadM$sex) <- c("Male","Female"); table(covs_subadM$sex)
covs_subadM$sex <- relevel(covs_subadM$sex, ref="Female"); table(covs_subadM$sex)
# juvenile Female
covs_juvF <- data.frame(sex=as.factor(rep("Female",lengthout)), 
                        TL=rep(meanTL_juven,lengthout), 
                        cosyd=cos(2*pi*seq(1,365,1)/365), 
                        sinyd=sin(2*pi*seq(1,365,1)/365))
levels(covs_juvF$sex) <- c("Female","Male"); table(covs_juvF$sex)
# juvenile Male
covs_juvM <- data.frame(sex=as.factor(rep("Male",lengthout)), 
                        TL=rep(meanTL_juven,lengthout), 
                        cosyd=cos(2*pi*seq(1,365,1)/365), 
                        sinyd=sin(2*pi*seq(1,365,1)/365))
levels(covs_juvM$sex) <- c("Male","Female"); table(covs_juvM$sex)
covs_juvM$sex <- relevel(covs_juvM$sex, ref="Female"); table(covs_juvM$sex)

dim(covs_adF); head(covs_adF)
dim(covs_adM); head(covs_adM)
dim(covs_subadF); head(covs_subadF)
dim(covs_subadM); head(covs_subadM)
dim(covs_juvF); head(covs_juvF)
dim(covs_juvM); head(covs_juvM)

# Rename model for brevity and extract model info
m <- wsm6a 
nbStates <- length(wsm6a$stateNames)
betaMLE <- m$mle$beta

# Design matrix
desMat_adF <- model.matrix(m$conditions$formula, data = covs_adF)
desMat_adM <- model.matrix(m$conditions$formula, data = covs_adM)
desMat_subadF <- model.matrix(m$conditions$formula, data = covs_subadF)
desMat_subadM <- model.matrix(m$conditions$formula, data = covs_subadM)
desMat_juvF <- model.matrix(m$conditions$formula, data = covs_juvF)
desMat_juvM <- model.matrix(m$conditions$formula, data = covs_juvM)

probs_adF <- stationary(m, covs=desMat_adF)[[1]]
probs_adM <- stationary(m, covs=desMat_adM)[[1]]
probs_subadF <- stationary(m, covs=desMat_subadF)[[1]]
probs_subadM <- stationary(m, covs=desMat_subadM)[[1]]
probs_juvF <- stationary(m, covs=desMat_juvF)[[1]]
probs_juvM <- stationary(m, covs=desMat_juvM)[[1]]

# Covariance matrix of estimates
Sigma <- ginv(m$mod$hessian)
dim(Sigma)

# Indices corresponding to regression coefficients in m$mod$estimate
gamInd <- 7:16
# to get the right index, count the number of parameters you 
# are estimating for the density distributions, that is where your 
# gammas will start. then count the number of elements that contribute to 
# your tpm (coefficients in your linear predictor) and that's where 
# your gammas will start.

# Quantile for confidence interval (1.96 for 95% CI)
quantSup <- 1.96

# Input: beta, Output: delta
# for differentiation in delta method below
get_stat <- function(beta, covs, nbStates, i) {
  gamma <- moveHMM:::trMatrix_rcpp(nbStates, beta, covs)[,,1] # get gamma from beta
  solve(t(diag(nbStates) - gamma + 1), rep(1, nbStates))[i] # get delta from gamma
}

# Loop over states
lci_adF <- matrix(NA, lengthout, nbStates); lci_adM <- matrix(NA, lengthout, nbStates)
uci_adF <- matrix(NA, lengthout, nbStates); uci_adM <- matrix(NA, lengthout, nbStates)
lci_subadF <- matrix(NA, lengthout, nbStates); lci_subadM <- matrix(NA, lengthout, nbStates)
uci_subadF <- matrix(NA, lengthout, nbStates); uci_subadM <- matrix(NA, lengthout, nbStates)
lci_juvF <- matrix(NA, lengthout, nbStates); lci_juvM <- matrix(NA, lengthout, nbStates)
uci_juvF <- matrix(NA, lengthout, nbStates); uci_juvM <- matrix(NA, lengthout, nbStates)

for(state in 1:nbStates) {
  # Get gradient of get_stat function
  dN_adF <- t(apply(desMat_adF, 1, function(x)
    grad(get_stat, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state)))
  dN_adM <- t(apply(desMat_adM, 1, function(x)
    grad(get_stat, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state)))
  dN_subadF <- t(apply(desMat_subadF, 1, function(x)
    grad(get_stat, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state)))
  dN_subadM <- t(apply(desMat_subadM, 1, function(x)
    grad(get_stat, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state)))
  dN_juvF <- t(apply(desMat_juvF, 1, function(x)
    grad(get_stat, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state)))
  dN_juvM <- t(apply(desMat_juvM, 1, function(x)
    grad(get_stat, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state)))
  
  # Standard errors from delta method formula
  se_adF <- t(apply(dN_adF, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
  se_adM <- t(apply(dN_adM, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
  se_subadF <- t(apply(dN_subadF, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
  se_subadM <- t(apply(dN_subadM, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
  se_juvF <- t(apply(dN_juvF, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
  se_juvM <- t(apply(dN_juvM, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
  
  # Lower and upper bounds of confidence interval
  lci_adF[,state] <- plogis(qlogis(probs_adF[,state]) - quantSup*se_adF/(probs_adF[,state]-probs_adF[,state]^2))
  uci_adF[,state] <- plogis(qlogis(probs_adF[,state]) + quantSup*se_adF/(probs_adF[,state]-probs_adF[,state]^2))
  
  lci_adM[,state] <- plogis(qlogis(probs_adM[,state]) - quantSup*se_adM/(probs_adM[,state]-probs_adM[,state]^2))
  uci_adM[,state] <- plogis(qlogis(probs_adM[,state]) + quantSup*se_adM/(probs_adM[,state]-probs_adM[,state]^2))

  lci_subadF[,state] <- plogis(qlogis(probs_subadF[,state]) - quantSup*se_subadF/(probs_subadF[,state]-probs_subadF[,state]^2))
  uci_subadF[,state] <- plogis(qlogis(probs_subadF[,state]) + quantSup*se_subadF/(probs_subadF[,state]-probs_subadF[,state]^2))
  
  lci_subadM[,state] <- plogis(qlogis(probs_subadM[,state]) - quantSup*se_subadM/(probs_subadM[,state]-probs_subadM[,state]^2))
  uci_subadM[,state] <- plogis(qlogis(probs_subadM[,state]) + quantSup*se_subadM/(probs_subadM[,state]-probs_subadM[,state]^2))

  lci_juvF[,state] <- plogis(qlogis(probs_juvF[,state]) - quantSup*se_subadF/(probs_subadF[,state]-probs_subadF[,state]^2))
  uci_juvF[,state] <- plogis(qlogis(probs_juvF[,state]) + quantSup*se_subadF/(probs_subadF[,state]-probs_subadF[,state]^2))
  
  lci_juvM[,state] <- plogis(qlogis(probs_juvM[,state]) - quantSup*se_juvM/(probs_juvM[,state]-probs_juvM[,state]^2))
  uci_juvM[,state] <- plogis(qlogis(probs_juvM[,state]) + quantSup*se_juvM/(probs_juvM[,state]-probs_juvM[,state]^2))
}

# Plot state probs and confidence intervals
quartz()
pal <- c("firebrick", "royalblue")
plot(NA, xlim = c(1,365), ylim = c(0, 1))
for(state in 1:nbStates) {
  points(seq(1,365,1), probs_adF[,state], type = "l", col = pal[state])
  points(seq(1,365,1), lci_adF[,state], type = "l", lty = 2, col = pal[state])
  points(seq(1,365,1), uci_adF[,state], type = "l", lty = 2, col = pal[state])
}

quartz()
pal <- c("firebrick", "royalblue")
plot(NA, xlim = c(1,365), ylim = c(0, 1))
for(state in 1:nbStates) {
  points(seq(1,365,1), probs_adM[,state], type = "l", col = pal[state])
  points(seq(1,365,1), lci_adM[,state], type = "l", lty = 2, col = pal[state])
  points(seq(1,365,1), uci_adM[,state], type = "l", lty = 2, col = pal[state])
}

quartz()
pal <- c("firebrick", "royalblue")
plot(NA, xlim = c(1,365), ylim = c(0, 1))
for(state in 1:nbStates) {
  points(seq(1,365,1), probs_subadF[,state], type = "l", col = pal[state])
  points(seq(1,365,1), lci_subadF[,state], type = "l", lty = 2, col = pal[state])
  points(seq(1,365,1), uci_subadF[,state], type = "l", lty = 2, col = pal[state])
}

quartz()
pal <- c("firebrick", "royalblue")
plot(NA, xlim = c(1,365), ylim = c(0, 1))
for(state in 1:nbStates) {
  points(seq(1,365,1), probs_subadM[,state], type = "l", col = pal[state])
  points(seq(1,365,1), lci_subadM[,state], type = "l", lty = 2, col = pal[state])
  points(seq(1,365,1), uci_subadM[,state], type = "l", lty = 2, col = pal[state])
}

#################### sanity check - compare with your "home made" plots of the stationary distributions
CIbeta(m)
quartz()
# winter
plotStationary(m, covs=data.frame(sex="M", TL=420, yearday=300), plotCI=TRUE)
plotStationary(wsm6, covs=data.frame(sex="M", TL=320, yearday=300), plotCI=TRUE)
plotStationary(wsm6, covs=data.frame(sex="M", TL=264, yearday=300), plotCI=TRUE)
#
plotStationary(wsm6, covs=data.frame(sex="F", TL=505, yearday=300), plotCI=TRUE)
plotStationary(wsm6, covs=data.frame(sex="F", TL=399, yearday=300), plotCI=TRUE)
plotStationary(wsm6, covs=data.frame(sex="F", TL=264, yearday=300), plotCI=TRUE)

# summer
plotStationary(wsm6, covs=data.frame(sex="M", TL=420, yearday=30), plotCI=TRUE)
plotStationary(wsm6, covs=data.frame(sex="M", TL=320, yearday=30), plotCI=TRUE)
plotStationary(wsm6, covs=data.frame(sex="M", TL=264, yearday=30), plotCI=TRUE)

plotStationary(wsm6, covs=data.frame(sex="F", TL=505, yearday=30), plotCI=TRUE)
plotStationary(wsm6, covs=data.frame(sex="F", TL=399, yearday=30), plotCI=TRUE)
plotStationary(wsm6, covs=data.frame(sex="F", TL=264, yearday=30), plotCI=TRUE)
###################

head(wshark_locs_states)
table(wshark_locs_states$viterbi_state)
stateNames <- c("Resident", "Transient")

# 2) ACROSS THE RANGE OF SHARK SIZE
# The range of TL is bigger for females than males
wshark_locs_states %>% filter(sex=="Male") %>% summarise(minTL=min(TL), maxTL=max(TL)) %>% as.numeric -> TLrangeM 
wshark_locs_states %>% filter(sex=="Female") %>% summarise(minTL=min(TL), maxTL=max(TL)) %>% as.numeric -> TLrangeF
 
lengthoutM <- TLrangeM[2]-TLrangeM[1]+1
lengthoutF <- TLrangeF[2]-TLrangeF[1]+1

summer <- 40 # julian day 40 (9 Feb for non-leap years)
winter <- 220 # julian day 220 (8 Aug for non-leap years)

# late summer Female
covs_sumF <- data.frame(sex=as.factor(rep("Female",lengthoutF)), 
                        TL=seq(TLrangeF[1],TLrangeF[2],1), 
                        cosyd=cos(2*pi*summer/365), 
                        sinyd=sin(2*pi*summer/365))
levels(covs_sumF$sex) <- c("Female","Male"); table(covs_sumF$sex)
# late summer Male
covs_sumM <- data.frame(sex=as.factor(rep("Male",lengthoutM)), 
                        TL=seq(TLrangeM[1],TLrangeM[2],1), 
                        cosyd=cos(2*pi*summer/365), 
                        sinyd=sin(2*pi*summer/365))
levels(covs_sumM$sex) <- c("Male","Female"); table(covs_sumM$sex)
covs_sumM$sex <- relevel(covs_sumM$sex, ref="Female"); table(covs_sumM$sex)
# late winter Female
covs_winF <- data.frame(sex=as.factor(rep("Female",lengthoutF)), 
                        TL=seq(TLrangeF[1],TLrangeF[2],1), 
                        cosyd=cos(2*pi*winter/365), 
                        sinyd=sin(2*pi*winter/365))
levels(covs_winF$sex) <- c("Female","Male"); table(covs_winF$sex)
# late winter Male
covs_winM <- data.frame(sex=as.factor(rep("Male",lengthoutM)), 
                        TL=seq(TLrangeM[1],TLrangeM[2],1), 
                        cosyd=cos(2*pi*winter/365), 
                        sinyd=sin(2*pi*winter/365))
levels(covs_winM$sex) <- c("Male","Female"); table(covs_winM$sex)
covs_winM$sex <- relevel(covs_winM$sex, ref="Female"); table(covs_winM$sex)

dim(covs_sumF); head(covs_sumF)
dim(covs_sumM); head(covs_sumM)
dim(covs_winF); head(covs_winF)
dim(covs_winM); head(covs_winM)

# Design matrix
desMat_sumF <- model.matrix(m$conditions$formula, data = covs_sumF)
desMat_sumM <- model.matrix(m$conditions$formula, data = covs_sumM)
desMat_winF <- model.matrix(m$conditions$formula, data = covs_winF)
desMat_winM <- model.matrix(m$conditions$formula, data = covs_winM)

probs_sumF <- stationary(m, covs=desMat_sumF)[[1]]
probs_sumM <- stationary(m, covs=desMat_sumM)[[1]]
probs_winF <- stationary(m, covs=desMat_sumF)[[1]]
probs_winM <- stationary(m, covs=desMat_sumM)[[1]]

# Covariance matrix of estimates
Sigma <- ginv(m$mod$hessian)
dim(Sigma)

# Indices corresponding to regression coefficients in m$mod$estimate
gamInd <- 7:16
# to get the right index, count the number of parameters you 
# are estimating for the density distributions, that is where your 
# gammas will start. then count the number of elements that contribute to 
# your tpm (coefficients in your linear predictor) and that's where 
# your gammas will start.

# Quantile for confidence interval (1.96 for 95% CI)
quantSup <- 1.96

# Input: beta, Output: delta
# for differentiation in delta method below
get_stat <- function(beta, covs, nbStates, i) {
  gamma <- moveHMM:::trMatrix_rcpp(nbStates, beta, covs)[,,1] # get gamma from beta
  solve(t(diag(nbStates) - gamma + 1), rep(1, nbStates))[i] # get delta from gamma
}

# Loop over states
lci_sumF <- matrix(NA, lengthoutF, nbStates); lci_sumM <- matrix(NA, lengthoutM, nbStates)
uci_sumF <- matrix(NA, lengthoutF, nbStates); uci_sumM <- matrix(NA, lengthoutM, nbStates)
lci_winF <- matrix(NA, lengthoutF, nbStates); lci_winM <- matrix(NA, lengthoutM, nbStates)
uci_winF <- matrix(NA, lengthoutF, nbStates); uci_winM <- matrix(NA, lengthoutM, nbStates)

for(state in 1:nbStates) {
  # Get gradient of get_stat function
  dN_sumF <- t(apply(desMat_sumF, 1, function(x)
    grad(get_stat, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state)))
  dN_sumM <- t(apply(desMat_sumM, 1, function(x)
    grad(get_stat, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state)))
  dN_winF <- t(apply(desMat_winF, 1, function(x)
    grad(get_stat, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state)))
  dN_winM <- t(apply(desMat_winM, 1, function(x)
    grad(get_stat, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state)))
  
  # Standard errors from delta method formula
  se_sumF <- t(apply(dN_sumF, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
  se_sumM <- t(apply(dN_sumM, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
  se_winF <- t(apply(dN_winF, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
  se_winM <- t(apply(dN_winM, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
  
  # Lower and upper bounds of confidence interval
  lci_sumF[,state] <- plogis(qlogis(probs_sumF[,state]) - quantSup*se_sumF/(probs_sumF[,state]-probs_sumF[,state]^2))
  uci_sumF[,state] <- plogis(qlogis(probs_sumF[,state]) + quantSup*se_sumF/(probs_sumF[,state]-probs_sumF[,state]^2))
  
  lci_sumM[,state] <- plogis(qlogis(probs_sumM[,state]) - quantSup*se_sumM/(probs_sumM[,state]-probs_sumM[,state]^2))
  uci_sumM[,state] <- plogis(qlogis(probs_sumM[,state]) + quantSup*se_sumM/(probs_sumM[,state]-probs_sumM[,state]^2))
  
  lci_winF[,state] <- plogis(qlogis(probs_winF[,state]) - quantSup*se_winF/(probs_winF[,state]-probs_winF[,state]^2))
  uci_winF[,state] <- plogis(qlogis(probs_winF[,state]) + quantSup*se_winF/(probs_winF[,state]-probs_winF[,state]^2))
  
  lci_winM[,state] <- plogis(qlogis(probs_winM[,state]) - quantSup*se_winM/(probs_winM[,state]-probs_winM[,state]^2))
  uci_winM[,state] <- plogis(qlogis(probs_winM[,state]) + quantSup*se_winM/(probs_winM[,state]-probs_winM[,state]^2))
}

# Plot state probs and confidence intervals
quartz()
pal <- c("firebrick", "royalblue")
plot(NA, xlim = TLrangeF, ylim = c(0, 1))
for(state in 1:nbStates) {
  points(seq(TLrangeF[1],TLrangeF[2],length.out = lengthoutF), probs_sumF[,state], type = "l", col = pal[state])
  points(seq(TLrangeF[1],TLrangeF[2],length.out = lengthoutF), lci_sumF[,state], type = "l", lty = 2, col = pal[state])
  points(seq(TLrangeF[1],TLrangeF[2],length.out = lengthoutF), uci_sumF[,state], type = "l", lty = 2, col = pal[state])
}

quartz()
pal <- c("firebrick", "royalblue")
plot(NA, xlim = TLrangeM, ylim = c(0, 1))
for(state in 1:nbStates) {
  points(seq(TLrangeM[1],TLrangeM[2],length.out = lengthoutM), probs_sumM[,state], type = "l", col = pal[state])
  points(seq(TLrangeM[1],TLrangeM[2],length.out = lengthoutM), lci_sumM[,state], type = "l", lty = 2, col = pal[state])
  points(seq(TLrangeM[1],TLrangeM[2],length.out = lengthoutM), uci_sumM[,state], type = "l", lty = 2, col = pal[state])
}

quartz()
pal <- c("firebrick", "royalblue")
plot(NA, xlim = TLrangeF, ylim = c(0, 1))
for(state in 1:nbStates) {
  points(seq(TLrangeF[1],TLrangeF[2],length.out = lengthoutF), probs_winF[,state], type = "l", col = pal[state])
  points(seq(TLrangeF[1],TLrangeF[2],length.out = lengthoutF), lci_winF[,state], type = "l", lty = 2, col = pal[state])
  points(seq(TLrangeF[1],TLrangeF[2],length.out = lengthoutF), uci_winF[,state], type = "l", lty = 2, col = pal[state])
}

quartz()
pal <- c("firebrick", "royalblue")
plot(NA, xlim = TLrangeM, ylim = c(0, 1))
for(state in 1:nbStates) {
  points(seq(TLrangeM[1],TLrangeM[2],length.out = lengthoutM), probs_winM[,state], type = "l", col = pal[state])
  points(seq(TLrangeM[1],TLrangeM[2],length.out = lengthoutM), lci_winM[,state], type = "l", lty = 2, col = pal[state])
  points(seq(TLrangeM[1],TLrangeM[2],length.out = lengthoutM), uci_winM[,state], type = "l", lty = 2, col = pal[state])
}

###################################################################################
############ Create plot: SEX AND LENGTH VS DAY OF THE YEAR
###################################################################################

# Create dataframe for Male and Female by maturity, for plotting purposes
adF <- data.frame(sex=covs_adF$sex, TL=covs_adF$TL, yearday=seq(1,365,1),
                   Resident_low=lci_adF[,1], Resident_mle=probs_adF[,1], Resident_upp=uci_adF[,1],
                   Transient_low=lci_adF[,2], Transient_mle=probs_adF[,2], Transient_upp=uci_adF[,2])

adM <- data.frame(sex=covs_adM$sex, TL=covs_adM$TL, yearday=seq(1,365,1),
                   Resident_low=lci_adM[,1], Resident_mle=probs_adM[,1], Resident_upp=uci_adM[,1],
                   Transient_low=lci_adM[,2], Transient_mle=probs_adM[,2], Transient_upp=uci_adM[,2])

subadF <- data.frame(sex=covs_subadF$sex, TL=covs_subadF$TL, yearday=seq(1,365,1),
                   Resident_low=lci_subadF[,1], Resident_mle=probs_subadF[,1], Resident_upp=uci_subadF[,1],
                   Transient_low=lci_subadF[,2], Transient_mle=probs_subadF[,2], Transient_upp=uci_subadF[,2])

subadM <- data.frame(sex=covs_subadM$sex, TL=covs_subadM$TL, yearday=seq(1,365,1),
                   Resident_low=lci_subadM[,1], Resident_mle=probs_subadM[,1], Resident_upp=uci_subadM[,1],
                   Transient_low=lci_subadM[,2], Transient_mle=probs_subadM[,2], Transient_upp=uci_subadM[,2])

juvF <- data.frame(sex=covs_juvF$sex, TL=covs_juvF$TL, yearday=seq(1,365,1),
                     Resident_low=lci_juvF[,1], Resident_mle=probs_juvF[,1], Resident_upp=uci_juvF[,1],
                     Transient_low=lci_juvF[,2], Transient_mle=probs_juvF[,2], Transient_upp=uci_juvF[,2])

juvM <- data.frame(sex=covs_juvM$sex, TL=covs_juvM$TL, yearday=seq(1,365,1),
                     Resident_low=lci_juvM[,1], Resident_mle=probs_juvM[,1], Resident_upp=uci_juvM[,1],
                     Transient_low=lci_juvM[,2], Transient_mle=probs_juvM[,2], Transient_upp=uci_juvM[,2])

alpha.trans <- 0.2
mycols <- viridis_pal(begin=0.05, end=0.65, option="D")(2)
ymax <- 1.05; xmax <- 365
label_text_size <- 3
# Summer period: 1st Dec is day 335, 28th Feb is day 59
summer_rect1 <- data.frame(rxmin=1, rxmax=59, rymin=0, rymax=1)
summer_rect2 <- data.frame(rxmin=335, rxmax=365, rymin=0, rymax=1)
summer_label1 <- data.frame(xl=30, yl=1.05, ltext="summer")
summer_label2 <- data.frame(xl=340, yl=1.05, ltext="summer")
# Winter period: 1st Jun is day 152, 31st Aug is day 243
winter_rect <- data.frame(rxmin=152, rxmax=243, rymin=0, rymax=1)
winter_label <- data.frame(xl=197, yl=1.05, ltext="winter")

## Bands for the confidence intervals 
## ADULT FEMALES with legend
quartz()
adF_withLeg <- ggplot(adF) +
  
  # Season rectangles
  geom_rect(data=summer_rect1, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  geom_rect(data=summer_rect2, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  geom_rect(data=winter_rect, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  
  # Resident
  geom_line(aes(x=yearday, y=Resident_mle, colour="Resident")) + 
  geom_ribbon(aes(x=yearday, ymin=Resident_low, ymax=Resident_upp, 
                  fill="Resident"), alpha=alpha.trans) +
  # Transient
  geom_line(aes(x=yearday, y=Transient_mle, colour="Transient")) +
  geom_ribbon(aes(x=yearday, ymin=Transient_low, ymax=Transient_upp, 
                  fill="Transient"), alpha=alpha.trans) +
  
  # Lines showing season
  #geom_vline(xintercept = summer, linetype="dashed", colour="grey70") +
  #geom_vline(xintercept = winter, linetype="dashed", colour="grey70") +
  
  # Season labels
  geom_label(data=summer_label1, aes(x=xl, y=yl, label=ltext), size = label_text_size, label.size = NA) +
  geom_label(data=summer_label2, aes(x=xl, y=yl, label=ltext), size = label_text_size, label.size = NA) +
  geom_label(data=winter_label, aes(x=xl, y=yl, label=ltext), size = label_text_size, label.size = NA) +
  
  # Legends
  scale_colour_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                      labels=c("Resident","Transient")) + 
  scale_fill_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                    labels=c("Resident","Transient")) +
  ylim(0,ymax) + 
  xlim(0,xmax) +
  xlab("Day of the year") + 
  ylab("Stationary state probability") + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "bottom",
        legend.justification = c(0.5,0.99),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Legend plot
Leg <- cowplot::get_legend(adF_withLeg)
Leg_plot <- ggdraw(Leg) 

## ADULT FEMALES without legend
adF_plot <- ggplot(adF) +
  
  # Season rectangles
  geom_rect(data=summer_rect1, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  geom_rect(data=summer_rect2, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  geom_rect(data=winter_rect, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  
  # Resident
  geom_line(aes(x=yearday, y=Resident_mle, colour="Resident")) + 
  geom_ribbon(aes(x=yearday, ymin=Resident_low, ymax=Resident_upp, 
                  fill="Resident"), alpha=alpha.trans) +
  # Transient
  geom_line(aes(x=yearday, y=Transient_mle, colour="Transient")) +
  geom_ribbon(aes(x=yearday, ymin=Transient_low, ymax=Transient_upp, 
                  fill="Transient"), alpha=alpha.trans) +

  # Season labels
  geom_label(data=summer_label1, aes(x=xl, y=yl, label=ltext), size = label_text_size, label.size = NA) +
  geom_label(data=summer_label2, aes(x=xl, y=yl, label=ltext), size = label_text_size, label.size = NA) +
  geom_label(data=winter_label, aes(x=xl, y=yl, label=ltext), size = label_text_size, label.size = NA) +
  
  # Legends
  scale_colour_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                      labels=c("Resident","Transient")) + 
  scale_fill_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                    labels=c("Resident","Transient")) +
  ylim(0,ymax) + 
  xlim(0,xmax) +
  xlab("") + ylab("") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## ADULT MALES 
adM_plot <- ggplot(adM) +
  
  # Season rectangles
  geom_rect(data=summer_rect1, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  geom_rect(data=summer_rect2, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  geom_rect(data=winter_rect, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  
  # Resident
  geom_line(aes(x=yearday, y=Resident_mle, colour="Resident")) + 
  geom_ribbon(aes(x=yearday, ymin=Resident_low, ymax=Resident_upp, 
                  fill="Resident"), alpha=alpha.trans) +
  # Transient
  geom_line(aes(x=yearday, y=Transient_mle, colour="Transient")) +
  geom_ribbon(aes(x=yearday, ymin=Transient_low, ymax=Transient_upp, 
                  fill="Transient"), alpha=alpha.trans) +

  # Season labels
  geom_label(data=summer_label1, aes(x=xl, y=yl, label=ltext), size = label_text_size, label.size = NA) +
  geom_label(data=summer_label2, aes(x=xl, y=yl, label=ltext), size = label_text_size, label.size = NA) +
  geom_label(data=winter_label, aes(x=xl, y=yl, label=ltext), size = label_text_size, label.size = NA) +
  
  # Legends
  scale_colour_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                      labels=c("Resident","Transient")) + 
  scale_fill_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                    labels=c("Resident","Transient")) +
  ylim(0,ymax) + 
  xlim(0,xmax) +
  xlab("") + ylab("") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## SUBADULT FEMALES 
subadF_plot <- ggplot(subadF) +
  
  # Season rectangles
  geom_rect(data=summer_rect1, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  geom_rect(data=summer_rect2, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  geom_rect(data=winter_rect, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  
  # Resident
  geom_line(aes(x=yearday, y=Resident_mle, colour="Resident")) + 
  geom_ribbon(aes(x=yearday, ymin=Resident_low, ymax=Resident_upp, 
                  fill="Resident"), alpha=alpha.trans) +
  # Transient
  geom_line(aes(x=yearday, y=Transient_mle, colour="Transient")) +
  geom_ribbon(aes(x=yearday, ymin=Transient_low, ymax=Transient_upp, 
                  fill="Transient"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                      labels=c("Resident","Transient")) + 
  scale_fill_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                    labels=c("Resident","Transient")) +
  ylim(0,ymax) + 
  xlim(0,xmax) +
  xlab("") + ylab("") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## SUBADULT MALES 
subadM_plot <- ggplot(subadM) +
  
  # Season rectangles
  geom_rect(data=summer_rect1, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  geom_rect(data=summer_rect2, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  geom_rect(data=winter_rect, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  
  # Resident
  geom_line(aes(x=yearday, y=Resident_mle, colour="Resident")) + 
  geom_ribbon(aes(x=yearday, ymin=Resident_low, ymax=Resident_upp, 
                  fill="Resident"), alpha=alpha.trans) +
  # Transient
  geom_line(aes(x=yearday, y=Transient_mle, colour="Transient")) +
  geom_ribbon(aes(x=yearday, ymin=Transient_low, ymax=Transient_upp, 
                  fill="Transient"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                      labels=c("Resident","Transient")) + 
  scale_fill_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                    labels=c("Resident","Transient")) +
  ylim(0,ymax) + 
  xlim(0,xmax) +
  xlab("") + ylab("") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## JUVENILE FEMALES 
juvF_plot <- ggplot(juvF) +
  
  # Season rectangles
  geom_rect(data=summer_rect1, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  geom_rect(data=summer_rect2, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  geom_rect(data=winter_rect, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  
  # Resident
  geom_line(aes(x=yearday, y=Resident_mle, colour="Resident")) + 
  geom_ribbon(aes(x=yearday, ymin=Resident_low, ymax=Resident_upp, 
                  fill="Resident"), alpha=alpha.trans) +
  # Transient
  geom_line(aes(x=yearday, y=Transient_mle, colour="Transient")) +
  geom_ribbon(aes(x=yearday, ymin=Transient_low, ymax=Transient_upp, 
                  fill="Transient"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                      labels=c("Resident","Transient")) + 
  scale_fill_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                    labels=c("Resident","Transient")) +
  ylim(0,ymax) + 
  xlim(0,xmax) +
  xlab("") + ylab("") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## JUVENILE MALES 
juvM_plot <- ggplot(juvM) +
  
  # Season rectangles
  geom_rect(data=summer_rect1, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  geom_rect(data=summer_rect2, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  geom_rect(data=winter_rect, aes(xmin=rxmin, xmax=rxmax, ymin=rymin, ymax=rymax), 
            fill="grey70", alpha=0.2) +
  
  # Resident
  geom_line(aes(x=yearday, y=Resident_mle, colour="Resident")) + 
  geom_ribbon(aes(x=yearday, ymin=Resident_low, ymax=Resident_upp, 
                  fill="Resident"), alpha=alpha.trans) +
  # Transient
  geom_line(aes(x=yearday, y=Transient_mle, colour="Transient")) +
  geom_ribbon(aes(x=yearday, ymin=Transient_low, ymax=Transient_upp, 
                  fill="Transient"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                      labels=c("Resident","Transient")) + 
  scale_fill_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                    labels=c("Resident","Transient")) +
  ylim(0,ymax) + 
  xlim(0,xmax) +
  xlab("") + ylab("") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# create maturity labels
Ad_lab <- ggdraw() + draw_label("Adult", x = 0.4, y = 0.45, angle = -90,
                                  vjust = 6, hjust = 1.5, size = 15); Ad_lab
Subad_lab <- ggdraw() + draw_label("Subadult", x = 0.4, y = 0.37, angle = -90,
                                  vjust = 6, hjust = 1.2, size = 15); Subad_lab
Juv_lab <- ggdraw() + draw_label("Juvenile", x = 0.4, y = 0.39, angle = -90,
                                   vjust = 6, hjust = 1.2, size = 15); Juv_lab

# create common x and y labels
y.lab <- textGrob("Stationary state probabilities", 
                  vjust = 1.5, hjust = 0.4,
                  gp=gpar(fontsize=18), rot=90)

x.lab <- textGrob("Day of the year", 
                  vjust = -8, hjust = 0.9,
                  gp=gpar(fontsize=18))

# create title plots
Title_plotF <- ggdraw() + draw_text("Females", x = 0.85, y = 0.3,
                                    vjust = 2.3, hjust = 1.3, size = 15); Title_plotF
Title_plotM <- ggdraw() + draw_text("Males", x = 0.8, y = 0.3,
                                    vjust = 2.3, hjust = 1.3, size = 15); Title_plotM

# add labels to plot
lay <- rbind(c(1,1,1,2,2,2,3,3),
             c(5,5,5,6,6,6,7,7),
             c(5,5,5,6,6,6,7,7),
             c(8,8,8,9,9,9,10,10),
             c(8,8,8,9,9,9,10,10),
             c(11,11,11,12,12,12,13,13),
             c(11,11,11,12,12,12,13,13),
             c(14,14,14,14,14,14,14,15)
             )

yearday_comp_plot <- grid.arrange(Title_plotF, Title_plotM, nullGrob(),
                          adF_plot, adM_plot, Ad_lab,
                          subadF_plot, subadM_plot, Subad_lab,
                          juvF_plot, juvM_plot, Juv_lab, Leg_plot, 
                          layout_matrix=lay)

quartz()
yearday_comp_plotL <- grid.arrange(arrangeGrob(yearday_comp_plot, left = y.lab, bottom = x.lab))                         

ggsave(filename=here::here("figures","sex_vs_yearday_stationary.jpg"), 
       plot=yearday_comp_plotL,
       width=23, height=25, units="cm",dpi=700)

###################################################################################
############ Create plot: SEX AND SEASON VS LENGTH
###################################################################################

# Create dataframe for Male and Female by season, for plotting purposes
sumF <- data.frame(sex=covs_sumF$sex, TL=covs_sumF$TL, yearday=summer,
                   Resident_low=lci_sumF[,1], Resident_mle=probs_sumF[,1], Resident_upp=uci_sumF[,1],
                   Transient_low=lci_sumF[,2], Transient_mle=probs_sumF[,2], Transient_upp=uci_sumF[,2])

sumM <- data.frame(sex=covs_sumM$sex, TL=covs_sumM$TL, yearday=summer,
                   Resident_low=lci_sumM[,1], Resident_mle=probs_sumM[,1], Resident_upp=uci_sumM[,1],
                   Transient_low=lci_sumM[,2], Transient_mle=probs_sumM[,2], Transient_upp=uci_sumM[,2])

winF <- data.frame(sex=covs_winF$sex, TL=covs_winF$TL, yearday=winter,
                   Resident_low=lci_winF[,1], Resident_mle=probs_winF[,1], Resident_upp=uci_winF[,1],
                   Transient_low=lci_winF[,2], Transient_mle=probs_winF[,2], Transient_upp=uci_winF[,2])

winM <- data.frame(sex=covs_winM$sex, TL=covs_winM$TL, yearday=winter,
                   Resident_low=lci_winM[,1], Resident_mle=probs_winM[,1], Resident_upp=uci_winM[,1],
                   Transient_low=lci_winM[,2], Transient_mle=probs_winM[,2], Transient_upp=uci_winM[,2])

alpha.trans <- 0.2
mycols <- viridis_pal(begin=0.05, end=0.65, option="D")(2)
ymax <- 1.05; xxmin <- 225; xmax <- 510
label_text_size <- 3
# Female TL size classes: Juvenile <300, Sub-adult 300-480, Adult >480
fmjuv_rect <- data.frame(rxmin=225, rxmax=300, rymin=0, rymax=1)
fsubad_rect <- data.frame(rxmin=300, rxmax=480, rymin=0, rymax=1)
fad_rect <- data.frame(rxmin=480, rxmax=510, rymin=0, rymax=1)
fjuv_label <- data.frame(xl=265, yl=1.05, ltext="juvenile")
fsubad_label <- data.frame(xl=390, yl=1.05, ltext="sub-adult")
fad_label <- data.frame(xl=495, yl=1.05, ltext="adult")
# Female TL size classes: Juvenile <300, Sub-adult 300-360, Adult >360
msubad_rect <- data.frame(rxmin=300, rxmax=360, rymin=0, rymax=1)
mad_rect <- data.frame(rxmin=360, rxmax=510, rymin=0, rymax=1)
mjuv_label <- data.frame(xl=265, yl=1.05, ltext="juvenile")
msubad_label <- data.frame(xl=330, yl=1.05, ltext="sub-adult")
mad_label <- data.frame(xl=435, yl=1.05, ltext="adult")
# Fill colour for rectangles
jvfill <- "grey70"
safill <- "white"
adfill <- "grey70"

## Bands for the confidence intervals 
## FEMALES IN SUMMER with legend
quartz()
sumF_withLeg <- ggplot(sumF) +

  # Resident
  geom_line(aes(x=TL, y=Resident_mle, colour="Resident")) + 
  geom_ribbon(aes(x=TL, ymin=Resident_low, ymax=Resident_upp, 
                  fill="Resident"), alpha=alpha.trans) +
  # Transient
  geom_line(aes(x=TL, y=Transient_mle, colour="Transient")) +
  geom_ribbon(aes(x=TL, ymin=Transient_low, ymax=Transient_upp, 
                  fill="Transient"), alpha=alpha.trans) +
  
  # Legends
  scale_colour_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                      labels=c("Resident","Transient")) + 
  scale_fill_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                    labels=c("Resident","Transient")) +
  ylim(0,ymax) + xlim(xxmin,xmax) +
  xlab("Length (cm)") + 
  ylab("Stationary state probability") + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "bottom",
        legend.justification = c(0.55,1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Legend plot
Leg <- cowplot::get_legend(sumF_withLeg)
Leg_plot <- ggdraw(Leg) #+ draw_label("Landmarks absent", x = 1, y = 0.95,
#            vjust = 1, hjust = 1, size = 25); Leg_plot

## FEMALES IN SUMMER without legend
sumF_plot <- ggplot(sumF) +
  
  # Resident
  geom_line(aes(x=TL, y=Resident_mle, colour="Resident")) + 
  geom_ribbon(aes(x=TL, ymin=Resident_low, ymax=Resident_upp, 
                  fill="Resident"), alpha=alpha.trans) +
  # Transient
  geom_line(aes(x=TL, y=Transient_mle, colour="Transient")) +
  geom_ribbon(aes(x=TL, ymin=Transient_low, ymax=Transient_upp, 
                  fill="Transient"), alpha=alpha.trans) +
  
  # Legends
  scale_colour_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                      labels=c("Resident","Transient")) + 
  scale_fill_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                    labels=c("Resident","Transient")) +
  ylim(0,ymax) + xlim(xxmin,xmax) +
  xlab("") + ylab("") + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## MALES IN SUMMER 
sumM_plot <- ggplot(sumM) +
  
  # Resident
  geom_line(aes(x=TL, y=Resident_mle, colour="Resident")) + 
  geom_ribbon(aes(x=TL, ymin=Resident_low, ymax=Resident_upp, 
                  fill="Resident"), alpha=alpha.trans) +
  # Transient
  geom_line(aes(x=TL, y=Transient_mle, colour="Transient")) +
  geom_ribbon(aes(x=TL, ymin=Transient_low, ymax=Transient_upp, 
                  fill="Transient"), alpha=alpha.trans) +
  
  # Legends
  scale_colour_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                      labels=c("Resident","Transient")) + 
  scale_fill_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                    labels=c("Resident","Transient")) +
  ylim(0,ymax) + xlim(xxmin,xmax) +
  xlab("") + ylab("") + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## FEMALES IN WINTER 
winF_plot <- ggplot(winF) +
  
  # Resident
  geom_line(aes(x=TL, y=Resident_mle, colour="Resident")) + 
  geom_ribbon(aes(x=TL, ymin=Resident_low, ymax=Resident_upp, 
                  fill="Resident"), alpha=alpha.trans) +
  # Transient
  geom_line(aes(x=TL, y=Transient_mle, colour="Transient")) +
  geom_ribbon(aes(x=TL, ymin=Transient_low, ymax=Transient_upp, 
                  fill="Transient"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                      labels=c("Resident","Transient")) + 
  scale_fill_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                    labels=c("Resident","Transient")) +
  ylim(0,1) + xlim(225,510) +
  xlab("") + ylab("") + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## MALES IN WINTER 
winM_plot <- ggplot(winM) +

  # Resident
  geom_line(aes(x=TL, y=Resident_mle, colour="Resident")) + 
  geom_ribbon(aes(x=TL, ymin=Resident_low, ymax=Resident_upp, 
                  fill="Resident"), alpha=alpha.trans) +
  # Transient
  geom_line(aes(x=TL, y=Transient_mle, colour="Transient")) +
  geom_ribbon(aes(x=TL, ymin=Transient_low, ymax=Transient_upp, 
                  fill="Transient"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                      labels=c("Resident","Transient")) + 
  scale_fill_manual(name="States", values = c("Resident" = mycols[1], "Transient" = mycols[2]),
                    labels=c("Resident","Transient")) +
  ylim(0,1) + xlim(225,510) +
  xlab("") + ylab("") + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# create season labels
Sum_lab <- ggdraw() + draw_label("Summer", x = 0.4, y = 0.45, angle = -90,
                                 vjust = 4, hjust = 1, size = 15); Sum_lab
Win_lab <- ggdraw() + draw_label("Winter", x = 0.4, y = 0.5, angle = -90,
                                 vjust = 4, hjust = 1, size = 15); Win_lab

# create common x and y labels
y.lab <- textGrob("Stationary state probabilities", 
                  vjust = 1.5, hjust = 0.4,
                  gp=gpar(fontsize=18), rot=90)

x.lab <- textGrob("Length (cm)", 
                  vjust = -10, hjust = 0.7,
                  gp=gpar(fontsize=18))

# create title plots
Title_plotSum <- ggdraw() + draw_text("Females", x = 0.65, y = 0.3,
                                      vjust = 3.2, hjust = 1, size = 15); Title_plotSum
Title_plotWin <- ggdraw() + draw_text("Males", x = 0.65, y = 0.3,
                                      vjust = 3.2, hjust = 1, size = 15); Title_plotWin

# add labels to plot
lay <- rbind(c(1,1,1,2,2,2,3),
             c(4,4,4,5,5,5,6),
             c(4,4,4,5,5,5,6),
             c(7,7,7,8,8,8,9),
             c(7,7,7,8,8,8,9),
             c(10,10,10,10,10,10,11))

TL_comp_plot <- grid.arrange(Title_plotSum, Title_plotWin, nullGrob(),
                             sumF_plot, sumM_plot, Sum_lab,
                             winF_plot, winM_plot, Win_lab,
                             Leg_plot, nullGrob(),
                             layout_matrix=lay)

#quartz()
TL_comp_plotL <- grid.arrange(arrangeGrob(TL_comp_plot, left = y.lab, bottom = x.lab))                         

ggsave(filename=here::here("figures","sex_vs_length_stationary.jpg"), 
       plot=TL_comp_plotL,
       width=25, height=25, units="cm",dpi=700)


#
# END OF SCRIPT
#
