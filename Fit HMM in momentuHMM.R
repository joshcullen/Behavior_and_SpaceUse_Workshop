
#############################################
### Fit discrete-time hidden Markov model ###
#############################################

### Additional notes to discuss:
# * {momentuHMM} can estimate multiple discrete behavioral states using any number of covariates; this includes ancillary biologging data and evironmental variables
# * Options for BRW and CRW models, including activity centers that have attractive or repulsive forces
# * Option to perform multiple imputation via wrapper function for {crawl} methods, which can provide options to analyze tracks w/ measurement error and/or irregular time series indirectly
# * Can include random effects by ID on TPM
# * Can accomodate data streams/behaviors at multiple time scales (hierarchical HMM)

library(tidyverse)
library(lubridate)
library(momentuHMM)  #v1.5.4
library(sf)  #v1.0.7
library(rnaturalearth)
library(tictoc)
library(plotly)
# library(future)  #needed to properly run foieGras::osar() in parallel


### Load data ###

dat <- read.csv('Processed_data/SSM_mp6hr_FDN Cmydas tracks.csv')

glimpse(dat)
summary(dat)




### Wrangle data for analysis using {momentuHMM} ###

# Convert all 'Date' to datetime format
dat <- dat %>%
  mutate(date = as_datetime(date)) %>%
  rename(ID = id)

dat2 <- prepData(dat, type = 'UTM', coordNames = c('x','y'))  #can also be done using lat/long
plot(dat2)


# Plot time series and distributions of step lengths
ggplot(dat2, aes(step)) +
  geom_histogram(binwidth = 0.5) +
  theme_bw()

ggplot(dat2, aes(date, step)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")
#for 2 states, initial mean values probably around 1 and 15 km whereas SD probably around 0.5 and 5


# Plot time series and distributions of turning angles
ggplot(dat2, aes(angle)) +
  geom_histogram(binwidth = pi/8) +
  theme_bw()

ggplot(dat2, aes(date, angle)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")
#for 2 states, initial mean values probably both around 0 rad whereas concentration probably around 0.99 and 0.2





### Fit HMM for 2 states using step lengths and turning angles ###

# initial step distribution natural scale parameters
stepPar0 <- c(1, 15, 0.5, 5) # (mu_1, mu_2, sd_1, sd_2)

# initial angle distribution natural scale parameters
anglePar0 <- c(0, 0, 0.2, 0.99) # (mean_1, mean_2, concentration_1, concentration_2)


tic()
fit_hmm_2states <- fitHMM(data = dat2,
                          nbStates = 2,
                          dist = list(step = "gamma", angle = "wrpcauchy"),
                          Par0 = list(step = stepPar0, angle = anglePar0),
                          formula = ~ 1,
                          estAngleMean = list(angle=TRUE),
                          stateNames = c('Encamped', 'Migratory'),
                          retryFits = 30)
toc()  #took 4.5 min to run

fit_hmm_2states
plot(fit_hmm_2states)
plotPR(fit_hmm_2states, ncores = 5)  #plot of pseudo-residuals show that there are likely some problems





### Fit HMM for 3 states using step lengths and turning angles ###

# Step lengths
ggplot(dat2, aes(date, step)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")
#try means at 1, 10, and 20 km; SDs at 0.5, 2, and 5 km

# Turning angles
ggplot(dat2, aes(date, angle)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")
#try means at 0 rad; concentrations at 0.1, 0.5, 0.99

# initial step distribution natural scale parameters
stepPar0 <- c(1, 10, 20, 0.5, 2, 5) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3)

# initial angle distribution natural scale parameters
anglePar0 <- c(0, 0, 0, 0.1, 0.5, 0.99) # (mean_1, mean_2, mean_3, concentration_1, concentration_2, concentration_3)


tic()
fit_hmm_3states <- fitHMM(data = dat2,
                          nbStates = 3,
                          dist = list(step = "gamma", angle = "wrpcauchy"),
                          Par0 = list(step = stepPar0, angle = anglePar0),
                          formula = ~ 1,
                          estAngleMean = list(angle=TRUE),
                          stateNames = c('Encamped', 'ARS', 'Migratory'),
                          retryFits = 30)
toc()  #took 19 min to run

fit_hmm_3states
plot(fit_hmm_3states)
plotPR(fit_hmm_3states, ncores = 5)  #looks about the same as before, but QQ plots look better





### Fit HMM for 2 states including effect of daily SST and yday ###

