
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
library(rerddapXtracto)
library(cmocean)
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

set.seed(2022)
tic()
fit_hmm_3states <- fitHMM(data = dat2,
                          nbStates = 3,
                          dist = list(step = "gamma", angle = "wrpcauchy"),
                          Par0 = list(step = stepPar0, angle = anglePar0),
                          formula = ~ 1,
                          estAngleMean = list(angle=TRUE),
                          stateNames = c('Encamped', 'ARS', 'Migratory'),
                          retryFits = 30)
toc()  #took 19.5 min to run

fit_hmm_3states

plot(fit_hmm_3states)  #no noticeable impact of including SST, although slight effect of yday; stationary states show some patterns
plotStates(fit_hmm_3states)  #given how quickly the states fluctuate between Encamped and ARS, there likely is only 1 slow state present
timeInStates(fit_hmm_3states)
plotPR(fit_hmm_3states, ncores = 5)  #looks about the same as before





### Fit HMM for 2 states including effect of daily SST and yday ###

## Extract sea surface temperature

#Available at url: https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.html

# Extract monthly SST
xpos <- dat2$lon
ypos <- dat2$lat
tpos <- dat2$date
sstInfo <- rerddap::info('jplMURSST41')

tic()
sst <- rxtracto(sstInfo, parameter = 'analysed_sst',
                xcoord = xpos, ycoord = ypos, tcoord = tpos,
                xlen = 0.2, ylen = 0.2, progress_bar = TRUE)
toc()  #takes 16 min to run

plotTrack(sst, xpos, ypos, tpos, plotColor = 'thermal')

# Add SST and year-day (yday) to dataset

dat2$sst <- sst$`mean analysed_sst`
dat2$yday <- yday(dat2$date)




# formula for transition probabilities
formula <- ~ sst * cosinor(yday, period = 365)

# initial parameters (obtained from nested model m1)
Par0_covar <- getPar0(model = fit_hmm_2states, formula = formula)



tic()
fit_hmm_2states_covar1 <- fitHMM(data = dat2,
                                nbStates = 2,
                                dist = list(step = "gamma", angle = "wrpcauchy"),
                                Par0 = Par0_covar$Par,
                                beta0 = Par0_covar$beta,
                                formula = formula,
                                estAngleMean = list(angle=TRUE),
                                stateNames = c('Encamped', 'Migratory'),
                                retryFits = 30)
toc()  #took 17 min to run

fit_hmm_2states_covar1
plot(fit_hmm_2states_covar1, plotCI = TRUE, plotStationary = TRUE)  #no noticeable impact of including SST, although slight effect of yday; stationary states show some patterns
plotStates(fit_hmm_2states_covar1)
stationary(fit_hmm_2states_covar1)  #gets stationary state probs plotted above per obs
# plotStationary(fit_hmm_2states_covar1, plotCI = TRUE)
timeInStates(fit_hmm_2states_covar1)
plotPR(fit_hmm_2states_covar1, ncores = 5)  #looks about the same as before








### Also include effect of covars on step lengths and turning angles

# formulas for parameters of state-dependent distributions
DM <- list(step = list(mean = formula,
                       sd = formula),
           angle = list(mean = formula,
                        concentration = formula))

# initial parameters (obtained from nested model; i.e., previous model)
Par0_covar2 <- getPar0(model = fit_hmm_2states_covar1, formula = formula, DM = DM)


tic()
fit_hmm_2states_covar2 <- fitHMM(data = dat2,
                                 nbStates = 2,
                                 dist = list(step = "gamma", angle = "wrpcauchy"),
                                 Par0 = Par0_covar2$Par,
                                 beta0 = Par0_covar2$beta,
                                 formula = formula,
                                 DM = DM,
                                 estAngleMean = list(angle=TRUE),
                                 stateNames = c('Encamped', 'Migratory'),
                                 retryFits = 10)
toc()  #took 12.5 min to run

fit_hmm_2states_covar2
plot(fit_hmm_2states_covar2, plotCI = TRUE, plotStationary = TRUE)  #affects of covariates are greater on SL and TA vs TPM; stationary states show some patterns
plotStates(fit_hmm_2states_covar2)
stationary(fit_hmm_2states_covar2)  #gets stationary state probs plotted above per obs
# plotStationary(fit_hmm_2states_covar2, plotCI = TRUE)
timeInStates(fit_hmm_2states_covar2)
plotPR(fit_hmm_2states_covar2, ncores = 5)  #looks about the same as before




### Compare among models ###

AIC(fit_hmm_2states, fit_hmm_3states, fit_hmm_2states_covar1, fit_hmm_2states_covar2)
AICweights(fit_hmm_2states, fit_hmm_3states, fit_hmm_2states_covar1, fit_hmm_2states_covar2)






### Fit HMM for 3 states that evaluates SL, TA, and NSD

# Calculate net-squared displacement manually
# calc_disp <- function(data, x, y) {
#   data$disp <- sqrt((data[,"x"] - data[,"x"][1])^2 + (data[,"y"] - data[,"y"][1])^2)
#
#   return(data)
# }
#
# dat3 <- dat2 %>%
#   split(.$ID) %>%
#   purrr::map(., calc_disp, x, y) %>%
#   bind_rows()
#
#
# foo <- dat3 %>%
#   split(.$ID) %>%
#   purrr::map(., ~{.x %>%
#       mutate(disp_diff = c(diff(disp), NA))
#   }) %>%
#   bind_rows()
#
# foo <- foo %>%
#   mutate(state = case_when(disp < 6 ~ 'res_breed',
#                            disp_diff > 2 ~ 'mig',
#                            disp > 6 & disp_diff < 2 ~ 'res_forag'))
#
#
# ggplot(foo, aes(date, disp)) +
#   geom_path(aes(group = ID, color = state)) +
#   theme_bw() +
#   facet_wrap(~ID, scales = "free_x")
#
#
# ggplot(foo, aes(disp, fill = state)) +
#   geom_histogram(alpha = 0.6) +
#   theme_bw()
#
# ggplot(foo, aes(step, fill = state)) +
#   geom_histogram(alpha = 0.6) +
#   theme_bw()
#
# ggplot(foo, aes(angle, fill = state)) +
#   geom_histogram(alpha = 0.6) +
#   theme_bw()
#
#
# foo %>%
#   group_by(state) %>%
#   summarize(mean.step = mean(step, na.rm = T),
#             sd.step = sd(step, na.rm = T),
#             mean.disp = mean(disp, na.rm = T),
#             sd.disp = sd(disp, na.rm = T))
#
#
# # Step lengths
# ggplot(dat3, aes(date, step)) +
#   geom_line() +
#   theme_bw() +
#   facet_wrap(~ID, scales = "free_x")
# #try means at 1, 1, and 20 km; SDs at 0.5, 1, and 5 km
#
# # Turning angles
# ggplot(dat3, aes(date, angle)) +
#   geom_line() +
#   theme_bw() +
#   facet_wrap(~ID, scales = "free_x")
# #try means at 0 rad; concentrations at 0.1, 0.1, 0.99
#
# # Displacement
# ggplot(dat3, aes(disp)) +
#   geom_histogram(binwidth = 50) +
#   theme_bw()
#
# ggplot(dat3, aes(date, disp)) +
#   geom_line() +
#   theme_bw() +
#   facet_wrap(~ID, scales = "free_x")
# #try means at 0, 500, and 250 km; SDs at 2, 200, and 100 km
#
#
#
#
# # initial step distribution natural scale parameters
# stepPar0 <- c(0.5, 1.5, 15, 1, 1.5, 8) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3)
#
# # initial angle distribution natural scale parameters
# anglePar0 <- c(0, 0, 0, 0.1, 0.1, 0.99) # (mean_1, mean_2, mean_3, concentration_1, concentration_2, concentration_3)
#
# #initial displacement distribution natural scale parameters
# whichzero <- which(dat3$disp == 0)
# propzero <- length(whichzero)/nrow(dat3)
# zeromass0 <- c(propzero, 1e-9, 1e-9)        #for zero distances by state
# dispPar0 <- c(1, 570, 450, 1, 213, 250, zeromass0) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3, proportion of zeroes likely present per state)
#
#
#
# set.seed(2022)
# tic()
# fit_hmm_3states_3covars <- fitHMM(data = dat3,
#                                   nbStates = 3,
#                                   dist = list(step = "gamma", angle = "wrpcauchy", disp = "gamma"),
#                                   Par0 = list(step = stepPar0, angle = anglePar0, disp = dispPar0),
#                                   formula = ~ 1,
#                                   estAngleMean = list(angle=TRUE),
#                                   stateNames = c('Resident_Breeding','Resident_Foraging','Migratory'))
# toc()  #took 25 min to run
#
# fit_hmm_3states_3covars
#
# plot(fit_hmm_3states_3covars)
# plotStates(fit_hmm_3states_3covars)  #given how quickly the states fluctuate between Encamped and ARS, there likely is only 1 slow state present
# timeInStates(fit_hmm_3states_3covars)
# plotPR(fit_hmm_3states_3covars, ncores = 5)  #looks about the same as before



