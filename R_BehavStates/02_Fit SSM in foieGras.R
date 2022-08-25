
#################################################################################
### Fit continuous-time state-space model while accounting for location error ###
#################################################################################

library(tidyverse)
library(lubridate)
library(foieGras)  #v1.0-7
library(sf)  #v1.0.7
library(rnaturalearth)
library(tictoc)
library(plotly)
# library(future)  #needed to properly run foieGras::osar() in parallel


#### Load data ####

dat <- read.csv('Processed_data/Cleaned_FDN Cmydas tracks.csv')

glimpse(dat); str(dat)
summary(dat)


#### Wrangle data for analysis using {foieGras} ####

# Convert all 'Quality' values to "G" for FastGPS data and 'Date' to datetime format
dat <- dat %>%
  mutate(Quality = ifelse(Type == 'FastGPS', 'G', Quality),
         Date = as_datetime(Date))


# Rename columns for {foieGras}
dat2<- dat %>%
  rename(id = Ptt, date = Date, lc = Quality, lon = Longitude, lat = Latitude,
         eor = Error.Ellipse.orientation, smaj = Error.Semi.major.axis,
         smin = Error.Semi.minor.axis) %>%
  dplyr::select(id, date, lc, lon, lat, smaj, smin, eor)  #reorders and subsets the columns

glimpse(dat2)



#### Inspect time steps of transmissions for making predictions ####

tmp <- dat2 %>%
  split(.$id) %>%
  purrr::map(., ~mutate(.x,
                        dt = difftime(c(date[-1], NA),
                                      date,
                                      units = "secs") %>%
                          as.numeric())
             ) %>%
  bind_rows()


ggplot(tmp, aes(date, dt)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~id, scales = "free_x")
# some outliers present, but nothing to worry about for now



# Determine primary time step
ggplot(tmp) +
  geom_histogram(aes(dt), binwidth = 3600) +
  theme_bw() +
  xlim(0, 3600*24)


tmp %>%
  group_by(id) %>%
  summarize(mean = mean(dt, na.rm = TRUE),
            median = median(dt, na.rm = TRUE))
# Mean/median time step is ~1 hr




#################
#### Run SSM ####
#################

# Change `id` to character to avoid problems during model runs
dat2$id <- as.character(dat2$id)



#### Account for location error at observed irregular time interval ####

# Estimate 'true' locations on irregular sampling interval (by setting `time.step = NA`)
tic()
fit_crw_fitted <- fit_ssm(dat2, vmax = 3, model = "crw", time.step = NA,
                          control = ssm_control(verbose = 1))
toc()  #took 30 sec where time.step = NA

print(fit_crw_fitted)


# Viz the filtered outliers (gold), raw observations (blue), and estimated locations (red), along with associated uncertainty (red shading)
plot(fit_crw_fitted, what = "fitted", type = 1, ask = TRUE)
plot(fit_crw_fitted, what = "fitted", type = 2, ask = TRUE)

foieGras::map(fit_crw_fitted,
    what = "fitted",
    by.id = TRUE)


# Estimate behavioral state (i.e., move persistence; gamma)
# Joint move persistence model ('jmpm') uses hierarchical approach across IDs
tic()
fit_crw_jmpm_fitted <- fit_mpm(fit_crw_fitted, what = "fitted", model = "jmpm",
                              control = mpm_control(verbose = 1))
toc()  #took 48 sec to fit

print(fit_crw_jmpm_fitted)
plot(fit_crw_jmpm_fitted)


# Grab results and plot
res_crw_fitted<- join(ssm = fit_crw_fitted,
                      mpm = fit_crw_jmpm_fitted,
                      what.ssm = "fitted")


# Compare raw tracks vs fitted tracks (for adults tagged at Fernando de Noronha)
brazil<- ne_countries(scale = 50, country = "Brazil", returnclass = 'sf')

ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat3, aes(lon, lat, group = id), color = 'black') +  #raw tracks
  geom_path(data = res_crw_fitted, aes(lon, lat, group = id), color = "blue") +  #modeled tracks
  theme_bw() +
  facet_wrap(~id) +
  coord_sf(xlim = c(-42, -32), ylim = c(-8, -1))


# Viz modeled tracks together
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_path(data = res_crw_fitted, aes(lon, lat, group = id, color = id), size = 0.75, alpha = 0.8) +
    scale_color_viridis_d() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)

# Viz modeled tracks together w/ behavior plotted
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_point(data = res_crw_fitted, aes(lon, lat, group = id, color = g), size = 0.75, alpha = 0.8) +
    scale_color_viridis_c(option = "inferno") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)


#interactively explore using projected coordinates (World Mercator Projection; EPSG:3395, units = km)
res_crw_fitted %>%
  dplyr::select(id, date, x, y, s, g) %>%
  data.frame() %>%  #current version of shiny_tracks won't work w/ tibble format or when a column has all NAs
  bayesmove::shiny_tracks(., "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs")


#interactively explore using lat/long (EPSG:4326)
res_crw_fitted %>%
  dplyr::select(-c(x, y, s.se)) %>%
  rename(x = lon, y = lat) %>%
  dplyr::select(id, date, x, y, s, g) %>%
  data.frame() %>%  #current version of shiny_tracks won't work w/ tibble format or when a column has all NAs
  bayesmove::shiny_tracks(., 4326)


# Check model fit w/ diagnostic "one-step-ahead residuals"; takes long time!
# diag_crw <- osar(fit_crw_fitted)
#
# plot(diag_crw, type = "ts")  #time series of residuals
# plot(diag_crw, type = "qq")  #Q-Q plot
# plot(diag_crw, type = "acf")  #Autocorrelation function plot




#### Account for location error at regularized time interval with correlated random walk ####

# Estimate 'true' locations on regular sampling interval of 8 hrs
tic()
fit_crw_8hr <- fit_ssm(dat2, vmax = 3, model = "crw", time.step = 8,
                          control = ssm_control(verbose = 1))
toc()  #took 34 sec to fit model

print(fit_crw_8hr)


# Viz the filtered outliers (gold), raw observations (blue), and estimated locations (red), along with associated uncertainty (red shading)
plot(fit_crw_8hr, what = "predicted", type = 1, ask = TRUE)
plot(fit_crw_8hr, what = "predicted", type = 2, ask = TRUE)

# plot(fit_crw_8hr, what = "fitted", type = 1, ask = TRUE)
# plot(fit_crw_8hr, what = "fitted", type = 2, ask = TRUE)




# Estimate behavioral state (i.e., move persistence; gamma)
# Individual move persistence model ('mpm') estimates behavioral states separately across IDs
tic()
fit_crw_mpm_8hr <- fit_mpm(fit_crw_8hr, what = "predicted", model = "mpm",
                              control = mpm_control(verbose = 1))
toc()  #took 6 sec to fit
#if issues w/ model not converging, try changing time step of SSM and re-running

print(fit_crw_mpm_8hr)
plot(fit_crw_mpm_8hr)



# Grab results and plot
res_crw_8hr<- join(ssm = fit_crw_8hr,
                      mpm = fit_crw_mpm_8hr,
                      what.ssm = "predicted")


# Compare raw tracks vs fitted tracks (for adults tagged at Fernando de Noronha)
ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat3, aes(lon, lat, group = id), color = 'black') +  #raw tracks
  geom_path(data = res_crw_8hr, aes(lon, lat, group = id), color = "blue") +  #modeled tracks
  theme_bw() +
  facet_wrap(~id) +
  coord_sf(xlim = c(-42, -32), ylim = c(-8, -1))

# Viz modeled tracks together
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_path(data = res_crw_8hr, aes(lon, lat, group = id, color = id), size = 0.75, alpha = 0.8) +
    scale_color_viridis_d() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)

# Viz modeled tracks together w/ behavior plotted
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_point(data = res_crw_8hr, aes(lon, lat, group = id, color = g), size = 0.75, alpha = 0.8) +
    scale_color_viridis_c(option = "inferno") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)


# Compare predicted tracks against tracks fitted at irregular time step
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_path(data = res_crw_fitted, aes(lon, lat, group = id, color = id), size = 0.75, alpha = 0.3) +
    geom_path(data = res_crw_8hr, aes(lon, lat, group = id, color = id), size = 0.75, alpha = 0.8) +
    scale_color_viridis_d() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)





#### Account for location error at regularized time interval w/ move persistence model ####

# Estimate 'true' locations on regular sampling interval of 8 hrs using 'move persistence' model
tic()
fit_mp_8hr <- fit_ssm(dat2, vmax = 3, model = "mp", time.step = 8,
                       control = ssm_control(verbose = 1))
toc()  #took 1.25 min to fit model

print(fit_mp_8hr)


# Viz the filtered outliers (gold), raw observations (blue), and estimated locations (red), along with associated uncertainty (red shading)
plot(fit_mp_8hr, what = "predicted", type = 1, ask = TRUE)
plot(fit_mp_8hr, what = "predicted", type = 2, ask = TRUE)
plot(fit_mp_8hr, what = "predicted", type = 3, ask = TRUE)

# plot(fit_mp_8hr, what = "fitted", type = 1, ask = TRUE)
# plot(fit_mp_8hr, what = "fitted", type = 2, ask = TRUE)
# plot(fit_mp_8hr, what = "fitted", type = 3, ask = TRUE)




# Grab results and plot
res_mp_8hr_pred<- grab(fit_mp_8hr, what = "predicted")


# Viz modeled tracks together
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_path(data = res_mp_8hr_pred, aes(lon, lat, group = id, color = id), size = 0.75,
              alpha = 0.8) +
    scale_color_viridis_d() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)

# Viz modeled tracks together w/ behavior plotted
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_point(data = res_mp_8hr_pred, aes(lon, lat, group = id, color = g), size = 0.75,
               alpha = 0.8) +
    scale_color_viridis_c(option = "inferno") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)



### Compare estimates across the different approaches ###

# Behavioral states
ggplot() +
  geom_line(data = res_crw_fitted, aes(date, g, color = "CRW_Irregular")) +
  geom_line(data = res_crw_8hr, aes(date, g, color = "CRW_8hr")) +
  geom_line(data = res_mp_8hr_pred, aes(date, g, color = "MP_8hr")) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  theme_bw() +
  facet_wrap(~id, scales = "free_x")

## Of these 3 different models, the CRW model w/ 8 hr prediction time step seems to give best broad scale state estimates


# Track relocations
ggplot() +
  geom_path(data = res_crw_fitted, aes(x, y, color = "CRW_Irregular")) +
  geom_path(data = res_crw_8hr, aes(x, y, color = "CRW_8hr")) +
  geom_path(data = res_mp_8hr_pred, aes(x, y, color = "MP_8hr")) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  theme_bw() +
  facet_wrap(~id, scales = "free")

## Of these 3 approaches, the move persistence model seems to provide the best estimates when the individual appears to be stationary  (i.e., fewer unnecessary looping movements)


### Overall verdict: for relatively fine-scale movements (time step > 1 hr & < 1 d), the CRW model seems to provide the best behavioral state estimates. Whlie all 3 models provided very similar estimates for the 'true' animal movements, the MP model seems to provide the best estimates, particularly when the animal appears to be encamped within a given site


### Since the other behavioral state models do not account for location error and the analysis of step lengths and turning angles (if used) must be at a regular time interval, we will export the tracks from the move persistence model for further use.




### Export fitted tracks ###

write.csv(res_mp_8hr_pred, "Processed_data/SSM_mp8hr_FDN Cmydas tracks.csv", row.names = FALSE)

save(fit_crw_8hr, fit_crw_mpm_8hr, file = "Processed_data/SSM_model_fits.RData")
