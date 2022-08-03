
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
library(future)  #needed to properly run foieGras::osar() in parallel


### Load data ###

dat <- read.csv('Processed_data/Cleaned_FDN Cmydas tracks.csv')

glimpse(dat)
summary(dat)


### Wrangle data for analysis using {foieGras} ###

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



### Inspect time steps of transmissions for making predictions ###

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
# it looks like there are some very large gaps in the tracks for PTTs 41587, 41588; let's remove these points that occur before (or after) depending on the ID


# Find rows where 1st dt of an ID is > 7 d (604800 secs)
ind <- tmp %>%
  group_by(id) %>%
  filter(dt > 3600*24*7) %>%
  slice(1)

# Remove obs before or after a gap of 7 days
dat3 <- dat2 %>%
  filter(!(id == 41587 & date <= ind$date[1])) %>%
  filter(!(id == 41588 & date >= ind$date[2]))

# Determine primary time step
ggplot(tmp) +
  geom_histogram(aes(dt), binwidth = 3600) +
  theme_bw() +
  xlim(0,3600*24)

table(tmp$dt) %>%
  sort(decreasing = TRUE)

tmp %>%
  filter(!(id == 41587 & date <= ind$date[1])) %>%
  filter(!(id == 41588 & date >= ind$date[2])) %>%
  group_by(id) %>%
  summarize(mean = mean(dt, na.rm = TRUE),
            median = median(dt, na.rm = TRUE))
# Mean/median time step is ~1 hr





### Run SSM ###

# Change `id` to character to avoid problems during model runs
dat3$id <- as.character(dat3$id)



### Using correlated random walk (CRW) model

# Estimate 'true' locations on irregular sampling interval (by setting `time.step = NA`)
tic()
fit_crw_fitted <- fit_ssm(dat3, vmax = 3, model = "crw", time.step = NA,
                          control = ssm_control(verbose = 1))
toc()  #took 2 min where time.step = NA

print(fit_crw_fitted)


# Viz the filtered outliers (gold), raw observations (blue), and estimated locations (red), along with associated uncertainty (red shading)
plot(fit_crw_fitted, what = "fitted", type = 1, ask = TRUE)
plot(fit_crw_fitted, what = "fitted", type = 2, ask = TRUE)

map(fit_crw_fitted, what = "fitted", aes = aes_lst(mp_pal = hcl.colors(n=100, "viridis")))


# Estimate behavioral state (i.e., move persistence; gamma)
# Joint move persistence model ('jmpm') uses hierarchical approach across IDs
tic()
fit_crw_jmpm_fitted <- fit_mpm(fit_crw_fitted, what = "fitted", model = "jmpm",
                              control = mpm_control(verbose = 1))
toc()  #took 4.5 min to fit

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





# Estimate 'true' locations on regular sampling interval of 6 hrs
tic()
fit_crw_6hr <- fit_ssm(dat3, vmax = 3, model = "crw", time.step = 6,
                          control = ssm_control(verbose = 1))
toc()  #took 2 min to fit model

print(fit_crw_6hr)


# Viz the filtered outliers (gold), raw observations (blue), and estimated locations (red), along with associated uncertainty (red shading)
plot(fit_crw_6hr, what = "predicted", type = 1, ask = TRUE)
plot(fit_crw_6hr, what = "predicted", type = 2, ask = TRUE)

# plot(fit_crw_6hr, what = "fitted", type = 1, ask = TRUE)
# plot(fit_crw_6hr, what = "fitted", type = 2, ask = TRUE)




# Estimate behavioral state (i.e., move persistence; gamma)
# Individual move persistence model ('mpm') estimates behavioral states separately across IDs
tic()
fit_crw_mpm_6hr <- fit_mpm(fit_crw_6hr, what = "predicted", model = "mpm",
                              control = mpm_control(verbose = 1))
toc()  #took 49 sec to fit

print(fit_crw_mpm_6hr)  #model for PTT 205537 didn't converge
plot(fit_crw_mpm_6hr)



# Grab results and plot
res_crw_6hr<- join(ssm = fit_crw_6hr,
                      mpm = fit_crw_mpm_6hr,
                      what.ssm = "predicted")


# Compare raw tracks vs fitted tracks (for adults tagged at Fernando de Noronha)
ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat3, aes(lon, lat, group = id), color = 'black') +  #raw tracks
  geom_path(data = res_crw_6hr, aes(lon, lat, group = id), color = "blue") +  #modeled tracks
  theme_bw() +
  facet_wrap(~id) +
  coord_sf(xlim = c(-42, -32), ylim = c(-8, -1))

# Viz modeled tracks together
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_path(data = res_crw_6hr, aes(lon, lat, group = id, color = id), size = 0.75, alpha = 0.8) +
    scale_color_viridis_d() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)

# Viz modeled tracks together w/ behavior plotted
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_point(data = res_crw_6hr, aes(lon, lat, group = id, color = g), size = 0.75, alpha = 0.8) +
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
    geom_path(data = res_crw_6hr, aes(lon, lat, group = id, color = id), size = 0.75, alpha = 0.8) +
    scale_color_viridis_d() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)




# Estimate 'true' locations on regular sampling interval of 6 hrs using 'move persistence' model
tic()
fit_mp_6hr <- fit_ssm(dat3, vmax = 3, model = "mp", time.step = 6,
                       control = ssm_control(verbose = 1))
toc()  #took 5 min to fit model

print(fit_mp_6hr)


# Viz the filtered outliers (gold), raw observations (blue), and estimated locations (red), along with associated uncertainty (red shading)
plot(fit_mp_6hr, what = "predicted", type = 1, ask = TRUE)
plot(fit_mp_6hr, what = "predicted", type = 2, ask = TRUE)
plot(fit_mp_6hr, what = "predicted", type = 3, ask = TRUE)

# plot(fit_mp_6hr, what = "fitted", type = 1, ask = TRUE)
# plot(fit_mp_6hr, what = "fitted", type = 2, ask = TRUE)
# plot(fit_mp_6hr, what = "fitted", type = 3, ask = TRUE)




# Grab results and plot
res_mp_6hr<- grab(fit_mp_6hr, what = "predicted")


# Viz modeled tracks together
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_path(data = res_mp_6hr, aes(lon, lat, group = id, color = id), size = 0.75, alpha = 0.8) +
    scale_color_viridis_d() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)

# Viz modeled tracks together w/ behavior plotted
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_point(data = res_mp_6hr, aes(lon, lat, group = id, color = g), size = 0.75, alpha = 0.8) +
    scale_color_viridis_c(option = "inferno") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)



### Compare estimates across the different approaches ###

# Behavioral states
ggplot() +
  geom_line(data = res_crw_fitted, aes(date, g, color = "CRW_Irregular")) +
  geom_line(data = res_crw_6hr, aes(date, g, color = "CRW_6hr")) +
  geom_line(data = res_mp_6hr, aes(date, g, color = "MP_6hr")) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  theme_bw() +
  facet_wrap(~id, scales = "free_x")

## Of these 3 different models, the CRW model w/ 6 hr prediction time step seems to give best broad scale state estimates


# Track relocations
ggplot() +
  geom_path(data = res_crw_fitted, aes(x, y, color = "CRW_Irregular")) +
  geom_path(data = res_crw_6hr, aes(x, y, color = "CRW_6hr")) +
  geom_path(data = res_mp_6hr, aes(x, y, color = "MP_6hr")) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  theme_bw() +
  facet_wrap(~id, scales = "free")

## Of these 3 approaches, the move persistence model seems to provide the best estimates when the individual appears to be stationary  (i.e., fewer unnecessary looping movements)


### Overall verdict: for relatively fine-scale movements (time step > 1 hr & < 1 d), the CRW model seems to provide the best behavioral state estimates. Whlie all 3 models provided very similar estimates for the 'true' animal movements, the MP model seems to provide the best estimates, particularly when the animal appears to be encamped within a given site


### Since the other behavioral state models do not account for location error and the analysis of step lengths and turning angles (if used) must be at a regular time interval, we will export the tracks from the move persistence model for further use.




### Export fitted tracks ###

write.csv(res_mp_6hr, "Processed_data/SSM_mp6hr_FDN Cmydas tracks.csv", row.names = FALSE)

