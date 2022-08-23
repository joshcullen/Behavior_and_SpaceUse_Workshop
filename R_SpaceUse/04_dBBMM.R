
########################################################
### Calculate dynamic Brownian Bridge Movement Model ###
########################################################

library(tidyverse)
library(lubridate)
library(move)  #v4.0.6
library(sf)  #v1.0.7
library(rnaturalearth)
library(plotly)
library(tictoc)
library(MetBrewer)




### Load data ###

dat <- read.csv('Processed_data/SSM_mp6hr_FDN Cmydas tracks.csv')

glimpse(dat)
summary(dat)




### Wrangle and prep data for dBBMM estimation ###

dat <- dat %>%
  mutate(date = as_datetime(date),  #convert to datetime format
         SE = ifelse(x.se > y.se, y.se, x.se))  #create new column with smallest location SE

# Viz time steps per ID
# dat <- dat %>%
#   split(.$id) %>%
#   purrr::map(., ~mutate(.x,
#                         dt = difftime(c(date[-1], NA),
#                                       date,
#                                       units = "hours") %>%  #calculate time step per ID
#                           as.numeric())
#   ) %>%
#   bind_rows()
#
#
#   ggplot(dat) +
#   geom_point(aes(date, dt)) +
#   geom_hline(aes(yintercept = median(dt, na.rm = T)), color = "red", size = 1) +
#   theme_bw() +
#   facet_wrap(~ id, scales = "free")
#
# median(dat$dt, na.rm = TRUE)  #0.88 hrs
# dat %>%
#   group_by(id) %>%
#   summarize(median.dt = median(dt, na.rm = TRUE))
# let's round and say that typical time step across all IDs is 1 hr


# Viz histogram of smallest location SEs
hist(dat$SE)  #units are in km; vast majority (99.6%) under 2 km


# Split data into list by ID
dat2 <- dat %>%
  split(.$id)






### Run dBBMM model ###

# Create 'move' object
dat.list <- vector("list", length(dat2))  #to store dBBMM results
contours<- vector("list", length(dat2))  #to store resulting 50% and 95% UD contours


# Estimate separately by ID
for (i in 1:length(dat.list)) {
  print(paste("ID", dat2[[i]]$id[1]))  #print current ID

  dat.mov <- move(x = dat2[[i]]$x, y = dat2[[i]]$y, time = dat2[[i]]$date, data = dat2[[i]],
                  proj = CRS("+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs"),
                  animal = dat2[[i]]$id)

  # Conditionally define extent; necessary for turtles that don't migrate
  x.ext <- diff(dat.mov@bbox[1,])
  rast.ext <- ifelse(x.ext < 100, 3, 0.3)


  ## Remove variance from segments w/ large time gaps to prevent errors when running dBBMM model

  # calculate the dynamic brownian motion variance of the gappy track
  # tic()
  # dbbv <- brownian.motion.variance.dyn(dat.mov, location.error = dat2[[i]]$SE, window.size = (4*7)+1,
  #                                      margin = 9)
  # toc()

  # ignore all variance values when time interval is > 24 hrs
  # dbbv@interest[unlist(timeLag(dat.mov, "hours")) > 24] <- FALSE



  ## Run dBBMM (this will take a little while to run)

  tic()
  dat.list[[i]] <- brownian.bridge.dyn(object = dat.mov, raster = 0.25, location.error = "SE",
                                       margin = 9, window.size = 29, ext = rast.ext)
  toc()


  ## Extract 50 and 95% contours of space-use
  res <- raster2contour(dat.list[[i]], levels = c(0.5, 0.95))

  contours[[i]] <- st_as_sf(res)

  if (st_geometry_type(contours[[i]])[1] == 'LINESTRING') {
    contours[[i]] <- st_cast(contours[[i]], "POLYGON")
  } else {
    contours[[i]] <- st_cast(contours[[i]], "MULTIPOLYGON") %>%
      st_make_valid()  #fixes issue w/ negative areas being calculated
  }

}


# Quick plot of UDs at defined kernel volume (in raster form)
ud95 <- getVolumeUD(dat.list[[1]])
ud95[ud95 > 0.95] <- NA
plot(ud95, main="UD95")

ud50 <- getVolumeUD(dat.list[[1]])
ud50[ud50 > 0.50] <- NA
plot(ud50, main="UD50")


names(contours) <- names(dat2)
map(contours, st_area)
contours <- map2(.x = contours, .y = names(contours),
                 .f = ~{.x %>%
                     mutate(id = .y)})  #add ID name

contours2 <- do.call(rbind, contours)





## Let's viz the results

ggplot() +
  geom_sf(data = contours2, aes(color = level), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  facet_wrap(~ level, nrow = 2)

ggplot() +
  geom_sf(data = contours2, aes(color = level), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  facet_wrap(~ id)


# Highlight a few of the different IDs

#205540
ggplot() +
  geom_path(data = dat2[['205540']], aes(x, y), size = 0.5, alpha = 0.5) +
  geom_sf(data = contours2 %>%
            filter(id == 205540), aes(color = level), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  labs(title = 'ID 205540') +
  theme_bw()

#205542
ggplot() +
  geom_path(data = dat2[['205542']], aes(x, y), size = 0.5, alpha = 0.25) +
  geom_sf(data = contours2 %>%
            filter(id == 205542), aes(color = level), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  labs(title = 'ID 205542') +
  theme_bw()

#205537
ggplot() +
  geom_path(data = dat2[['205537']], aes(x, y), size = 0.5, alpha = 0.5) +
  geom_sf(data = contours2 %>%
            filter(id == 205537), aes(color = level), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  labs(title = 'ID 205537') +
  theme_bw()

#226073
ggplot() +
  geom_path(data = dat2[['226073']], aes(x, y), size = 0.5, alpha = 0.5) +
  geom_sf(data = contours2 %>%
            filter(id == 226073), aes(color = level), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  labs(title = 'ID 226073') +
  theme_bw()




# Plot all tracks w/ UD contours
ggplot() +
  geom_path(data = dat, aes(x, y, color = factor(id), group = id), size = 0.5, alpha = 0.5) +
  geom_sf(data = contours2, fill = "transparent", size = 0.75, color = 'black') +
  scale_color_viridis_d(option = 'turbo') +
  theme_bw() +
  facet_wrap(~ level, nrow = 2)
