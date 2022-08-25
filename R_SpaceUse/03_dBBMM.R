
################################################################
### Calculate dynamic Brownian Bridge Movement Model (dBBMM) ###
################################################################

library(tidyverse)
library(lubridate)
library(move)  #v4.0.6
library(sf)  #v1.0.7
library(rnaturalearth)
library(plotly)
library(tictoc)
library(MetBrewer)




#### Load data ####

dat <- read.csv('Processed_data/SSM_mp8hr_FDN Cmydas tracks.csv')

glimpse(dat)
summary(dat)




#### Wrangle and prep data for dBBMM estimation ####

dat <- dat %>%
  mutate(date = as_datetime(date),  #convert to datetime format
         SE = ifelse(x.se > y.se, y.se, x.se))  #create new column with smallest location SE



# Viz histogram of smallest location SEs
hist(dat$SE)  #units are in km; vast majority (97%) under 1 km


# Split data into list by ID
dat2 <- dat %>%
  split(.$id)






#### Run dBBMM model ####

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
plot(ud95, main = paste("ID", names(dat2)[[1]], ": 95% UD"))

ud50 <- getVolumeUD(dat.list[[1]])
ud50[ud50 > 0.50] <- NA
plot(ud50, main = paste("ID", names(dat2)[[1]], ": 50% UD"))


names(contours) <- names(dat2)
map(contours, st_area)
contours <- map2(.x = contours, .y = names(contours),
                 .f = ~{.x %>%
                     mutate(id = .y)})  #add ID name

contours2 <- do.call(rbind, contours) %>%
  rename(UD.level = level)





## Let's viz the results

brazil<- ne_countries(scale = 50, country = "Brazil", returnclass = 'sf') %>%
  st_transform(crs = "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs")

ggplot() +
  geom_sf(data = brazil) +
  geom_sf(data = contours2, aes(color = UD.level), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  facet_wrap(~ UD.level, nrow = 2) +
  coord_sf(xlim = c(min(dat$x) - 50, max(dat$x) + 50),
           ylim = c(min(dat$y) - 50, max(dat$y) + 50))

ggplot() +
  geom_sf(data = brazil) +
  geom_sf(data = contours2, aes(color = UD.level), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  facet_wrap(~ id) +
  coord_sf(xlim = c(min(dat$x) - 50, max(dat$x) + 50),
           ylim = c(min(dat$y) - 50, max(dat$y) + 50))





#### Highlight IDs individually ####

#205537
ggplot() +
  geom_path(data = dat2[['205537']], aes(x, y), size = 0.75, alpha = 0.5) +
  geom_sf(data = contours2 %>%
            filter(id == 205537), aes(color = UD.level), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  labs(title = 'ID 205537') +
  theme_bw()

#205540
ggplot() +
  geom_path(data = dat2[['205540']], aes(x, y), size = 0.75, alpha = 0.5) +
  geom_sf(data = contours2 %>%
            filter(id == 205540), aes(color = UD.level), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  labs(title = 'ID 205540') +
  theme_bw()

#205544
ggplot() +
  geom_path(data = dat2[['205544']], aes(x, y), size = 0.75, alpha = 0.25) +
  geom_sf(data = contours2 %>%
            filter(id == 205544), aes(color = UD.level), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  labs(title = 'ID 205544') +
  theme_bw()

#226072
ggplot() +
  geom_path(data = dat2[['226072']], aes(x, y), size = 0.5, alpha = 0.5) +
  geom_sf(data = contours2 %>%
            filter(id == 226072), aes(color = UD.level), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  labs(title = 'ID 226072') +
  theme_bw()

#41614
ggplot() +
  geom_path(data = dat2[['41614']], aes(x, y), size = 0.5, alpha = 0.5) +
  geom_sf(data = contours2 %>%
            filter(id == 41614), aes(color = UD.level), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  labs(title = 'ID 41614') +
  theme_bw()




# Plot all tracks w/ UD contours
ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat, aes(x, y, color = factor(id), group = id), size = 0.5, alpha = 0.5) +
  geom_sf(data = contours2, fill = "transparent", size = 0.75, color = 'black') +
  scale_color_viridis_d(option = 'turbo') +
  theme_bw() +
  facet_wrap(~ UD.level, nrow = 2) +
  coord_sf(xlim = c(min(dat$x) - 50, max(dat$x) + 50),
           ylim = c(min(dat$y) - 50, max(dat$y) + 50))






#### Export datasets for easy loading ####

save(contours2, file = "Processed_data/dBBMM_fits.RData")
