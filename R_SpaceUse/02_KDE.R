
################################################
### Calculate kernel density estimates (KDE) ###
################################################

library(tidyverse)
library(lubridate)
library(amt)  #v0.1.7
library(sf)  #v1.0.7
library(rnaturalearth)
library(plotly)
library(tictoc)
library(MetBrewer)
library(units)



### Load data ###

dat <- read.csv('Processed_data/SSM_mp_FDN_irreg Cmydas tracks.csv')

glimpse(dat)
summary(dat)


### Wrangle and prep data for KDE ###

dat <- dat %>%
  mutate(date = as_datetime(date))

dat.track <- make_track(dat, x, y, date, crs = "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs",
                        all_cols = TRUE)

# Create template raster for KDE
trast <- make_trast(dat.track, res = 4)





### Calculate KDE for all IDs ###

## Href (reference bandwidth method)
dat.kde.ref <- hr_kde(dat.track, trast = trast, h = hr_kde_ref(dat.track), levels = c(0.5, 0.95))
dat.kde.ref
plot(dat.kde.ref, col = c("red", "blue"))

kde.href.contours <- hr_isopleths(dat.kde.ref)


brazil<- ne_countries(scale = 50, country = "Brazil", returnclass = 'sf')

# Plot returned raster layer (may not be ideal in some instances)
ggplot() +
  geom_tile(data = raster::as.data.frame(dat.kde.ref$ud, xy = TRUE), aes(x, y, fill = layer)) +
  geom_sf(data = st_transform(brazil, crs = "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs")) +
  geom_point(data = dat, aes(x, y), alpha = 0.25, size = 1, color = "chartreuse") +
  scale_fill_viridis_c(option = "inferno") +
  theme_bw() +
  coord_sf(xlim = c(min(dat$x) - 50, max(dat$x) + 50),
           ylim = c(min(dat$y) - 50, max(dat$y) + 50))


# Plot isopleth contours as requested levels
ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat, aes(lon, lat, group = id), alpha = 0.25, size = 0.3) +
  geom_sf(data = kde.href.contours, aes(color = factor(level)), fill = 'transparent', size = 0.75) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))





## Hpi (plug-in method)
h_pi <- hr_kde_pi(dat.track, rescale = "xvar")
dat.kde.pi <- hr_kde(dat.track, trast = trast, h = h_pi, levels = c(0.5, 0.95))
dat.kde.pi
plot(dat.kde.pi, col = c("red", "blue"))

kde.hpi.contours <- hr_isopleths(dat.kde.pi)


# Plot isopleth contours as requested levels
ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat, aes(lon, lat, group = id), alpha = 0.25, size = 0.3) +
  geom_sf(data = kde.hpi.contours, aes(color = factor(level)), fill = 'transparent', size = 0.75) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))







### Calculate KDE per ID ###


## Href

# Need to turn data.frame into list to map hr_kde() function and then recombine
tic()
dat.id.kde.href <- dat.track %>%
  split(.$id) %>%
  map(~hr_kde(.x,
              trast = make_trast(.x, res = 0.5),
              h = hr_kde_ref(.x),
              levels = c(0.5, 0.95))
      ) %>%
  map(hr_isopleths) %>%
  do.call(rbind, .)
toc()  #takes 2 min to run
BRRR::skrrrahh('ross1')

dat.id.kde.href <- dat.id.kde.href %>%
  mutate(id = rownames(.), .before = level) %>%
  mutate(id = str_replace(id, "\\..$", ""))  #remove decimal and extra number


ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat, aes(lon, lat, group = id), alpha = 0.25, size = 0.3) +
  geom_sf(data = dat.id.kde.href, aes(color = factor(level)), fill = 'transparent', size = 0.75) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  coord_sf(xlim = c(-44, -30), ylim = c(-9, 0)) +
  facet_wrap(~ id)





## Hpi

# Need to turn data.frame into list to map hr_kde() function and then recombine
dat.id.list <- dat.track %>%
  split(.$id) %>%
  map(~{.x %>%
      mutate(pattern = ifelse(min(lon) > -34, 'Resident', 'Migratory'))
    })

id.pattern <- dat.id.list %>%
  map(pluck, 'pattern') %>%
  map(~{.x[1]})

res.ind <- which(id.pattern == 'Resident')
mig.ind <- which(id.pattern == 'Migratory')


# Run for Migratory IDs
tic()
dat.mig.kde.hpi <- dat.id.list[mig.ind] %>%
  map(~hr_kde(.x,
              trast = make_trast(.x, res = 1),
              h = hr_kde_pi(.x, rescale = 'xvar'),
              levels = c(0.5, 0.95))
  ) %>%
  map(hr_isopleths) %>%
  do.call(rbind, .)
toc()  #takes 22 sec to run
BRRR::skrrrahh('khaled3')


# Run for Resident IDs
dat.res.kde.hpi <- dat.id.list[res.ind] %>%
  map(~hr_kde(.x,
              trast = make_trast(.x, res = 0.1),
              h = hr_kde_pi(.x, rescale = 'xvar'),
              levels = c(0.5, 0.95))
  ) %>%
  map(hr_isopleths) %>%
  do.call(rbind, .)
BRRR::skrrrahh('liljon')

dat.id.kde.hpi <- rbind(dat.mig.kde.hpi, dat.res.kde.hpi) %>%
  mutate(id = rownames(.), .before = level) %>%
  mutate(id = str_replace(id, "\\..$", ""))  #remove decimal and extra number


ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat, aes(lon, lat, group = id), alpha = 0.25, size = 0.3) +
  geom_sf(data = dat.id.kde.hpi, aes(color = factor(level)), fill = 'transparent', size = 0.75) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  coord_sf(xlim = c(-42, -32), ylim = c(-8, -2)) +
  facet_wrap(~ id)







### Compare KDE estimates between methods for bandwidth estimation ###

# Collate both ID-level data.frames together
dat.id.kde.href$method <- 'href'
dat.id.kde.hpi$method <- 'hpi'

dat.id.kde <- rbind(dat.id.kde.href, dat.id.kde.hpi)

ggplot(dat.id.kde, aes(factor(level), area, color = method)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, position = position_dodge(0.55)) +
  geom_point(alpha = 0.7, position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.55)) +
  scale_color_met_d("Hokusai3") +
  theme_bw()

