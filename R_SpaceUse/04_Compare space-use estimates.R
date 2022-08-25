
#################################################
### Compare space-use estimates among methods ###
#################################################

library(tidyverse)
library(lubridate)
library(amt)  #v0.1.7
library(move)  #v4.0.6
library(sf)  #v1.0.7
library(rnaturalearth)
library(MetBrewer)
library(units)



#### Load the model results from each method ####

load("Processed_data/MCP_fits.RData")
load("Processed_data/KDE_fits.RData")
load("Processed_data/dBBMM_fits.RData")

dat <- read.csv('Processed_data/SSM_mp8hr_FDN Cmydas tracks.csv')


#### Wrangle model results to match up properly ####

# Transform project for MCP to match other output
mcp <- dat.id.mcp %>%
  st_transform(crs = st_crs(contours2)) %>%
  mutate(method = 'MCP') %>%
  dplyr::select(id, level, method, geometry)

# Rename other output
kde.href <- dat.id.kde.href %>%
  mutate(method = 'KDE_href') %>%
  dplyr::select(id, level, method, geometry)
kde.hpi <- dat.id.kde.hpi %>%
  mutate(method = 'KDE_hpi') %>%
  dplyr::select(id, level, method, geometry)
dbbmm <- contours2 %>%
  rename(level = UD.level) %>%
  mutate(method = 'dBBMM') %>%
  dplyr::select(id, level, method, geometry)




#### Visually compare UDs among methods ####

brazil<- ne_countries(scale = 50, country = "Brazil", returnclass = 'sf') %>%
  st_transform(crs = "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs") %>%
  dplyr::select(-level)

ud.all <- rbind(mcp, kde.href, kde.hpi, dbbmm)


# Show all IDs, methods, and levels
ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat, aes(x, y, group = id), size = 0.5, alpha = 0.5) +
  geom_sf(data = ud.all, aes(color = method), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  coord_sf(xlim = c(min(dat$x) - 50, max(dat$x) + 50),
           ylim = c(min(dat$y) - 50, max(dat$y) + 50)) +
  facet_grid(id ~ level)
#large disparity among methods for migrating IDs

# Zoom in on resident IDs
ggplot() +
  geom_path(data = dat %>%
              filter(id %in% c(205544,226072)), aes(x, y, group = id), size = 0.5, alpha = 0.15) +
  geom_sf(data = ud.all %>%
            filter(id %in% c(205544,226072)), aes(color = method), fill = "transparent", size = 0.5) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  facet_grid(id ~ level)
#all methods relatively comparable here






#### Compare estimated area of space-use ####

ud.all$area <- st_area(ud.all)
ud.all$strategy <- ifelse(ud.all$id %in% c(205544,226072), 'Resident', 'Migratory')


# Compare by method and UD level
ggplot(ud.all, aes(factor(level), area, color = method)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, position = position_dodge(0.55)) +
  geom_point(alpha = 0.7, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.55)) +
  scale_color_met_d("Hokusai3") +
  theme_bw()


# Compare by movement strategy
ggplot(ud.all %>%
         filter(level == 0.95), aes(strategy, area, color = method)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, position = position_dodge(0.55)) +
  geom_point(alpha = 0.7, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.55)) +
  scale_color_met_d("Hokusai3") +
  theme_bw()

# Focus on Residents
ggplot(ud.all %>%
         filter(strategy == 'Resident'), aes(factor(level), area, color = method)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, position = position_dodge(0.55)) +
  geom_point(alpha = 0.7, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.55)) +
  scale_color_met_d("Hokusai3") +
  theme_bw()
