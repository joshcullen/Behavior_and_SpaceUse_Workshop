
################################################
### Calculate minimum convex polygons (MCPs) ###
################################################

library(tidyverse)
library(lubridate)
library(amt)
library(sf)
library(rnaturalearth)
library(plotly)


#### Load data ####

dat <- read.csv('Processed_data/SSM_mp8hr_FDN Cmydas tracks.csv')

glimpse(dat)
summary(dat)


#### Wrangle and prep data for estimation of MCPs ####

dat <- dat %>%
  mutate(date = as_datetime(date))

dat.track <- make_track(dat, lon, lat, date, crs = 4326, all_cols = TRUE)



#### Calculate MCPs for all IDs ####

dat.mcp <- hr_mcp(dat.track, levels = c(0.5, 0.95, 1))
dat.mcp
plot(dat.mcp, col = c('red','green','blue'))


brazil<- ne_countries(scale = 50, country = "Brazil", returnclass = 'sf')

ggplot() +
  # geom_sf(data = brazil) +
  geom_point(data = dat, aes(lon, lat), alpha = 0.25, size = 1) +
  geom_sf(data = dat.mcp$mcp, aes(color = factor(level)), fill = 'transparent', linewidth = 0.75) +
  scale_color_viridis_d(direction = -1) +
  theme_bw() #+
  # coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))






#### Calculate MCPs per ID (at 95% level) ####

# Need to turn data.frame into list to map hr_mcp() function and then recombine
dat.id.mcp <- dat.track %>%
  split(.$id) %>%
  map(hr_mcp, levels = 0.95) %>%
  map(pluck, 1) %>%
  do.call(rbind, .)

dat.id.mcp <- dat.id.mcp %>%
  mutate(id = rownames(.), .before = level)


ggplotly(
  ggplot() +
    # geom_sf(data = brazil) +
    geom_point(data = dat, aes(lon, lat, color = factor(id)), alpha = 0.25, size = 1) +
    geom_sf(data = dat.id.mcp, aes(color = id), fill = 'transparent', linewidth = 0.75) +
    scale_color_viridis_d() +
    theme_bw() #+
    # coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)




#### Export datasets for easy loading ####

save(dat.id.mcp, file = "Processed_data/MCP_fits.RData")
