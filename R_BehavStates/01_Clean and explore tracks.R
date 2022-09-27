
##############################
### Clean and explore data ###
##############################

library(tidyverse)
library(lubridate)
library(plotly)
library(bayesmove)
library(rnaturalearth)
library(sf)


#### Load data ####

dat <- read.csv("Raw_data/FDN Cmydas tracks_subset.csv")

# Explore data summaries
glimpse(dat)
summary(dat)
n_distinct(dat$Ptt)

## Structure of dataset ##
#PTT: animal ID (Platform Transmitter Terminal; PTT)
#Date: datetime of each relocation
#Longitude: longitude of relocation
#Latitude: latitude of relocation
#Type: label defining whether observation was transmitted via Fastloc GPS or Argos doppler-shift
#Quality: the number of satellites for FastGPS observations and the location class (LC) for Argos observations
#Error.radius: the radius of the estimated error for each observation (in meters)
#Error.Semi.major.axis: the distance (in meters) along the semi-major axis of the estimated error ellipse
#Error.Semi.minor.axis: the distance (in meters) along the semi-minor axis of the estimated error ellipse
#Error.Ellipse.orientation: the angle (in degrees) defining the orientation of the estimated error ellipse


# Inspect columns of pertinent variables
table(dat$Ptt, useNA = "ifany")
table(dat$Type, useNA = "ifany")
table(dat$Quality, useNA = "ifany")

#After inspection of the data, we'll need to convert the 'Date' column into a datetime format and remove any observations where Quality is either 'Z' or NA



#### Filter data based on 'Quality' ####

dat2 <- dat %>%
  filter(Quality != "Z" & !is.na(Quality))

table(dat2$Quality, useNA = "ifany")  #Quality looks good now



#### Convert 'Date' to datetime format ####

dat2$Date[1:10]  #see what current format looks like

dat2 <- dat2 %>%
  mutate(Date = as_datetime(Date, format = '%Y-%m-%d %H:%M:%S'))
dat2$Date[1:10]  #inspect modified format




#### Visualize tracks to determine if any anomalies ####

# Change Ptt to character so treated as discrete quantities
dat2$Ptt <- as.character(dat2$Ptt)

ggplot(dat2, aes(Longitude, Latitude, color = Ptt)) +
  geom_path(aes(group = Ptt), size = 0.25) +
  scale_color_viridis_d() +
  theme_bw()
#there seems to be some testing points remaining from Wildlife Computers HQ in Seattle


# Let's add some land layers to give better context
world <- rnaturalearth::ne_countries(scale = 50, continent = c("North America","South America"),
                                     returnclass = "sf")

ggplot() +
  geom_sf(data = world) +
  geom_path(data = dat2, aes(Longitude, Latitude, group = Ptt, color = Ptt), size = 0.25) +
  scale_color_viridis_d() +
  theme_bw() +
  coord_sf(xlim = c(min(dat2$Longitude) - 1, max(dat2$Longitude) + 1),
           ylim = c(min(dat2$Latitude) - 1, max(dat2$Latitude) + 1))


# Let's make this map interactive so we can see a little better
plotly::ggplotly(
  ggplot() +
    geom_sf(data = world) +
    geom_path(data = dat2, aes(Longitude, Latitude, group = Ptt, color = Ptt), size = 0.25) +
    scale_color_viridis_d() +
    theme_bw() +
    coord_sf(xlim = c(min(dat2$Longitude) - 1, max(dat2$Longitude) + 1),
             ylim = c(min(dat2$Latitude) - 1, max(dat2$Latitude) + 1))
)
#there appear to be a handful of IDs that have spiked trajectories or unlikely movements that we'll need to filter


# Let's use a Shiny app from {bayesmove} to explore a little further

dat2 %>%
  rename(id = Ptt, date = Date, x = Longitude, y = Latitude) %>%
  bayesmove::shiny_tracks(epsg = 4326)



#### Remove locations likely before tag deployment ####

# Filter by study extent
dat3<- dat2 %>%
  filter(Latitude < -3 & Longitude > -40)

ggplot() +
  geom_path(data = dat3, aes(Longitude, Latitude, group = Ptt, color = Ptt), size = 0.25) +
  scale_color_viridis_d() +
  theme_bw()
# looks better; unfiltered spikes will be removed later

# Inspect time series plots of lat and long
ggplot() +
  geom_line(data = dat3, aes(Date, Longitude, color = Ptt)) +
  scale_color_viridis_d() +
  theme_bw() +
  facet_wrap(~Ptt, scales = "free")

ggplot() +
  geom_line(data = dat3, aes(Date, Latitude, color = Ptt)) +
  scale_color_viridis_d() +
  theme_bw() +
  facet_wrap(~Ptt, scales = "free")
# nothing else needs to be dealt with at this point




#### Export cleaned data ####

write.csv(dat3, "Processed_data/Cleaned_FDN Cmydas tracks.csv", row.names = FALSE)
