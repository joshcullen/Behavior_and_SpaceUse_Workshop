
####################################
### Calculate autocorrelated KDE ###
####################################

library(tidyverse)
library(lubridate)
library(ctmm)  #v1.0.0
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

