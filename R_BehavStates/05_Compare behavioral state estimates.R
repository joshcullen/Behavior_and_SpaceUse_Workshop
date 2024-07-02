
#######################################################
### Compare behavioral state estimates among models ###
#######################################################

library(tidyverse)
library(lubridate)
library(aniMotum)
library(momentuHMM)
library(bayesmove)
library(sf)
library(rnaturalearth)
library(plotly)



#### Load the model results from each method ####

load("Processed_data/SSM_model_fits.RData")
load("Processed_data/HMM_data_and_model_fits.RData")
load("Processed_data/bayesmove_model_fits.RData")



#### Wrangle model results to compile into data.frame ####

ssm_res<- join(ssm = fit_crw_8hr,
               mpm = fit_crw_mpm_8hr,
               what.ssm = "predicted")

hmm_res <- dat3 %>%
  mutate(state = viterbi(fit_hmm_3states_3vars))

bayes_res <- dat.out



#### Compare state-dependent distributions from HMM and M4 (bayesmove) ####

# HMM
plot(fit_hmm_3states_3vars, plotTracks = FALSE)


# M4
behav.res.seg2 <- behav.res.seg %>%
  mutate(behav1 = case_when(behav == 1 ~ 'Breeding',
                           behav == 2 ~ 'Foraging',
                           behav == 3 ~ 'Migratory',
                           TRUE ~ behav)) %>%
  filter(!behav %in% c(4:7)) %>%
  mutate(across(behav1, factor, levels = c('Breeding','Foraging','Migratory')))


ggplot(behav.res.seg2, aes(x = bin.vals, y = prop, fill = behav1)) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(size = 10)) +
  scale_fill_viridis_d(guide = 'none') +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  facet_grid(behav1 ~ var, scales = "free_x")




#### Viz time series of state estimates ####

# SSM
plot(fit_crw_mpm_8hr)

# HMM
plotStates(fit_hmm_3states_3vars)

# M4
ggplot(theta.estim.long) +
  geom_area(aes(x=date, y=prop, fill = behavior), color = "black", linewidth = 0.25,
            position = "fill") +
  labs(x = "\nTime", y = "Proportion of Behavior\n") +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id, scales = "free_x")



#### Viz map of behavioral state estimates across methods ####

brazil <- ne_countries(scale = 50, country = "brazil", returnclass = 'sf') %>%
  st_transform(crs = "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs")



## Compare w/ focal PTT 205537

# SSM
ssm_res_205537 <- filter(ssm_res, id == 205537)

ggplotly(
  ggplot() +
    # geom_sf(data = brazil) +
    geom_path(data = ssm_res_205537, aes(x=x, y=y), color="grey60", size=0.25) +
    geom_point(data = ssm_res_205537, aes(x, y, color=g), size=1.5, alpha=0.7) +
    geom_point(data = ssm_res_205537 %>%
                 slice(which(row_number() == 1)), aes(x, y), color = "green", pch = 21,
               size = 3, stroke = 1.25) +
    geom_point(data = ssm_res_205537 %>%
                 slice(which(row_number() == n())), aes(x, y), color = "red", pch = 24,
               size = 3, stroke = 1.25) +
    scale_color_viridis_c("SSM Behavior") +
    labs(x = "Easting", y = "Northing") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) #+
    # coord_sf(xlim = c(min(ssm_res_205537$x - 50), max(ssm_res_205537$x + 50)),
    #          ylim = c(min(ssm_res_205537$y - 50), max(ssm_res_205537$y + 50)))
)


# HMM
hmm_res_205537 <- hmm_res %>%
  filter(ID == 205537) %>%
  mutate(state1 = case_when(state == 1 ~ 'Breeding',
                            state == 2 ~ 'Foraging',
                            state == 3 ~ 'Migratory'))

ggplotly(
  ggplot() +
    # geom_sf(data = brazil) +
    geom_path(data = hmm_res_205537, aes(x=x, y=y), color="grey60", size=0.25) +
    geom_point(data = hmm_res_205537, aes(x, y, color=state1), size=1.5, alpha=0.7) +
    geom_point(data = hmm_res_205537 %>%
                 slice(which(row_number() == 1)), aes(x, y), color = "green", pch = 21,
               size = 3, stroke = 1.25) +
    geom_point(data = hmm_res_205537 %>%
                 slice(which(row_number() == n())), aes(x, y), color = "red", pch = 24,
               size = 3, stroke = 1.25) +
    scale_color_viridis_d("HMM Behavior") +
    labs(x = "Easting", y = "Northing") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) #+
    # coord_sf(xlim = c(min(hmm_res_205537$x - 50), max(hmm_res_205537$x + 50)),
    #          ylim = c(min(hmm_res_205537$y - 50), max(hmm_res_205537$y + 50)))
)


# M4
bayes_res_205537 <- bayes_res %>%
  filter(id == 205537)

ggplotly(
  ggplot() +
    # geom_sf(data = brazil) +
    geom_path(data = bayes_res_205537, aes(x=x, y=y), color="grey60", size=0.25) +
    geom_point(data = bayes_res_205537, aes(x, y, color=behav), size=1.5, alpha=0.7) +
    geom_point(data = bayes_res_205537 %>%
                 slice(which(row_number() == 1)), aes(x, y), color = "green", pch = 21,
               size = 3, stroke = 1.25) +
    geom_point(data = bayes_res_205537 %>%
                 slice(which(row_number() == n())), aes(x, y), color = "red", pch = 24,
               size = 3, stroke = 1.25) +
    scale_color_viridis_d("Bayesian M4 Behavior") +
    labs(x = "Easting", y = "Northing") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) #+
    # coord_sf(xlim = c(min(bayes_res_205537$x - 50), max(bayes_res_205537$x + 50)),
    #          ylim = c(min(bayes_res_205537$y - 50), max(bayes_res_205537$y + 50)))
)

ggplotly(
  ggplot() +
    # geom_sf(data = brazil) +
    geom_path(data = bayes_res_205537, aes(x=x, y=y), color="grey60", size=0.25) +
    geom_point(data = bayes_res_205537, aes(x, y, color=Foraging), size=1.5, alpha=0.7) +
    geom_point(data = bayes_res_205537 %>%
                 slice(which(row_number() == 1)), aes(x, y), color = "green", pch = 21,
               size = 3, stroke = 1.25) +
    geom_point(data = bayes_res_205537 %>%
                 slice(which(row_number() == n())), aes(x, y), color = "red", pch = 24,
               size = 3, stroke = 1.25) +
    scale_color_viridis_c("Bayesian M4 Behavior", option = 'inferno', end = 0.95) +
    labs(x = "Easting", y = "Northing") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) #+
    # coord_sf(xlim = c(min(bayes_res_205537$x - 50), max(bayes_res_205537$x + 50)),
    #          ylim = c(min(bayes_res_205537$y - 50), max(bayes_res_205537$y + 50)))
)
