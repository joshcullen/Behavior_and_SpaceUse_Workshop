
######################################################################################
### Fit non-parametric Bayesian movement models to track segments and observations ###
######################################################################################

### Additional notes to discuss:
# * {bayesmove} can also estimate multiple discrete behavioral states using any number of covariates; this includes ancillary biologging data and evironmental variables
# * {bayesmove} does not require the user to specify the number of states they'd like to fit; just enter the maximum number they think would be possible to detect
# * Option to pre-specify breakpoints for the RJMCMC to consider for the segmentation model
# * All data streams must be discretized into bins before analysis
# * The mixture model resembles an HMM, but does not have a Markov property (i.e., state estimates are independent of each other)
# * Can more easily accommodate data streams that are difficult to select a parametric PDF for (and set intial params for) compared to HMM
# * Built-in parallelization of models, so runtime is much quicker "out-of-the-box" than some other models

library(tidyverse)
library(lubridate)
library(bayesmove)  #v0.2.1
library(sf)  #v1.0.7
library(rnaturalearth)
library(plotly)
library(furrr)
library(future)


### Load data ###

dat <- read.csv('Processed_data/SSM_mp6hr_FDN Cmydas tracks.csv')

glimpse(dat)
summary(dat)




### Wrangle data for analysis using {bayesmove} ###

# Convert all 'date' to datetime format
dat <- dat %>%
  mutate(date = as_datetime(date))

dat2<- prep_data(dat = dat, coord.names = c("x","y"), id = "id")
head(dat2)
#since x and y are in km, steps and NSD are in km
# calculates step lengths, turning angles, net-squared displacement (NSD), and time step (dt)


# Let's double-check that all time-steps are at 6 hrs (21600 s)
table(dat2$dt)  #yes


# Since we don't need to filter out obs at other time intervals, we still need to add required variables to data.frame
dat2 <- dat2 %>%
  group_by(id) %>%  #need to number rows separately for each ID
  mutate(time1 = 1:n(),
         obs = 1:n()) %>%
  ungroup()

#verify that it worked properly
dat2 %>%
  dplyr::select(id, date, time1, obs) %>%   #select only a few cols since tibble hides time1 and obs
  split(.$id) %>%
  head()


# For direct comparison w/ HMM results, create displacement variable
dat2$disp <- sqrt(dat2$NSD)






### Discretize data streams for models ###

# Viz density plots of each data stream
ggplot(dat2) +
  geom_density(aes(step), fill = "cadetblue") +
  theme_bw()

ggplot(dat2) +
  geom_density(aes(angle), fill = "firebrick") +
  theme_bw()

ggplot(dat2) +
  geom_density(aes(disp), fill = "goldenrod") +
  # scale_x_continuous(breaks = seq(0, 1000, by = 200)) +
  theme_bw()



# Define bin limits (and number of bins)

# turning angle (naturally constrained in [0,2*pi] or [-pi,+pi])
angle.bin.lims <- seq(from = -pi, to = pi, by = pi/4)  #8 bins

# step length (must be positive, but no upper bound)
step.bin.lims <- c(seq(from = 0, to = 5, length = 6), max(dat2$step, na.rm = TRUE))  #6 bins

# displacement (must be positive, but no upper bound)
disp.bin.lims <- c(seq(from = 0, to = 750, by = 250), max(dat2$disp, na.rm = TRUE))  #4 bins


angle.bin.lims
step.bin.lims
disp.bin.lims


# Discretize data streams
dat.disc <- discrete_move_var(dat2,
                              lims = list(step.bin.lims, angle.bin.lims, disp.bin.lims),
                              varIn = c("step","angle","disp"),
                              varOut = c("SL","TA","Disp"))



# Viz histograms of discretized data streams
ggplot(dat.disc) +
  geom_bar(aes(SL), fill = "cadetblue") +
  theme_bw()

ggplot(dat.disc) +
  geom_bar(aes(TA), fill = "firebrick") +
  theme_bw()

ggplot(dat.disc) +
  geom_bar(aes(Disp), fill = "goldenrod") +
  theme_bw()








### Fit observation-level mixture model to estimate states ###

# Only retain columns of discretized data streams
dat.disc.sub<- dat.disc %>%
  dplyr::select(SL, TA, Disp) %>%
  data.frame()   #cluster_obs() function crashes if trying to use 'tibble'


set.seed(123)

# Define model params
alpha = 0.1  #prior on Dirichlet distribution
ngibbs = 10000  #number of Gibbs sampler iterations
nburn = ngibbs/2  #number of burn-in iterations
nmaxclust = 7  #number of maximum possible states (clusters) present

# Run model
dat.res.obs<- cluster_obs(dat = dat.disc.sub, alpha = alpha, ngibbs = ngibbs, nmaxclust = nmaxclust,
                      nburn = nburn)
# took 9 min to run


# Inspect traceplot of log-likelihood
plot(dat.res.obs$loglikel, type = "l")
abline(v = nburn, col = "red", lwd = 2)


## Inspect and plot results
post.seq<- (nburn + 1):ngibbs  #posterior samples

theta<- dat.res.obs$theta[post.seq,]
# colnames(theta)<- 1:ncol(theta)
theta1<- colMeans(theta)
# theta1<- sort(theta1, decreasing = TRUE)
cumsum(theta1)  #possibly 4 states present; represents > 90% of all obs



# Extract bin estimates for each possible state from the `phi` matrix of the model results
behav.res<- get_behav_hist(dat = dat.res.obs, nburn = nburn, ngibbs = ngibbs, nmaxclust = nmaxclust,
                           var.names = c("Step Length","Turning Angle","Displacement"))
behav.res$behav<- factor(behav.res$behav, levels = 1:nmaxclust)

# Plot state-dependent distributions
ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c(viridis::viridis(4), rep("grey35", 3)), guide = 'none') +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  scale_x_continuous(breaks = 1:8) +
  facet_grid(behav ~ var, scales = "free_x")
##actually looks like there's 4 states, but that 3 & 4 should be merged and state 5 is a transiting state



## Assign behavioral states to observations

# Using MAP estimate, threshold of 75% assignments from posterior, and most common state
z.post<- as.data.frame(dat.res.obs$z.posterior)
z.post$`3` <- z.post[,3] + z.post[,4]  #combine states 3 and 4
z.post <- z.post[,-c(3:4)]  #remove original states 3 and  4
z.post <- relocate(z.post, `3`, .after = V2)  #reorder columns for new state 3
names(z.post) <- 1:ncol(z.post)  #rename columns/states

z.post2<- t(apply(z.post, 1, function(x) x/sum(x)))  #calculate proportions of samples from posterior distribution assigned to each state
thresh<- 0.75  #user-defined threshold percentage for classifying a state
z.post.thresh<- apply(z.post2, 1, function(x) ifelse(max(x) > thresh, which(x > thresh), NA))
z.post.max<- apply(z.post2, 1, function(x) which.max(x))
z.map <- ifelse(dat.res.obs$z.MAP == 4, 3,
                ifelse(dat.res.obs$z.MAP == 5, 4, dat.res.obs$z.MA))

## Add states to data frame
dat.states<- dat.disc %>%
  mutate(z.map = z.map,
         z.post.thresh = z.post.thresh,
         z.post.max = z.post.max)

n.states<- 4
dat.states$z.map<- ifelse(dat.states$z.map > n.states, NA, dat.states$z.map)
dat.states$z.post.max<- ifelse(dat.states$z.post.max > n.states, NA, dat.states$z.post.max)





# Assign names to states
dat.states2<- dat.states %>%
  mutate(across(c('z.map','z.post.thresh','z.post.max'),
                ~case_when(. == 1 ~ "Breeding_ARS",
                           . == 2 ~ "Foraging",
                           . == 3 ~ "Breeding_Encamped",
                           . == 4 ~ "Migratory",
                           is.na(.) ~ "Unclassified")
  )) %>%
  mutate(across(c('z.map','z.post.thresh','z.post.max'),
                factor, levels = c('Breeding_Encamped','Breeding_ARS','Foraging',
                                   'Migratory','Unclassified')
  ))



# Inspect number of obs assigned to each state

dat.states2 %>%   # for estimates based on MAP estimate
  group_by(z.map) %>%
  tally() %>%
  mutate(prop = n/sum(n))

dat.states2 %>%   # for estimates based on threshold on posterior
  group_by(z.post.thresh) %>%
  tally() %>%
  mutate(prop = n/sum(n))

dat.states2 %>%   # for estimates based on mode of posterior
  group_by(z.post.max) %>%
  tally() %>%
  mutate(prop = n/sum(n))



# Map results
brazil <- ne_countries(scale = 50, country = "brazil", returnclass = 'sf') %>%
  st_transform(crs = "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs")

# Using MAP estimates
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil) +
    geom_path(data = dat.states2, aes(x, y, group = id), color="grey60", size=0.25) +
    geom_point(data = dat.states2, aes(x, y, fill=z.map), size=1.5, pch=21, alpha=0.7) +
    geom_point(data = dat.states2 %>%
                 group_by(id) %>%
                 slice(which(row_number() == 1)) %>%
                 ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
    geom_point(data = dat.states2 %>%
                 group_by(id) %>%
                 slice(which(row_number() == n())) %>%
                 ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
    scale_fill_manual("Behavior",
                      values = c(viridis::viridis(4), "grey50")) +
    labs(x = "Easting", y = "Northing", title = "MAP estimate") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) +
    guides(fill = guide_legend(label.theme = element_text(size = 12),
                               title.theme = element_text(size = 14))) +
    coord_sf(xlim = c(min(dat.states2$x - 50), max(dat.states2$x + 50)),
             ylim = c(min(dat.states2$y - 50), max(dat.states2$y + 50)))
)

# Using estimates w/ threshold on posterior
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil) +
    geom_path(data = dat.states2, aes(x, y, group = id), color="grey60", size=0.25) +
    geom_point(data = dat.states2, aes(x, y, fill=z.post.thresh), size=1.5, pch=21, alpha=0.7) +
    geom_point(data = dat.states2 %>%
                 group_by(id) %>%
                 slice(which(row_number() == 1)) %>%
                 ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
    geom_point(data = dat.states2 %>%
                 group_by(id) %>%
                 slice(which(row_number() == n())) %>%
                 ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
    scale_fill_manual("Behavior",
                      values = c(viridis::viridis(4), "grey50")) +
    labs(x = "Easting", y = "Northing", title = "Threshold on posterior") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) +
    guides(fill = guide_legend(label.theme = element_text(size = 12),
                               title.theme = element_text(size = 14))) +
    coord_sf(xlim = c(min(dat.states2$x - 50), max(dat.states2$x + 50)),
             ylim = c(min(dat.states2$y - 50), max(dat.states2$y + 50)))
)


# Using estimates w/ most common state from posterior
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil) +
    geom_path(data = dat.states2, aes(x, y, group = id), color="grey60", size=0.25) +
    geom_point(data = dat.states2, aes(x, y, fill=z.post.max), size=1.5, pch=21, alpha=0.7) +
    geom_point(data = dat.states2 %>%
                 group_by(id) %>%
                 slice(which(row_number() == 1)) %>%
                 ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
    geom_point(data = dat.states2 %>%
                 group_by(id) %>%
                 slice(which(row_number() == n())) %>%
                 ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
    scale_fill_manual("Behavior",
                      values = c(viridis::viridis(4), "grey50")) +
    labs(x = "Easting", y = "Northing", title = "Most common state from posterior") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) +
    guides(fill = guide_legend(label.theme = element_text(size = 12),
                               title.theme = element_text(size = 14))) +
    coord_sf(xlim = c(min(dat.states2$x - 50), max(dat.states2$x + 50)),
             ylim = c(min(dat.states2$y - 50), max(dat.states2$y + 50)))
)


dat.states2 %>%
  mutate(across(z.map:z.post.max, as.numeric)) %>%  #can only viz time series of numeric vars
  shiny_tracks(., epsg = "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs")








### Fit segment-level mixed-membership model to estimate states ###

# Convert data to list by ID
dat.list <- dat.disc %>%
  df_to_list(., "id")

# Only retain id and discretized data streams
dat.list.sub<- map(dat.list,
                   subset,
                   select = c(id, SL, TA, Disp))



# Run the segmentation model (unsupervised)
set.seed(123)

alpha<- 1  # hyperparameter for prior (Dirichlet) distribution
ngibbs<- 50000  # number of iterations for Gibbs sampler
nbins<- c(6,8,4)  # define number of bins per data stream (in order from dat.list.sub)

progressr::handlers(progressr::handler_progress(clear = FALSE))  #to initialize progress bar
future::plan(multisession, workers = availableCores() - 2)  #run MCMC chains in parallel

dat.res.seg1<- segment_behavior(data = dat.list.sub, ngibbs = ngibbs, nbins = nbins,
                           alpha = alpha)
future::plan(future::sequential)  #return to single core
# takes 1.5 min to run



# Trace-plots for the number of breakpoints per ID
traceplot(data = dat.res.seg1, type = "nbrks")

# Trace-plots for the log marginal likelihood (LML) per ID
traceplot(data = dat.res.seg1, type = "LML")  #appears to have converged for each track






# Run the segmentation model (semi-supervised via pre-specification of breakpoints)

#Calculate the difference in displacement between subsequent steps
dat.list <- dat.list %>%
  purrr::map(., ~{.x %>%
      mutate(disp_diff = c(diff(disp), NA))
  })

#Pre-define these migratory phases
dat.list <- dat.list %>%
  map(., ~{.x %>%
      mutate(phase = case_when(disp < 6 ~ 1,
                               disp_diff > 2 ~ 2,
                               disp > 6 & disp_diff < 2 ~ 3)
      )
  })


#Find breakpoints based on 'phase'
breaks<- map(dat.list, ~find_breaks(dat = ., ind = "phase"))
breaks  #since some IDs have 0 estimated breaks and model needs at least 1 for all IDs, provide 1 fake brkpt

ind <- which(lengths(breaks) == 0)
breaks[ind] <- 1
breaks


set.seed(123)

alpha<- 1  # hyperparameter for prior (Dirichlet) distribution
ngibbs<- 50000  # number of iterations for Gibbs sampler
nbins<- c(6,8,4)

progressr::handlers(progressr::handler_progress(clear = FALSE))  #to initialize progress bar
future::plan(multisession, workers = availableCores() - 2)  #run MCMC chains in parallel

dat.res.seg2<- segment_behavior(data = dat.list.sub, ngibbs = ngibbs, nbins = nbins,
                                alpha = alpha, breakpt = breaks)
future::plan(future::sequential)  #return to single core
# takes 3 min to run



# Trace-plots for the number of breakpoints per ID
traceplot(data = dat.res.seg2, type = "nbrks")

# Trace-plots for the log marginal likelihood (LML) per ID
traceplot(data = dat.res.seg2, type = "LML")  #appears to have converged for each track




# Determine MAP for selecting breakpoints
MAP.est1<- get_MAP(dat = dat.res.seg1$LML, nburn = ngibbs/2)
brkpts1<- get_breakpts(dat = dat.res.seg1$brkpts, MAP.est = MAP.est1)

MAP.est2<- get_MAP(dat = dat.res.seg2$LML, nburn = ngibbs/2)
brkpts2<- get_breakpts(dat = dat.res.seg2$brkpts, MAP.est = MAP.est2)

# How many breakpoints estimated per ID?
apply(brkpts1[,-1], 1, function(x) length(purrr::discard(x, is.na)))
apply(brkpts2[,-1], 1, function(x) length(purrr::discard(x, is.na)))

brkpts1
brkpts2
#looks like segmentation results were entirely unchanged by pre-specified brkpts; good sign that we reached best fit model


# Viz breakpoints w/ respect to data streams
plot_breakpoints(data = dat.list, as_date = TRUE, var_names = c("step","angle","disp"),
                 var_labels = c("Step Length (km)", "Turning Angle (rad)", "Displacement (km)"),
                 brkpts = brkpts1)



# Assign track segments to each ID
dat.seg<- assign_tseg(dat = dat.list, brkpts = brkpts1)

head(dat.seg)








# Cluster segments into behavioral states

#Select only id, tseg, and discretized data streams
dat.seg2<- dat.seg[,c("id","tseg","SL","TA","Disp")]

#Summarize observations by track segment
nbins<- c(6,8,4)
obs<- summarize_tsegs(dat = dat.seg2, nbins = nbins)
obs


set.seed(123)

# Prepare for Gibbs sampler
ngibbs<- 5000  #number of MCMC iterations for Gibbs sampler
nburn<- ngibbs/2  #number of iterations for burn-in
nmaxclust<- 7  #same as used for mixture model on observations
ndata.types<- length(nbins)  #number of data types

# Set priors for LDA clustering model
gamma1<- 0.1
alpha<- 0.1

# Run LDA model
dat.res.segclust<- cluster_segments(dat = obs, gamma1 = gamma1, alpha = alpha,
                                    ngibbs = ngibbs, nmaxclust = nmaxclust,
                                    nburn = nburn, ndata.types = ndata.types)
# takes ~1 min to run


# Check traceplot of log likelihood
plot(dat.res.segclust$loglikel, type='l', xlab = "Iteration", ylab = "Log Likelihood")
abline(v = nburn, col = "red", lwd = 2)


#Determine likely number of states from proportion assigned to each segment
theta.estim<- extract_prop(res = dat.res.segclust, ngibbs = ngibbs, nburn = nburn,
                           nmaxclust = nmaxclust)

theta.estim_df<- theta.estim %>%
  as.data.frame() %>%
  pivot_longer(., cols = 1:nmaxclust, names_to = "behavior", values_to = "prop") %>%
  modify_at("behavior", factor)
levels(theta.estim_df$behavior)<- 1:nmaxclust

ggplot(theta.estim_df, aes(behavior, prop)) +
  geom_boxplot(fill = "grey35", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "grey35", position = position_jitter(0.1),
              alpha = 0.3) +
  labs(x="\nBehavior", y="Proportion of Total Behavior\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
#likely at least 4, but possibly 5 states; need to viz state-dependent distributions


#Calculate mean proportions per behavior
(theta.means<- round(colMeans(theta.estim), digits = 3))

#Calculate cumulative sum
cumsum(theta.means)  #probably 4 or 5 states



# Extract bin estimates from phi matrix (for behavior distribs)
behav.res<- get_behav_hist(dat = dat.res.segclust, nburn = nburn, ngibbs = ngibbs,
                           nmaxclust = nmaxclust,
                           var.names = c("Step Length","Turning Angle","Displacement"))

# Plot histograms of proportion data
ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c(viridis::viridis(4), rep("grey35", 3)), guide = 'none') +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  scale_x_continuous(breaks = 1:8) +
  facet_grid(behav ~ var, scales = "free_x")




#Reformat proportion estimates for all track segments
theta.estim.long<- expand_behavior(dat = dat.seg, theta.estim = theta.estim, obs = obs, nbehav = 6,
                                   behav.names = c("Migratory", "Foraging1", "Breeding_ARS",
                                                   "Breeding_Encamped", "Foraging2", "Foraging3"),
                                   behav.order = c(4,3,1,2,5:6))

#Plot results
ggplot(theta.estim.long) +
  geom_area(aes(x=date, y=prop, fill = behavior), color = "black", size = 0.25,
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




## Assign states to segments and map

# Convert segmented dataset into list
dat.seg.list<- df_to_list(dat = dat.seg, ind = "id")

# Merge results with original data
dat.out<- assign_behavior(dat.orig = dat.seg,
                          dat.seg.list = dat.seg.list,
                          theta.estim.long = theta.estim.long,
                          behav.names = levels(theta.estim.long$behavior))

# Map dominant behavior for all IDs
ggplot() +
  geom_path(data = dat.out, aes(x=x, y=y), color="grey60", size=0.25) +
  geom_point(data = dat.out, aes(x, y, fill=behav), size=1.5, pch=21, alpha=0.7) +
  geom_point(data = dat.out %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat.out %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_fill_viridis_d("Behavior", na.value = "grey50") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(~id, scales = "free")



# Proportion of a given state (e.g., migratory and foraging)
ggplot() +
  geom_path(data = dat.out, aes(x, y, color = Migratory, group = id), size=0.5, alpha=0.7) +
  geom_point(data = dat.out %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat.out %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_distiller("Proportion\nMigratory", palette = "Spectral", na.value = "grey50") +
  labs(x = "Easting", y = "Northing", title = "Migratory") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

dat.out2 <- dat.out %>%
  mutate(Foraging = Foraging1 + Foraging2 + Foraging3)


ggplot() +
  geom_path(data = dat.out2, aes(x, y, color = Foraging, group = id), size=0.5, alpha=0.7) +
  geom_point(data = dat.out2 %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat.out2 %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_distiller("Proportion\nForaging", palette = "Spectral", na.value = "grey50") +
  labs(x = "Easting", y = "Northing", title = "Foraging") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
