# data accessing
library('here')      # for easier directory referencing, conflicts with lubridate::here
library('readr')     # to read in files as tibbles

# data editing
library('dplyr')     # makes data editing easier
library('tibble')    # a tibble is a fancy data.frame
library('lubridate') # makes working with dates smoother 

# model fitting
library('pammtools') # tools for Piecewise-exponential Additive Mixed Models
library('mgcv')      # to fit GAMs
library('brms')      # to fit censored Gamma HGAM
library('survival')  # survival analysis functions

# graphics
library('ggplot2')   # fancy plots
library('cowplot')   # ggplot in grids
library('gratia')    # pretty GAM plots
library('ggridges')  # for ridgeline plots
library('ggExtra')   # for marginal plots
library('ggthemes')  # for colorblind palette

source(here::here('functions/save.plt.R')) # to save plots easily

theme_set(theme_bw())

start.year <- 1800 # start year of plots

# plot labels
june.lab <- expression(Days~after~June~30^{th})
sept.lab <- expression(Days~after~September~30^{th})

# palette for plots
pal <- c('#4477AA', '#ff8c00', '#66CCEE', '#009900',
         '#CCBB44', '#EE6677', '#AA3377', '#BBBBBB')
plot(1:length(pal), col = pal, cex = 5, pch = 15)

# ice phenology data
ice <- readRDS(here::here('data', 'lake-ice-data.rds'))

# intermittent/continuous ice cover
NCC.data <- read.csv(here::here('data', 'LakeIceIncidencewithCharacteristics_Final.csv'))
set.seed(1)
NCC.data$y <- runif(1:nrow(NCC.data)) # y-axis values for plot

#### exploratory plots ####
# which lake has intermittent ice cover at lowest temperature
continuous <- subset(NCC.data, NCC.data$IntermittentIceCover == 'N')
mohansic <- continuous[which.max(continuous$MeanAnnualAirTemp_c), ]

# which lake has continuous ice cover at lowest temperature
intermittent <- subset(NCC.data, NCC.data$IntermittentIceCover == 'Y')
clark <- intermittent[which.min(intermittent$MeanAnnualAirTemp_c), ]

### intermittent ice
# log-lat
WorldData <-
  map_data('world') %>%
  filter(lat > 0) %>%
  fortify() %>%
  as_tibble()

plt.long.lag <- ggplot(NCC.data, aes(Longitude_dd, Latitude_dd, col = IntermittentIceCover)) +
  geom_map(map = WorldData, aes(group = group, map_id = region), data = WorldData,
           inherit.aes = FALSE, fill = 'grey75', col = 'black', size = 0.5) +
  coord_map('azequidistant', xlim = c(-180, 180), ylim = c(30, 90)) +
  geom_point() +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  scale_color_manual('Regular yearly ice over', values = pal, labels = c('Yes', 'No')) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = 'top'); plt.long.lag

# reconstruction of first regression tree branch from Sharma et al. 2019, doi:10.1038/s41558-018-0393-5
plt.temp.branch <-
  ggplot(NCC.data, aes(MeanAnnualAirTemp_c, y, group = lakecode,
                       col = IntermittentIceCover)) +
  geom_vline(xintercept = 8.4, lty = 2) +
  geom_jitter() +
  labs(x = expression(paste('Mean annual air temperature (', degree, 'C)')), y = NULL) +
  scale_color_manual('Yearly ice over', values = pal,
                     labels = c('continuous', 'intermittent')) +
  scale_x_continuous(breaks = c(-10, 0, 8.4, 10), minor_breaks = c(-15, -5, 5),
                     limits = c(min(NCC.data$MeanAnnualAirTemp_c), 12)) +
  scale_y_continuous(breaks = NULL, minor_breaks = NULL) +
  theme(legend.position = 'bottom') +
  
  # add arrow for Mohansic
  annotate(geom = 'curve', x = mohansic$MeanAnnualAirTemp_c, y = mohansic$y + 0.25,
           xend = mohansic$MeanAnnualAirTemp_c - .1, yend = mohansic$y + .02,
           curvature = .3, arrow = arrow(length = unit(2, 'mm'))) +
  annotate(geom = 'text', x = mohansic$MeanAnnualAirTemp_c, y = mohansic$y + 0.3,
           label = 'Mohansic', hjust = 'center', col = pal[1]) +
  
  # add arrow for Lake Clark
  annotate(geom = 'curve', x = -15, y = clark$y,
           xend = clark$MeanAnnualAirTemp_c - .2, yend = clark$y,
           curvature = .3, arrow = arrow(length = unit(2, 'mm'))) +
  annotate(geom = 'text', x = -16, y = clark$y + .05, label = 'Lake Clark',
           hjust = 'center', col = pal[2]); plt.temp.branch
#save.plt(plt.temp.branch,'temp-branch.pdf')

# add marginal density plot
ggMarginal(plt.temp.branch, margins = 'x', fill = 'grey')

# freeze events by Year and DOY after june
ice.sub <- subset(ice, year >= start.year & year %% 5 == 0 & !is.na(froze.bool)) %>%
  mutate(On.DOY.jul.plot = if_else(froze.bool, On.DOY.jul, -1))
plt.dotplot.on <-
  ggplot(ice.sub, aes(On.DOY.jul.plot, year, col = froze.bool, alpha = froze.bool,
                      shape = froze.bool)) +
  geom_point() +
  scale_color_manual('Froze', values = c('TRUE' = 'black', 'FALSE' = 'red'), labels = c('No', 'Yes')) +
  scale_alpha_manual('Froze', values = c('TRUE' = 0.25, 'FALSE' = .5), labels = c('No', 'Yes')) +
  scale_shape_manual('Froze', values = c('TRUE' = 20, 'FALSE' = 17), labels = c('No', 'Yes')) +
  labs(x = june.lab, y = NULL) +
  theme(legend.position = 'left')

plt.dotplot.on.h <- ggMarginal(plt.dotplot.on + theme(legend.position = 'none'),
                               margins = 'y', fill = pal[1]) # density plot of the amount of dots
plt.dotplot.on.h

# thaw events by Year and DOY after june
plt.dotplot.off <- ggplot(ice.sub, aes(Off.DOY.oct, year)) +
  geom_point(alpha = 0.25, shape = 20) +
  labs(x = sept.lab, y = NULL)

plt.dotplot.off.h <- ggMarginal(plt.dotplot.off, margins = 'y', fill = pal[2]) # density plot of the amount of dots
plt.dotplot.off.h

plt.dotplot.on.off <- plot_grid(get_legend(plt.dotplot.on), plt.dotplot.on.h, plt.dotplot.off.h,
                                rel_widths = c(0.3, 1, 1), labels = c(NA, 'a.', 'b.'), nrow = 1)
plt.dotplot.on.off

# freeze events by lake only
plt.dotplot.lake <- ggplot(subset(ice, !is.na(froze.bool)),
                           aes(year, station, color = froze.bool), shape = 19) +
  geom_point(size = .25) +
  labs(x = NULL, y = NULL) +
  scale_y_discrete(breaks = NULL, expand = c(0.025, 0)) +
  xlim(c(start.year, NA)) +
  scale_color_manual('Froze', values = pal[2:1], labels = c('Yes', 'No')) +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  theme(legend.position = 'top'); plt.dotplot.lake

# density of freeze events
plt.dens.freeze <- ggplot(subset(ice, year %in% start.year:2005 & year %% 10 == 0),
                          aes(x = On.DOY.jul, y = factor(year))) +
  geom_density_ridges(fill = '#67a9cf') +
  labs(x = june.lab, y = NULL) +
  scale_y_discrete(expand = c(0, 1.5))
plt.dens.freeze

# density of thaw events
plt.dens.thaw <- ggplot(subset(ice, year %in% start.year:2005 & year %% 10 == 0),
                        aes(x = Off.DOY.oct, y = factor(year))) +
  geom_density_ridges(fill = '#ef8a62') +
  scale_y_discrete(expand = c(0, 1.5)) +
  labs(x = sept.lab, y = NULL)
plt.dens.thaw

# both density plots
plt.dens <- plot_grid(plt.dens.freeze, plt.dens.thaw, labels = c('a.', 'b.')); plt.dens

