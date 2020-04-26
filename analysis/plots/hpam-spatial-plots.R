# setup ####
# install packages if necessary:
#devtools::install_github('adibender/pammtools')
#install.packages('package-name')

# data accessing
library('here')      # for easier directory referencing, conflicts with lubridate::here
library('readr')     # to read in files as tibbles

# data processing
library('dplyr')     # makes data editing easier
library('tibble')    # a tibble is a fancy data.frame
library('lubridate') # makes working with dates smoother

# model fitting
library('pammtools') # tools for Piecewise-exponential Additive Mixed Models
library('mgcv')      # to fit GAMs

# graphics
library('ggplot2')   # fancy plots
library('cowplot')   # ggplot in grids
library('gratia')    # pretty GAM plots
source(here::here('functions/save.plt.R')) # to save plots easily

# Import data ####
# read in data
ice <- read_rds(here::here('data/lake-ice-data.rds')) %>%
  filter(Year >= 1950) %>%
  mutate(continent = if_else(long > -30, 'Eurasia', 'North America'))
ice.na <- filter(ice, continent == 'North America')
ice.eura <- filter(ice, continent == 'Eurasia')

# change to PED format
freeze.na <-
  select(ice.na, lake, station, Year, july.year, froze.bool, On.date, On.DOY, On.DOY.jul, long, lat) %>%
  as_ped(formula = Surv(time = On.DOY.jul,       # follow-up time
                        event = froze.bool) ~ .) # did the lake freeze? TRUE/FALSE
freeze.eura <-
  select(ice.eura, lake, station, Year, july.year, froze.bool, On.date, On.DOY, On.DOY.jul, long, lat) %>%
  as_ped(formula = Surv(time = On.DOY.jul,       # follow-up time
                        event = froze.bool) ~ .) # did the lake freeze? TRUE/FALSE
thaw.na <-
  select(ice.na, lake, station, Year, july.year, froze.bool, Off.date, Off.DOY, Off.DOY.oct, long, lat) %>%
  as_ped(formula = Surv(time = Off.DOY.oct,                            # follow-up time
                        event = froze.bool & !is.na(Off.DOY.oct)) ~ .) # lake thawed if it froze & thaw date not NA
thaw.eura <-
  select(ice.eura, lake, station, Year, july.year, froze.bool, Off.date, Off.DOY, Off.DOY.oct, long, lat) %>%
  as_ped(formula = Surv(time = Off.DOY.oct,                            # follow-up time
                        event = froze.bool & !is.na(Off.DOY.oct)) ~ .) # lake thawed if it froze & thaw date not NA

# palette for plots
pal <- c('#4477AA', '#ff8c00', '#66CCEE', '#009900',
         '#CCBB44', '#EE6677', '#AA3377', '#BBBBBB')
plot(1:length(pal), col = pal, cex = 5, pch = 15)

# labels for plots
june.lab <- expression(Days~after~June~30^{th})
sept.lab <- expression(Days~after~September~30^{th})

# change ggplot theme
theme_set(theme_bw())

# predictions ####
WorldData <-
  map_data('world') %>%
  filter(lat > 30) %>%
  fortify() %>%
  as_tibble()

# function to clip predictions
clip.pred <- function(pred, d, dist, max.y = 2005) {
  rbind(mutate(filter(pred, Year == 1950),
               too.far = exclude.too.far(filter(pred, Year == 1950)$long,
                                         filter(pred, Year == 1950)$lat,
                                         filter(d, Year == max.y)$long,
                                         filter(d, Year == max.y)$lat,
                                         dist)),
        mutate(filter(pred, Year == max.y),
               too.far = exclude.too.far(filter(pred, Year == max.y)$long,
                                         filter(pred, Year == max.y)$lat,
                                         filter(d, Year == max.y)$long,
                                         filter(d, Year == max.y)$lat,
                                         dist))) %>%
    filter(!too.far)
}

# import predictions 
preds <- read_rds('analysis/plots/predictions/list-of-hpam-predictions.rds')
na.f <- preds$pred.na.f %>% ungroup() %>% clip.pred(freeze.na, .5)
eura.f <- preds$pred.eura.f %>% ungroup() %>% clip.pred(freeze.eura, .5)
na.t <- preds$pred.na.t %>% ungroup() %>% clip.pred(thaw.na, .5)
eura.t <- preds$pred.eura.t %>% ungroup() %>% clip.pred(thaw.eura, .5)
rm(preds)

layout(matrix(1:4, ncol = 2))
hist(na.f$p)
hist(eura.f$p)
hist(na.t$p)
hist(eura.t$p)
layout(1)

# Spatio-temporal plots ####
# freeze dates
plt.f.na <-
  ggplot(WorldData) +
  facet_wrap(Year ~ ., ncol = 1) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'grey40', color = 'transparent') +
  geom_tile(aes(long, lat, fill = tend), na.f) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'transparent', color = 'grey') +
  coord_map('azequidistant') +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(breaks = NULL, expand = c(.2, 0)) +
  scale_y_continuous(breaks = NULL, expand = c(.2, 0)) +
  scale_fill_viridis_c(expression(atop(Mean~freezing~date,(days~after~June~30^{th}))),
                       option = 'B', direction = -1, values = ) +
  theme(legend.position = 'bottom')

plt.f.eura <-
  ggplot(WorldData) +
  facet_wrap(Year ~ ., ncol = 2) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'grey40', color = 'transparent') +
  geom_tile(aes(long, lat, fill = tend), eura.f) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'transparent', color = 'grey') +
  coord_map('azequidistant') +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(breaks = NULL, expand = c(.2, 0)) +
  scale_y_continuous(breaks = NULL, expand = c(.2, 0)) +
  scale_fill_viridis_c(expression(atop(Mean~freezing~date,(days~after~June~30^{th}))),
                       option = 'B', direction = -1) +
  theme(legend.position = 'bottom')

# thaw dates
plt.t.na <-
  ggplot(WorldData) +
  facet_wrap(Year ~ ., ncol = 1) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'grey40', color = 'transparent') +
  geom_tile(aes(long, lat, fill = tend), na.t) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'transparent', color = 'grey') +
  coord_map('azequidistant') +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(breaks = NULL, expand = c(.1, 0)) +
  scale_y_continuous(breaks = NULL, expand = c(.1, 0)) +
  scale_fill_viridis_c(expression(atop(Mean~thawing~date,(days~after~Sept.~30^{th}))),
                       option = 'B', direction = -1, values = ) +
  theme(legend.position = 'bottom')

plt.t.eura <-
  ggplot(WorldData) +
  facet_wrap(Year ~ ., ncol = 2) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'grey40', color = 'transparent') +
  geom_tile(aes(long, lat, fill = tend), eura.t) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'transparent', color = 'grey') +
  coord_map('azequidistant') +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(breaks = NULL, expand = c(.2, 0)) +
  scale_y_continuous(breaks = NULL, expand = c(.2, 0)) +
  scale_fill_viridis_c(expression(atop(Mean~thawing~date,(days~after~Sept.~30^{th}))),
                       option = 'B', direction = -1) +
  theme(legend.position = 'bottom')

plt.na <- plot_grid(plt.f.na, plt.t.na, labels = c('a.', 'b.'), ncol = 1)
plt.eura <- plot_grid(plt.f.eura, plt.t.eura, labels = c('a.', 'b.'), ncol = 1,
                      axis = 'r')

#save.plt(plt.na, dir = 'hpam-great-lakes.pdf', width = 6, height = 8)
#save.plt(plt.eura, dir = 'hpam-scandinavia.pdf', width = 5.5, height = 8)
