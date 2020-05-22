# setup ####
# install packages if necessary:
#devtools::install_github('adibender/pammtools')
#install.packages('package-name')

# data accessing
library('here')      # for easier directory referencing, conflicts with lubridate::here
library('readr')     # to read in files as tibbles

# data processing
library('dplyr')     # makes data wrangling easier
library('tidyr')     # makes data wrangling easier
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

# import predictions
organize.dates <- function(pred, d, dist = 0.15) {
  ungroup(pred) %>%
    mutate(stat = case_when(stat == 'lwr' ~ 'Lower',
                            stat == 'p' ~ 'Estimate',
                            stat == 'upr' ~ 'Upper'),
           stat = factor(stat, levels = c('Lower', 'Estimate', 'Upper'))) %>%
    
    # filter to avoid extrapolating too far
    filter(! exclude.too.far(long, lat, d$long, d$lat, dist))
  
}

na.f <-
  read_rds('analysis/plots/predictions/pam-freeze-na-spatial-predictions.rds') %>%
  organize.dates(freeze.na, 0.15)
eura.f <-
  read_rds('analysis/plots/predictions/pam-freeze-eura-spatial-predictions.rds') %>%
  organize.dates(freeze.eura, 0.15)
na.t <-
  read_rds('analysis/plots/predictions/pam-thaw-na-spatial-predictions.rds') %>%
  organize.dates(thaw.na, 0.15)
eura.t <-
  read_rds('analysis/plots/predictions/pam-thaw-eura-spatial-predictions.rds') %>%
  organize.dates(thaw.eura, 0.15)

# check if values are close to 0.5
layout(matrix(1:4, ncol = 2))
hist(na.f$value)
hist(eura.f$value)
hist(na.t$value)
hist(eura.t$value)
layout(1)

# check if estimated days are extreme
layout(matrix(1:4, ncol = 2))
hist(na.f$tend)
hist(eura.f$tend)
hist(na.t$tend)
hist(eura.t$tend)
layout(1)

# Spatio-temporal plots ####
# freeze dates
plt.f <-
  ggplot(WorldData) +
  facet_grid(stat ~ Year) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'grey40', color = 'transparent') +
  geom_tile(aes(long, lat, fill = tend), na.f) +
  geom_tile(aes(long, lat, fill = tend), eura.f) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'transparent', color = 'grey') +
  coord_map('azequidistant') +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(breaks = NULL, expand = c(.2, 0)) +
  scale_y_continuous(breaks = NULL, expand = c(.2, 0)) +
  scale_fill_viridis_c(expression(atop(Mean~freezing~date,(days~after~June~30^{th}))),
                       option = 'B', direction = -1) +
  theme(legend.position = 'bottom', text = element_text(size = 56),
        legend.key.width = unit(5, 'cm'), legend.key.height = unit(2, 'cm'))

# thaw dates
plt.t <-
  ggplot(WorldData) +
  facet_grid(stat ~ Year) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'grey40', color = 'transparent') +
  geom_tile(aes(long, lat, fill = tend), na.t) +
  geom_tile(aes(long, lat, fill = tend), eura.t) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'transparent', color = 'grey') +
  coord_map('azequidistant') +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(breaks = NULL, expand = c(.1, 0)) +
  scale_y_continuous(breaks = NULL, expand = c(.1, 0)) +
  scale_fill_viridis_c(expression(atop(Mean~thawing~date,(days~after~Sept.~30^{th}))),
                       option = 'B', direction = -1) +
  theme(legend.position = 'bottom', text = element_text(size = 56),
        legend.key.width = unit(5, 'cm'), legend.key.height = unit(2, 'cm'))

#save.plt(plt.f, dir = 'hpam-spatial-freeze.pdf', width = 32, height = 18)
#save.plt(plt.t, dir = 'hpam-spatial-thaw.pdf', width = 32, height = 18)

# change 1950-2010
diff.dates <- function(d) {
  d %>%
    select(long, lat, stat, Year, tend) %>%
    mutate(Year = paste0('y.', Year)) %>%
    pivot_wider(names_from = Year, values_from = tend) %>%
    mutate(diff = y.2010 - y.1950)
}

na.f.diff <- diff.dates(na.f)
eura.f.diff <- diff.dates(eura.f)
na.t.diff <- diff.dates(na.t)
eura.t.diff <- diff.dates(eura.t)

plt.change.f <- 
  ggplot(WorldData) +
  facet_wrap(stat ~ ., ncol = 1) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'grey40', color = 'transparent') +
  geom_tile(aes(long, lat, fill = diff), na.f.diff) +
  geom_tile(aes(long, lat, fill = diff), eura.f.diff) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'transparent', color = 'grey') +
  coord_map('azequidistant') +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(breaks = NULL, expand = c(.2, 0)) +
  scale_y_continuous(breaks = NULL, expand = c(.2, 0)) +
  scale_fill_distiller(expression(atop('Change in expected', 'freeze date (days)')),
                       type = 'div', palette = 5, limits = c(-200, 200)) +
  theme(legend.position = 'right')

plt.change.t <-
  ggplot(WorldData) +
  facet_wrap(stat ~ ., ncol = 1) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'grey40', color = 'transparent') +
  geom_tile(aes(long, lat, fill = diff), na.t.diff) +
  geom_tile(aes(long, lat, fill = diff), eura.t.diff) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'transparent', color = 'grey') +
  coord_map('azequidistant') +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(breaks = NULL, expand = c(.2, 0)) +
  scale_y_continuous(breaks = NULL, expand = c(.2, 0)) +
  scale_fill_distiller(expression(atop('Change in expected', 'thaw date (days)')),
                       type = 'div', palette = 5, limits = c(-110, 110),
                       direction = 1) +
  theme(legend.position = 'right')

plt.diff <- plot_grid(plt.change.f, plt.change.t, labels = c('a.', 'b.'))
#save.plt(plt.diff, dir = 'hpam-change-in-dates.pdf', width = 10, height = 9)
