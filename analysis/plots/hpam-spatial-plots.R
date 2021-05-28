# data accessing
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
library('sf')        # spatial features
#library('sp')        # spatial objects
library('spData')    # map data for sp package
library('raster')    # to clip rasters to inside map borders
source('functions/save.plt.R') # to save plots easily

# Import data ####
# read in data
ice <- read_rds('data/lake-ice-data.rds') %>%
  filter(year >= 1950) %>%
  mutate(continent = if_else(long > -30, 'Eurasia', 'North America'))
ice.na <- filter(ice, continent == 'North America')
ice.eura <- filter(ice, continent == 'Eurasia')

lakes <- dplyr::filter(ice, ! duplicated(station)) %>%
  dplyr::select(lake, station, long, lat)

# change to PED format
freeze.na <-
  dplyr::select(ice.na, lake, station, year, july.year, froze.bool, On.date, On.DOY, On.DOY.jul, long, lat) %>%
  as_ped(formula = Surv(time = On.DOY.jul,       # follow-up time
                        event = froze.bool) ~ .) # did the lake freeze? TRUE/FALSE
freeze.eura <-
  dplyr::select(ice.eura, lake, station, year, july.year, froze.bool, On.date, On.DOY, On.DOY.jul, long, lat) %>%
  as_ped(formula = Surv(time = On.DOY.jul,       # follow-up time
                        event = froze.bool) ~ .) # did the lake freeze? TRUE/FALSE
thaw.na <-
  dplyr::select(ice.na, lake, station, year, july.year, froze.bool, Off.date, Off.DOY, Off.DOY.oct, long, lat) %>%
  as_ped(formula = Surv(time = Off.DOY.oct,                            # follow-up time
                        event = froze.bool & !is.na(Off.DOY.oct)) ~ .) # lake thawed if it froze & thaw date not NA
thaw.eura <-
  dplyr::select(ice.eura, lake, station, year, july.year, froze.bool, Off.date, Off.DOY, Off.DOY.oct, long, lat) %>%
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

# estimate average freezing date ####
WorldData <-
  map_data('world') %>%
  filter(lat > 30) %>%
  fortify() %>%
  as_tibble()

# import predictions
organize.dates <- function(pred, d) {
  pred <-
    ungroup(pred) %>%
    mutate(stat = case_when(stat == 'lwr' ~ 'Lower 89% CI',
                            stat == 'p' ~ 'Estimate',
                            stat == 'upr' ~ 'Upper 89% CI'),
           stat = factor(stat,levels = c('Lower 89% CI', 'Estimate',
                                         'Upper 89% CI')))
  pred <- dplyr::filter(pred,
                        ! exclude.too.far(pred$long, pred$lat,
                                          lakes$long, lakes$lat,
                                          dist = 0.1))
  pred
}

na.f <-
  read_rds('analysis/plots/predictions/pred-na-f.rds') %>%
  organize.dates(freeze.na)
eura.f <-
  read_rds('analysis/plots/predictions/pred-eura-f.rds') %>%
  organize.dates(freeze.eura)
na.t <-
  read_rds('analysis/plots/predictions/pred-na-t.rds') %>%
  organize.dates(thaw.na)
eura.t <-
  read_rds('analysis/plots/predictions/pred-eura-t.rds') %>%
  organize.dates(thaw.eura)

# check if values are close to 0.5
layout(matrix(1:4, ncol = 2))
hist(na.f$value, breaks = 40)
hist(eura.f$value, breaks = 40)
hist(na.t$value, breaks = 40)
hist(eura.t$value, breaks = 40)
layout(1)

# check if estimated days are extreme
layout(matrix(1:4, ncol = 2))
hist(na.f$tend, breaks = 40)
hist(eura.f$tend, breaks = 40)
hist(na.t$tend, breaks = 40)
hist(eura.t$tend, breaks = 40)
layout(1)

#### Spatio-temporal plots
# change since 1950 ####
diff.dates <- function(d) {
  d %>%
    filter(stat == 'Estimate') %>%
    dplyr::select(long, lat, Year, tend) %>%
    mutate(Year = paste0('y.', Year)) %>%
    pivot_wider(names_from = 'Year', values_from = tend) %>%
    mutate(`1995` = y.1995 - y.1950,
           `2010` = y.2010 - y.1950
#           Significant = if_else(diff < 0,
#                                 # if negative difference signif if upr < lwr
#                                 `y.2010_Upper 89% CI` < `y.1950_Lower 89% CI`,
#                                 # if positive difference signif if lwr > upr
#                                 `y.2010_Lower 89% CI` > `y.1950_Upper 89% CI`)
    ) %>%
    pivot_longer(cols = c('1995', '2010'), names_to = 'year',
                 values_to = 'diff') %>%
    dplyr::select(year, long, lat, diff)
  # pivot_longer(c('Estimate', 'significant'), values_to = 'diff',
  #              names_to = 'stat') %>%
  # mutate(stat = case_when(stat == 'Estimate' ~ 'Estimated difference',
  #                         stat == 'significant' ~ '89% CI Difference'))
}

na.f.diff <- diff.dates(na.f)
eura.f.diff <- diff.dates(eura.f)
na.t.diff <- diff.dates(na.t)
eura.t.diff <- diff.dates(eura.t)

plt.diff <- function(event = c('f', 't'), polar.proj = FALSE, SIGNIF = FALSE) {
  if(event == 'f') {
    eura <- eura.f.diff 
    na <- na.f.diff
    
    DIR <- -1
    
    fill.lab <- 'Change in freeze date\nsince 1950 (days)'
    
    LIMS <- 50 * c(-1, 1)
    
  } else if(event == 't') {
    eura <- eura.t.diff 
    na <- na.t.diff
    
    DIR <- 1
    
    fill.lab <- 'Change in thaw date\nsince 1950 (days)'
    
    LIMS <- 50 * c(-1, 1)
    
  } else  stop('`event` must be "f" or "t"')
  
  if(SIGNIF) {
    na <- filter(na, Significant)
    eura <- filter(eura, Significant)
  }
  
  p <- 
    ggplot(WorldData) +
    facet_grid(. ~ year) +
    geom_map(map = WorldData, aes(group = group, map_id = region),
             fill = 'grey40', color = 'transparent') +
    geom_tile(aes(long, lat, fill = diff), na) +
    geom_tile(aes(long, lat, fill = diff), eura) +
    geom_map(map = WorldData, aes(group = group, map_id = region),
             fill = 'transparent', color = 'grey') +
    labs(x = NULL, y = NULL) +
    scale_x_continuous(breaks = NULL, expand = c(.2, 0)) +
    scale_y_continuous(breaks = NULL, expand = c(.2, 0)) +
    scale_fill_distiller(fill.lab, type = 'div', palette = 5,
                         limits = LIMS, direction = DIR,
                         na.value = 'black') +
    theme(legend.position = 'bottom')
  
  if(polar.proj) p <- p + coord_map('azequidistant', ylim = c(30, 90))
  
  p +
    theme(legend.position = 'bottom', text = element_text(size = 30),
          legend.key.width = unit(2, 'cm'), legend.key.height = unit(1, 'cm'))
}

plt.change.f <- plt.diff('f', polar.proj = TRUE, SIGNIF = FALSE)
plt.change.t <- plt.diff('t', polar.proj = TRUE, SIGNIF = FALSE)

plt.change <- plot_grid(plt.change.f, plt.change.t, labels = c('a.', 'b.'),
                        ncol = 1, label_size = 30)
#save.plt(plt.change, dir = 'hpam-change-in-dates.pdf', width = 17, height = 17)

# histograms of change
ggplot(rbind(na.f.diff, eura.f.diff)) +
  facet_grid(year ~ .) +
  geom_histogram(aes(diff, fill = diff > 0)) +
  scale_fill_brewer(type = 'qual', palette = 6, direction = -1)

ggplot(rbind(na.t.diff, eura.t.diff)) +
  facet_grid(year ~ .) +
  geom_histogram(aes(diff, fill = diff < 0)) +
  scale_fill_brewer(type = 'qual', palette = 6, direction = -1)

# check locations
plt.diff('f') +
  geom_point(aes(long, lat), color = 'green',
             data = dplyr::filter(ice.na, ! duplicated(lake))) +
  geom_point(aes(long, lat), color = 'green',
             data = dplyr::filter(ice.eura, ! duplicated(lake)))

plt.diff('t') +
  facet_wrap(NULL) +
  geom_point(aes(long, lat), color = 'green',
             data = dplyr::filter(ice.na, ! duplicated(lake))) +
  geom_point(aes(long, lat), color = 'green',
             data = dplyr::filter(ice.eura, ! duplicated(lake)))

# estimated dates ####
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
  scale_fill_viridis_c(expression(atop(Mean~freezing~date,
                                       (days~after~June~30^{th}))),
                       option = 'B', direction = 1) +
  theme(legend.position = 'bottom', text = element_text(size = 30),
        legend.key.width = unit(2.5, 'cm'), legend.key.height = unit(1, 'cm'))

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
  theme(legend.position = 'bottom', text = element_text(size = 30),
        legend.key.width = unit(2.5, 'cm'), legend.key.height = unit(1, 'cm'))

#save.plt(plt.f, dir = 'hpam-spatial-freeze.pdf', width = 16, height = 18)
#save.plt(plt.t, dir = 'hpam-spatial-thaw.pdf', width = 16, height = 18)
