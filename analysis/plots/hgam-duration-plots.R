# setup ####
# data accessing
library('here')      # for easier directory referencing, conflicts with lubridate::here
library('readr')     # to read in files as tibbles

# data editing
library('dplyr')     # makes data editing easier
library('tibble')    # a tibble is a fancy data.frame

# model fitting
library('mgcv')      # to fit GAMs

# graphics
library('ggplot2')   # fancy plots
library('cowplot')   # ggplot in grids
library('gratia')    # pretty GAM plots
source(here::here('functions/save.plt.R')) # to save plots easily

# palette for plots
pal <- c('#4477AA', '#ff8c00', '#66CCEE', '#009900',
         '#CCBB44', '#EE6677', '#AA3377', '#BBBBBB')
plot(1:length(pal), col = pal, cex = 5, pch = 15)

# for specific years
years <- seq(1950, 2010, by = 10) # years for predictions

# change ggplot theme
theme_set(theme_bw())

# Import data ####
ice <- read_rds(here::here('data/lake-ice-data.rds')) %>%
  filter(Year >= 1950) %>%
  mutate(continent = if_else(long > -30, 'Eurasia', 'North America'))
ice.na <- filter(ice, continent == 'North America')
ice.eura <- filter(ice, continent == 'Eurasia')

# import models ####
# GAMs for duration
gam.dur.na <- read_rds('analysis/models/gam-dur-na-twlss2.rds')
gam.dur.eura <- read_rds('analysis/models/gam-dur-eura-twlss2.rds')

# plots ####
# large unexplained variance in duration for some individual lakes
plt.obs.fit <-
  plot_grid(ggplot(mapping = aes(gam.dur.na$model$duration,
                                 gam.dur.na$fitted.values[, 1])) +
              geom_hex(col = 'grey') +
              geom_segment(aes(x = 0, y = 0, xend = 350, yend = 350), col = 'red') +
              scale_fill_viridis_c('Number of data points', option = 'B', direction = -1) +
              theme(legend.position = 'right') +
              labs(x = 'Observed', y = 'Estimated'),
            ggplot(mapping = aes(gam.dur.eura$model$duration,
                                 gam.dur.eura$fitted.values[, 1])) +
              geom_hex(col = 'grey') +
              geom_segment(aes(x = 0, y = 0, xend = 300, yend = 300), col = 'red') +
              scale_fill_viridis_c('Number of data points', option = 'B', direction = -1) +
              theme(legend.position = 'right') +
              labs(x = 'Observed', y = 'Estimated'),
            labels = c('a.', 'b.'), ncol = 1)
#save.plt(plt.obs.fit, 'gam-dur-fit-obs.pdf', width = 7, height = 10)

# average trends
## north america
# new data for predictions
newd.dur.na <- expand.grid(Year = years,
                           long = seq(min(ice.na$long) - 5,
                                      max(ice.na$long) + 5,
                                      length.out = 150),
                           lat = seq(min(ice.na$lat) - 5,
                                     max(ice.na$lat) + 5,
                                     length.out = 150),
                           lake = 'NA')

# predict w/o random effects
pred.dur.na <- predict(gam.dur.na, newdata = newd.dur.na,
                       terms = c('s(Year)', 's(long,lat)',
                                 'ti(Year,lat,long)',
                                 's.1(Year)', 's.1(long,lat)',
                                 's.2(Year)', 's.2(long,lat)'),
                       se.fit = FALSE)
pred.dur.na <- as_tibble(pred.dur.na)
pred.dur.na <- mutate(pred.dur.na, mu = exp(V1))   # map from link scale onto response scale
pred.dur.na <- bind_cols(pred.dur.na, newd.dur.na) # bind pred + newd

# eurasia
newd.dur.eura <- expand.grid(Year = years,
                             long = seq(min(ice.eura$long) - 5,
                                        max(ice.eura$long) + 5,
                                        length.out = 150),
                             lat = seq(min(ice.eura$lat) - 5,
                                       max(ice.eura$lat) + 5,
                                       length.out = 150),
                             lake = 'NA')
pred.dur.eura <- predict(gam.dur.eura, newdata = newd.dur.eura,
                         terms = c('s(Year)', 's(long,lat)',
                                   'ti(Year,lat,long)',
                                   's.1(Year)', 's.1(long,lat)',
                                   's.2(Year)', 's.2(long,lat)'),
                         se.fit = FALSE)
pred.dur.eura <- as_tibble(pred.dur.eura)
pred.dur.eura <- mutate(pred.dur.eura, mu = exp(V1))
pred.dur.eura <- bind_cols(pred.dur.eura, newd.dur.eura)

# prevent predictions from being too far from the data
pred.dur.na <- mutate(pred.dur.na,
                      mu = if_else(mu < 366, mu, 366),
                      too.far = exclude.too.far(long, lat, ice.na$long, ice.na$lat, .15),
                      mu.clipped = if_else(too.far, NA_real_, mu))
pred.dur.eura <- mutate(pred.dur.eura,
                        mu = if_else(mu < 366, mu, 366),
                        too.far = exclude.too.far(long, lat, ice.eura$long, ice.eura$lat, .15),
                        mu.clipped = if_else(too.far, NA_real_, mu))

# maps
WorldData <-
  map_data('world') %>%
  filter(lat > 30) %>%
  fortify() %>%
  as_tibble()

plt.dur <- ggplot(mapping = aes(long, lat, fill = mu.clipped)) +
  facet_wrap(Year ~ ., ncol = 3) +
  geom_map(aes(long, lat, group = group), WorldData, fill = 'transparent',
               color = 'grey') +
  geom_tile(data = pred.dur.na) +
  geom_tile(data = pred.dur.eura) +
  coord_map('azequidistant', xlim = c(-180, 180), ylim = c(30, 90)) +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  scale_fill_viridis_c('Average duration of ice cover (days)', option = 'B', direction = -1,
                       na.value = 'transparent') +
  theme(legend.position = 'bottom', text = element_text(size = 56),
        legend.key.width = unit(5, 'cm'), legend.key.height = unit(1, 'cm'))
#ggsave('plots/duration.pdf', plt.dur, height = 32, width = 32)

# average change between 1950 and 2010
filter(pred.dur.na, Year %in% c(1950, 2010)) %>%
  mutate(lat = round(lat, -1)) %>%
  group_by(lat, Year) %>%
  summarise(m = mean(mu)) %>%
  tidyr::pivot_wider(names_from = Year, values_from = m) %>%
  transmute(difference = `2010` - `1950`)

filter(pred.dur.na, Year %in% c(1950, 2010)) %>%
  mutate(long = round(long, -1)) %>%
  group_by(long, Year) %>%
  summarise(m = mean(mu)) %>%
  tidyr::pivot_wider(names_from = Year, values_from = m) %>%
  transmute(difference = `2010` - `1950`)

filter(pred.dur.eura, Year %in% c(1950, 2010)) %>%
  mutate(lat = round(lat, -1)) %>%
  group_by(lat, Year) %>%
  summarise(m = mean(mu)) %>%
  tidyr::pivot_wider(names_from = Year, values_from = m) %>%
  transmute(difference = `2010` - `1950`)

filter(pred.dur.eura, Year %in% c(1950, 2010)) %>%
  mutate(long = round(long, -1)) %>%
  group_by(long, Year) %>%
  summarise(m = mean(mu)) %>%
  tidyr::pivot_wider(names_from = Year, values_from = m) %>%
  transmute(difference = `2010` - `1950`)

# this model can be useful to estimate the global average trends
# lake kallavesi
pred.dur.kallavesi <- tibble(Year = seq(1950, 2006, length.out = 400),
                             long = unique(filter(ice.eura, lake == 'LAKE KALLAVESI (4079)')$long),
                             lat = unique(filter(ice.eura, lake == 'LAKE KALLAVESI (4079)')$lat),
                             lake = 'LAKE KALLAVESI (4079)')
pred.dur.kallavesi <- bind_cols(pred.dur.kallavesi,
                                as.data.frame(predict(gam.dur.eura, newdata = pred.dur.kallavesi,
                                                      se.fit = TRUE))) %>%
  mutate(fit = fit.1, se.fit = se.fit.1,
         mu = exp(fit), lwr = exp(fit - 1.96 * se.fit), upr = exp(fit + 1.96 * se.fit))

ggplot() +
  #geom_point(aes(Year, duration), ice.eura, alpha = .025, na.rm = TRUE) +
  geom_point(aes(Year, duration), filter(ice.eura, lake == 'LAKE KALLAVESI (4079)'),
             col = pal[1]) +
  geom_ribbon(aes(Year, ymin = lwr, ymax = upr), pred.dur.kallavesi, alpha = 0.25, fill = pal[1]) +
  geom_line(aes(Year, mu), pred.dur.kallavesi, color = pal[1]) +
  ylab('Duration')

# ... but don't forget that the average is not representative for all lakes!
# lake stechlin
pred.dur.stech <- tibble(Year = seq(1950, 2006, length.out = 400),
                         long = unique(filter(ice.eura, lake == 'STECHLINSEE')$long),
                         lat = unique(filter(ice.eura, lake == 'STECHLINSEE')$lat),
                         lake = 'STECHLINSEE')
pred.dur.stech <- bind_cols(pred.dur.stech,
                            as.data.frame(predict(gam.dur.eura, newdata = pred.dur.stech,
                                                  se.fit = TRUE))) %>%
  mutate(fit = fit.1, se.fit = se.fit.1,
         mu = exp(fit), lwr = exp(fit - 1.96 * se.fit), upr = exp(fit + 1.96 * se.fit))

ggplot() +
  geom_point(aes(Year, duration), ice.eura, alpha = .025, na.rm = TRUE) +
  geom_point(aes(Year, duration), filter(ice.eura, lake == 'STECHLINSEE'), col = pal[1]) +
  geom_ribbon(aes(Year, ymin = lwr, ymax = upr), pred.dur.stech, alpha = 0.25, fill = pal[1]) +
  geom_line(aes(Year, mu), pred.dur.stech, color = pal[1]) +
  ylab('Duration')
