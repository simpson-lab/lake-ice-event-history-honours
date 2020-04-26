# setup ####
# install packages if necessary:
#devtools::install_github('adibender/pammtools')
#install.packages('package-name')

# data accessing
library('here')      # for easier directory referencing, conflicts with lubridate::here
library('readr')     # to read in files as tibbles

# data editing
library('dplyr')     # makes data editing easier
library('tibble')    # a tibble is a fancy data.frame

# model fitting
library('pammtools') # tools for Piecewise-exponential Additive Mixed Models
library('mgcv')      # to fit GAMs
library('brms')      # to fit censored Gamma HGAM
library('survival')  # survival analysis functions

# graphics
library('ggplot2')   # fancy plots
library('cowplot')   # ggplot in grids
library('gratia')    # pretty GAM plots
source(here::here('functions/save.plt.R')) # to save plots easily

# palette for plots
pal <- c('#4477AA', '#ff8c00', '#66CCEE', '#009900',
         '#CCBB44', '#EE6677', '#AA3377', '#BBBBBB')
plot(1:length(pal), col = pal, cex = 5, pch = 15)

# labels for plots
june.lab <- expression(Days~after~June~30^{th})
sept.lab <- expression(Days~after~September~30^{th})

# change ggplot theme
theme_set(theme_bw())

# Import data ####
ice <- read_rds(here::here('data/lake-ice-data.rds')) %>%
  filter(Year >= 1950) %>%
  mutate(continent = if_else(long > -30, 'Eurasia', 'North America'))
ice.na <- filter(ice, continent == 'North America')
ice.eura <- filter(ice, continent == 'Eurasia')

# convert to PED format
freeze.eura <-
  select(ice.eura, lake, station, Year, july.year, froze.bool, On.date, On.DOY, On.DOY.jul,
         long, lat) %>%
  as_ped(formula = Surv(time = On.DOY.jul,       # follow-up time
                        event = froze.bool) ~ .) # did the lake freeze? TRUE/FALSE
freeze.na <-
  select(ice.na, lake, station, Year, july.year, froze.bool, On.date, On.DOY, On.DOY.jul,
         long, lat) %>%
  as_ped(formula = Surv(time = On.DOY.jul,       # follow-up time
                        event = froze.bool) ~ .) # did the lake freeze? TRUE/FALSE
thaw.eura <-
  select(ice.eura, lake, station, Year, july.year, froze.bool, Off.date, Off.DOY,
         Off.DOY.oct, long, lat) %>%
  as_ped(formula = Surv(time = Off.DOY.oct,                            # follow-up time
                        event = froze.bool & !is.na(Off.DOY.oct)) ~ .) # lake thawed if it froze & thaw date not NA
thaw.na <-
  select(ice.na, lake, station, Year, july.year, froze.bool, Off.date, Off.DOY,
         Off.DOY.oct, long, lat) %>%
  as_ped(formula = Surv(time = Off.DOY.oct,                            # follow-up time
                        event = froze.bool & !is.na(Off.DOY.oct)) ~ .) # lake thawed if it froze & thaw date not NA

# import models ####
# PAMs for freeze and thaw
pam.freeze.na <- read_rds('analysis/models/pam-freeze-na4.rds')
pam.freeze.eura <- read_rds('analysis/models/pam-freeze-eura4.rds')
pam.thaw.na <- read_rds('analysis/models/pam-thaw-na4.rds')
pam.thaw.eura <- read_rds('analysis/models/pam-thaw-eura4.rds')

# GAMs for duration
gam.dur.na <- read_rds('analysis/models/gam-dur-na-twlss2.rds')
gam.dur.eura <- read_rds('analysis/models/gam-dur-eura-twlss2.rds')

# change in hazard ####
# North American lake
lake.na <- c('WASCANA LAKE') # need a vector to prevent error
long.na <- filter(ice.na, lake == lake.na, !duplicated(lake))$long
lat.na <- filter(ice.na, lake == lake.na, !duplicated(lake))$lat

# Eurasian lake
lake.eura <- c('LAKE KALLAVESI (4079)') # need a vector to prevent error
long.eura <- filter(ice.eura, lake == lake.eura, !duplicated(lake))$long
lat.eura <- filter(ice.eura, lake == lake.eura, !duplicated(lake))$lat

### freezing
# predictions grouped by year
f.na.tend <- 
  make_newdata(freeze.na,
               Year = mean(Year),
               tend = unique(tend),
               lake = lake.na,
               long = long.na,
               lat = lat.na) %>%
  group_by(Year) %>%
  add_hazard(pam.freeze.na) %>%
  ggplot(aes(tend)) +
  #geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = pal[1]) +
  geom_line(aes(y = hazard)) +
  labs(x = june.lab, y = expression(hat(lambda)(t))); f.na.tend

f.na.year <- 
  make_newdata(freeze.na,
               Year = seq_range(Year, 100),
               lake = lake.na,
               long = long.na,
               lat = lat.na) %>%
  group_by(Year) %>%
  add_hazard(pam.freeze.na)

# change in freeze/thaw dates ####
years <- c(1950, 1965, 1980, 1995, 2010) # years for predictions

# North America (Wascana Lake, Canada) ----
### freezing
# predictions grouped by year
pred.freeze.na.year <- 
  make_newdata(freeze.na,
               tend = unique(tend),
               Year = years,
               lake = lake.na,
               long = long.na,
               lat = lat.na) %>%
  group_by(Year) %>%
  add_hazard(pam.freeze.na) %>%
  #add_cumu_hazard(pam.freeze.na) %>%
  add_surv_prob(pam.freeze.na) %>%
  mutate(p.f = 1 - surv_prob,
         p.f.lwr = 1 - surv_lower,
         p.f.upr = 1 - surv_upper)

### thawing
pred.thaw.na.year <- 
  make_newdata(thaw.na,
               tend = unique(tend),
               Year = years,
               lake = lake.na,
               long = long.na,
               lat = lat.na) %>%
  group_by(Year) %>%
  add_hazard(pam.thaw.na) %>%
  #add_cumu_hazard(pam.thaw.na) %>%
  add_surv_prob(pam.thaw.na) %>%
  mutate(p.t = 1 - surv_prob,
         p.t.lwr = 1 - surv_lower,
         p.t.upr = 1 - surv_upper)

# Eurasia (Lake Kallavesi, Finland) ----
### freezing
# predictions grouped by year
pred.freeze.eura.year <- 
  make_newdata(freeze.eura,
               tend = unique(tend),
               Year = years,
               lake = lake.eura,
               long = long.eura,
               lat = lat.eura) %>%
  group_by(Year) %>%
  add_hazard(pam.freeze.eura) %>%
  #add_cumu_hazard(pam.freeze.eura) %>%
  add_surv_prob(pam.freeze.eura) %>%
  mutate(p.f = 1 - surv_prob,
         p.f.lwr = 1 - surv_lower,
         p.f.upr = 1 - surv_upper)

### thawing
pred.thaw.eura.year <- 
  make_newdata(thaw.eura,
               tend = unique(tend),
               Year = years,
               lake = lake.eura,
               long = long.eura,
               lat = lat.eura) %>%
  group_by(Year) %>%
  add_hazard(pam.thaw.eura) %>%
  #add_cumu_hazard(pam.thaw.eura) %>%
  add_surv_prob(pam.thaw.eura) %>%
  mutate(p.t = 1 - surv_prob,
         p.t.lwr = 1 - surv_lower,
         p.t.upr = 1 - surv_upper)

# plots ####
# freezing ----
pred.freeze.year <- bind_rows(mutate(pred.freeze.na.year, continent = 'North America'),
                              mutate(pred.freeze.eura.year, continent = 'Eurasia'))

freeze <-
  ggplot(pred.freeze.year, aes(tend, p.f, col = factor(Year))) +
  facet_grid(continent ~ .) +
  geom_ribbon(aes(x = tend, ymin = p.f.lwr, ymax = p.f.upr, fill = factor(Year)),
              alpha = 0.1, inherit.aes = FALSE) +
  geom_line() +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal) +
  labs(x = june.lab, y = 'Probability of being frozen') +
  theme(legend.position = 'top')

# thawing ----
pred.thaw.year <- bind_rows(mutate(pred.thaw.na.year, continent = 'North America'),
                            mutate(pred.thaw.eura.year, continent = 'Eurasia'))

thaw <-
  ggplot(pred.thaw.year, aes(tend, p.t, col = factor(Year))) +
  facet_grid(continent ~ .) +
  geom_ribbon(aes(x = tend, ymin = p.t.lwr, ymax = p.t.upr, fill = factor(Year)),
              alpha = 0.1, inherit.aes = FALSE) +
  geom_line() +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal) +
  labs(x = sept.lab, y = 'Probability being of ice-free') +
  theme(legend.position = 'top')

leg.ft <- get_legend(thaw)

plt.ft <- plot_grid(freeze + theme(legend.position = 'none'),
                    thaw + theme(legend.position = 'none'),
                    labels = c('a.', 'b.'), ncol = 2)
plt.ft <- plot_grid(leg.ft, plt.ft, rel_heights = c(0.1, 1), ncol = 1); plt.ft
#save.plt(plt.ft, 'p-freeze-thaw.pdf', height = 3.5)

# duration ----
# large unexplained variance in duration for some individual lakes
plot_grid(ggplot(mapping = aes(gam.dur.na$model$duration,
                               gam.dur.na$fitted.values[, 1])) +
            geom_hex(col = 'grey') +
            geom_segment(aes(x = 0, y = 0, xend = 350, yend = 350), col = 'red') +
            scale_fill_viridis_c(option = 'B', direction = -1) +
            theme(legend.position = 'bottom') +
            labs(x = 'Observed', y = 'Fitted', title = 'North America'),
          ggplot(mapping = aes(gam.dur.eura$model$duration,
                               gam.dur.eura$fitted.values[, 1])) +
            geom_hex(col = 'grey') +
            geom_segment(aes(x = 0, y = 0, xend = 300, yend = 300), col = 'red') +
            scale_fill_viridis_c(option = 'B', direction = -1) +
            theme(legend.position = 'bottom') +
            labs(x = 'Observed', y = 'Fitted', title = 'Eurasia'))
#save.plt(dir = 'gam-dur-fit-obs.pdf', height = 4.5)

# this model can be useful to estimate the global average trends
# lake kallavesi
pred.dur.eura <- tibble(Year = seq(1950, 2006, length.out = 400),
                        long = unique(filter(ice.eura, lake == 'LAKE KALLAVESI (4079)')$long),
                        lat = unique(filter(ice.eura, lake == 'LAKE KALLAVESI (4079)')$lat),
                        lake = 'LAKE KALLAVESI (4079)')
pred.dur.eura <- bind_cols(pred.dur.eura,
                           as.data.frame(predict(gam.dur.eura, newdata = pred.dur.eura,
                                                 terms = c('s(Year)', 's(Year,lake)'),
                                                 se.fit = TRUE))) %>%
  mutate(fit = fit.1, se.fit = se.fit.1,
         mu = exp(fit), lwr = exp(fit - 1.96 * se.fit), upr = exp(fit + 1.96 * se.fit))

ggplot() +
  #geom_point(aes(Year, duration), ice.eura, alpha = .025, na.rm = TRUE) +
  geom_point(aes(Year, duration), filter(ice.eura, lake == 'LAKE KALLAVESI (4079)'),
             col = pal[1]) +
  geom_ribbon(aes(Year, ymin = lwr, ymax = upr), pred.dur.eura, alpha = 0.25, fill = pal[1]) +
  geom_line(aes(Year, mu), pred.dur.eura, color = pal[1]) +
  ylab('Duration')

# ... but don't forget that the average is not representative for all lakes!
# lake stechlin
pred.dur.eura <- tibble(Year = seq(1950, 2006, length.out = 400),
                        long = unique(filter(ice.eura, lake == 'STECHLINSEE')$long),
                        lat = unique(filter(ice.eura, lake == 'STECHLINSEE')$lat),
                        lake = 'STECHLINSEE')
pred.dur.eura <- bind_cols(pred.dur.eura,
                           as.data.frame(predict(gam.dur.eura, newdata = pred.dur.eura,
                                                 terms = c('s(Year)', 's(Year,lake)'),
                                                 se.fit = TRUE))) %>%
  mutate(fit = fit.1, se.fit = se.fit.1,
         mu = exp(fit), lwr = exp(fit - 1.96 * se.fit), upr = exp(fit + 1.96 * se.fit))

ggplot() +
  geom_point(aes(Year, duration), ice.eura, alpha = .025, na.rm = TRUE) +
  geom_point(aes(Year, duration), filter(ice.eura, lake == 'STECHLINSEE'), col = pal[1]) +
  geom_ribbon(aes(Year, ymin = lwr, ymax = upr), pred.dur.eura, alpha = 0.25, fill = pal[1]) +
  geom_line(aes(Year, mu), pred.dur.eura, color = pal[1]) +
  ylab('Duration')
