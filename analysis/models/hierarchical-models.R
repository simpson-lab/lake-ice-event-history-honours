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

# change ggplot theme
theme_set(theme_bw())

# Import data ####
# read in data
ice <- read_rds(here::here('data/lake-ice-data.rds')) %>%
  filter(Year >= 1800)

# visualize sampling and location bias ####
plt.obs.dens <- ggplot(ice, aes(Year, lat)) +
  geom_point(size = 1e-4, color = 'white') + # invisible points for marginal density plot
  geom_hex(col = 'grey') +
  theme(legend.position = 'bottom') +
  ylab('Latitude (degrees North)') +
  scale_fill_viridis_c('Number of lakes', option = 'B', direction = -1)
ggExtra::ggMarginal(plt.obs.dens, type = 'histogram', fill = 'grey75')

# filter record to start in 1950
ice <- filter(ice, Year >= 1950)

# Data pre-processing ####
# split data by continent
ice$continent <- if_else(ice$long > -30, 'Eurasia', 'North America')
ice.na <- filter(ice, continent == 'North America')
ice.eura <- filter(ice, continent == 'Eurasia')

# plot locations
WorldData <-
  map_data('world') %>%
  filter(lat > 0) %>%
  fortify() %>%
  as_tibble()

ggplot(WorldData) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'grey75', col = 'black', size = 0.5) +
  geom_point(aes(long, lat, col = continent), filter(ice, !duplicated(station))) +
  coord_map('azequidistant', xlim = c(-180, 180), ylim = c(30, 90)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  labs(x = NULL, y = NULL) +
  scale_color_manual('Continent', values = pal) +
  theme(legend.position = 'top')

# format to Piece-wise Exponential Data (PED) ----
## Freezing
# if the lake froze: and the freeze date is avaliable =>    keep the row
#                  : but the freeze date is NA =>           remove the row in the FREEZE data
# if the lake did NOT freeze =>                             keep the row with DOY = NA
freeze.eura <-
  select(ice.eura, lake, station, Year, july.year, froze.bool, On.date, On.DOY, On.DOY.jul, long, lat) %>%
  as_ped(formula = Surv(time = On.DOY.jul,       # follow-up time
                        event = froze.bool) ~ .) # did the lake freeze? TRUE/FALSE
freeze.na <-
  select(ice.na, lake, station, Year, july.year, froze.bool, On.date, On.DOY, On.DOY.jul, long, lat) %>%
  as_ped(formula = Surv(time = On.DOY.jul,       # follow-up time
                        event = froze.bool) ~ .) # did the lake freeze? TRUE/FALSE

## Thawing
# if the lake did NOT freeze (and thus could not thaw) =>   remove the row from the THAW data
# if the lake froze and thawed: and the day is available => keep the row
#                               but the day is NA =>        remove the row from the THAW data
# if the lake froze but did not thaw =>                     keep the row with DOY as NA
thaw.eura <-
  select(ice.eura, lake, station, Year, july.year, froze.bool, Off.date, Off.DOY, Off.DOY.oct, long, lat) %>%
  as_ped(formula = Surv(time = Off.DOY.oct,                            # follow-up time
                        event = froze.bool & !is.na(Off.DOY.oct)) ~ .) # lake thawed if it froze & thaw date not NA

thaw.na <-
  select(ice.na, lake, station, Year, july.year, froze.bool, Off.date, Off.DOY, Off.DOY.oct, long, lat) %>%
  as_ped(formula = Surv(time = Off.DOY.oct,                            # follow-up time
                        event = froze.bool & !is.na(Off.DOY.oct)) ~ .) # lake thawed if it froze & thaw date not NA

# note that the number of rows goes up drastically!
nrow(ice.na)
nrow(freeze.na)

## column values:
# id          : ID for "patient" (i.e. year)
# tstart      : first day of censored period
# tend        : end of censored period 
# interval    : censored period (can be more than one day)
# offset      : accounts for hazard within interval
# ped_status  : 1/0 = dead/alive (1 if the event was observed)
# Year        : calendar year (e.g. 2001)
# july.year   : year starting in July (e.g. 2000-2001)
# On.date     : freezing data
# Off.date    : thawing date
# On.DOY      : freezing day of year
# Off.DOY     : thawing day of year

# Model fitting ####
# Freeze/Thaw ----
# Models:
# Model 1: s(tend) + s(year), s(lat, long)
# Model 2: Model 1 + ti(tend, year) + ti(tend, lat, long) + ti(year, lat, long)
# Model 3: Model 2 + s(Year, station, bs = 'fs')
# Model 4: Model 3 + s(tend, station, bs = 'fs')
# Model 5: Model 4 + s(station, bs = 're')

# notes:
#       bam() is equivalent to using pamm(engine = 'bam')
#       fitted rates numerically 0 occurred in pam.freeze.na3, pam.thaw.eura3, pam.thaw.eura4

# model the hazard of freezing
pam.freeze.na <- bam(ped_status ~                               # frozen? y/n
                       s(tend, bs = 'cr', k = 10) +             # smooth effect of DOY ("follow-up dates")
                       s(Year, bs = 'cr', k = 10) +             # smooth effect of year
                       s(tend, lake, bs = 'fs', k = 10) +       # deviations from tend smooth
                       s(Year, lake, bs = 'fs', k = 10) +       # deviations from year smooth
                       s(long, lat, bs = 'ds', k = 20) +        # smooth effect of space (splines on sphere)
                       #s(lake, bs = 're') +                    # random effect of lake
                       #s(station, bs = 're') +                 # random effect of station
                       ti(tend, Year, bs = 'cr', k = c(5, 5)) + # change in seasonal trend over the years
                       ti(tend, long, lat, bs = c('cr', 'ds'), d = c(1, 2), k = c(5, 5)) + # change in seasonal trend over space
                       ti(Year, long, lat, bs = c('cr', 'ds'), d = c(1, 2), k = c(5, 5)),  # change of spatial trend over the years
                     data = freeze.na,
                     family = poisson('log'),
                     offset = offset,  # accounts for censoring
                     method = 'fREML', # "fast REML"
                     discrete = TRUE,  # for discretization of covariate values and C code level parallelization
                     nthreads = 1,
                     control = gam.control(trace = FALSE, maxit = 1e5))
#saveRDS(pam.freeze.na, 'analysis/models/pam-freeze-na.rds')

pam.freeze.eura.fs <- readRDS('analysis/models/pam-freeze-eura-fs.rds')
pam.freeze.eura <- bam(ped_status ~                                         # frozen? y/n
                         s(tend, bs = 'cr', k = 10) +                    # smooth effect of DOY ("follow-up dates")
                         s(Year, bs = 'cr', k = 10) +                    # smooth effect of year
                         s(tend, lake, bs = 'fs', k = 10) +       # deviations from tend smooth
                         s(Year, lake, bs = 'fs', k = 10) +       # deviations from year smooth
                         s(long, lat, bs = 'ds', k = 20) +               # smooth effect of space (splines on sphere)
                         #s(lake, bs = 're') +                           # random effect of lake
                         #s(station, bs = 're') +                        # random effect of station
                         ti(tend, Year, bs = 'cr', k = c(5, 5)) +        # change in seasonal trend over the years
                         ti(tend, long, lat, bs = c('cr', 'ds'), d = c(1, 2), k = c(5, 5)) + # change in seasonal trend over space
                         ti(Year, long, lat, bs = c('cr', 'ds'), d = c(1, 2), k = c(5, 5)),  # change of spatial trend over the years
                       data = freeze.eura,
                       family = poisson('log'),
                       offset = offset,
                       method = 'fREML', # "fast REML"
                       discrete = TRUE,  # for discretization of covariate values and C code level parallelization
                       nthreads = 1,
                       control = gam.control(trace = TRUE, maxit = 1e5))
#saveRDS(pam.freeze.eura, 'analysis/models/pam-freeze-eura.rds')

# model the hazard of thawing
pam.thaw.na.fs <- readRDS('analysis/models/pam-thaw-na-fs.rds')
pam.thaw.na <- bam(ped_status ~                                      # frozen? y/n
                     s(tend, bs = 'cr', k = 10) +                    # smooth effect of DOY ("follow-up dates")
                     s(Year, bs = 'cr', k = 10) +                    # smooth effect of year
                     s(tend, lake, bs = 'fs', k = 10) +              # deviations from tend smooth
                     s(Year, lake, bs = 'fs', k = 10) +              # deviations from year smooth
                     s(long, lat, bs = 'ds', k = 20) +               # smooth effect of space (splines on sphere)
                     #s(lake, bs = 're') +                           # random effect of lake
                     #s(station, bs = 're') +                        # random effect of station
                     ti(tend, Year, bs = 'cr', k = c(5, 5)) +        # change in seasonal trend over the years
                     ti(tend, long, lat, bs = c('cr', 'ds'), d = c(1, 2), k = c(5, 5)) + # change in seasonal trend over space
                     ti(Year, long, lat, bs = c('cr', 'ds'), d = c(1, 2), k = c(5, 5)),  # change of spatial trend over the years
                   data = thaw.na,
                   family = poisson('log'),
                   offset = offset,
                   method = 'fREML', # "fast REML"
                   discrete = TRUE,  # for discretization of covariate values and C code level parallelization
                   nthreads = 1,
                   control = gam.control(trace = TRUE))
#saveRDS(pam.thaw.na, 'analysis/models/pam-thaw-na.rds')

pam.thaw.eura.fs <- readRDS('analysis/models/pam-thaw-eura-fs.rds')
pam.thaw.eura <- bam(ped_status ~                               # frozen? y/n
                       s(tend, bs = 'cr', k = 10) +             # smooth effect of DOY ("follow-up dates")
                       s(Year, bs = 'cr', k = 10) +             # smooth effect of year
                       s(tend, lake, bs = 'fs', k = 10) +       # deviations from tend smooth
                       s(Year, lake, bs = 'fs', k = 10) +       # deviations from year smooth
                       s(long, lat, bs = 'ds', k = 20) +        # smooth effect of space (splines on sphere)
                       #s(lake, bs = 're') +                    # random effect of lake
                       #s(station, bs = 're') +                 # random effect of station
                       ti(tend, Year, bs = 'cr', k = c(5, 5)) + # change in seasonal trend over the years
                       ti(tend, long, lat, bs = c('cr', 'ds'), d = c(1, 2), k = c(5, 5)) + # change in seasonal trend over space
                       ti(Year, long, lat, bs = c('cr', 'ds'), d = c(1, 2), k = c(5, 5)),  # change of spatial trend over the years
                     data = thaw.eura,
                     family = poisson('log'),
                     offset = offset,
                     method = 'fREML', # "fast REML"
                     discrete = TRUE,  # for discretization of covariate values and C code level parallelization
                     nthreads = 1,
                     control = gam.control(trace = TRUE))
#saveRDS(pam.thaw.eura, 'analysis/models/pam-thaw-eura.rds')

# duration of ice cover ----
## Tweedie models
gam.dur.na <- readRDS('analysis/models/gam-dur-na-tw.rds')
gam.dur.na <- bam(duration ~
                    s(Year, bs = 'cr', k = 20) +            # smooth effect of year
                    s(Year, lake, bs = 'fs', k = 10) +      # deviations from global smooth
                    s(long, lat, k = 20, bs = 'ds') +       # smooth effect of space
                    ti(Year, lat, long, bs = c('cr', 'ds'),
                       d = 1:2, k = c(10, 10)),             # change in s(year) over space
                  data = ice.na,
                  family = tw(link = 'log'), # because zeros are not allowed for Gamma family
                  method = 'fREML', # fast REML
                  discrete = TRUE,
                  nthreads = 4,
                  control = list(trace = TRUE))
#saveRDS(gam.dur.na, 'analysis/models/gam-dur-na-tw.rds')

gam.dur.eura <- readRDS('analysis/models/gam-dur-eura-tw.rds')
gam.dur.eura <- bam(duration ~
                      s(Year, bs = 'cr', k = 20) +            # smooth effect of year
                      s(Year, lake, bs = 'fs', k = 10) +      # deviations from global smooth
                      s(long, lat, k = 20, bs = 'ds') +       # smooth effect of space
                      ti(Year, lat, long, bs = c('cr', 'ds'),
                         d = 1:2, k = c(10, 10)),             # change in s(year) over space
                    data = ice.eura,
                    family = tw(link = 'log'), # because zeros are not allowed for Gamma family
                    method = 'fREML', # fast REML
                    discrete = TRUE,
                    nthreads = 4,
                    control = list(trace = TRUE))
#saveRDS(gam.dur.eura, 'analysis/models/gam-dur-eura-tw.rds')

# bad qqplots!
qq.tw.na <-
  qq_plot(gam.dur.na, method = 'simulate', n_simulate = 1000, level = 0.89) +
  theme(title = element_blank())

qq.tw.eura <-
  qq_plot(gam.dur.eura, method = 'simulate', n_simulate = 1000, level = 0.89) +
  theme(title = element_blank())

plt.dur.tw.qq <- plot_grid(NULL, qq.tw.na, NULL, qq.tw.eura, labels = c('a.', NA, 'b.', NA),
                           ncol = 2, rel_widths = c(0.075, 1))
#save.plt(plt.dur.tw.qq, 'tweedie-qqplots.pdf', width = 4.5)

## fit location-scale Tweedie models (can't use bam())
# twlss1: list(duration ~ # formula for the mean
#                        s(Year, bs = 'cr', k = 20) +              # smooth effect of year
#                          s(Year, lake, bs = 'fs', k = 10) +      # deviations from global smooth
#                          s(long, lat, k = 20, bs = 'ds') +       # smooth effect of space
#                          ti(Year, lat, long, bs = c('cr', 'ds'),
#                             d = 1:2, k = c(10, 10)),             # change in s(year) over space
#                        ~ # formula for the power
#                          s(lake, bs = 're'),
#                        ~ # formula for the scale
#                          s(lake, bs = 're'))
# twlss2:twlss1 + s(Year) + s(long, lat) for power and scale

gam.dur.na <- gam(list(duration ~
                         # formula for the mean
                         s(Year, bs = 'cr', k = 20) +            # smooth effect of year
                         s(Year, lake, bs = 'fs', k = 10) +      # deviations from global smooth
                         s(long, lat, k = 20, bs = 'ds') +       # smooth effect of space
                         ti(Year, lat, long, bs = c('cr', 'ds'),
                            d = 1:2, k = c(10, 10)),             # change in s(year) over space
                       
                       ~ # formula for the power
                         s(Year, bs = 'cr', k = 5) +
                         s(long, lat, k = 5, bs = 'ds') +
                         s(lake, bs = 're'),
                       
                       ~ # formula for the scale
                         s(Year, bs = 'cr', k = 5) +
                         s(long, lat, k = 5, bs = 'ds') +
                         s(lake, bs = 're')),
                  data = ice.na,
                  family = twlss(), # location scale Tweedie
                  method = 'REML',
                  control = gam.control(nthreads = 1, trace = FALSE, maxit = 1000))

# Eurasia
gam.dur.eura <- gam(list(duration ~
                           # formula for the mean
                           s(Year, bs = 'cr', k = 20) +            # smooth effect of year
                           s(Year, lake, bs = 'fs', k = 10) +      # deviations from global smooth
                           s(long, lat, k = 20, bs = 'ds') +       # smooth effect of space
                           ti(Year, lat, long, bs = c('cr', 'ds'),
                              d = 1:2, k = c(10, 10)),             # change in s(year) over space
                         
                         ~ # formula for the power
                           s(Year, bs = 'cr', k = 5) +
                           s(long, lat, k = 5, bs = 'ds') +
                           s(lake, bs = 're'),
                         
                         ~ # formula for the scale
                           s(Year, bs = 'cr', k = 5) +
                           s(long, lat, k = 5, bs = 'ds') +
                           s(lake, bs = 're')),
                    data = ice.eura,
                    family = twlss(), # location scale Tweedie
                    method = 'REML',
                    control = gam.control(nthreads = 1, trace = TRUE))
plot(gam.dur.eura, pages = 1, scale = 0, scheme = 3, main = 'Eurasia duration twlss')

## fit a hurdle Gamma model:
# linear change in hurdle over the years 
hurdle.dur.eura <- brm(bf(duration ~ t2(Year, lat, long, bs = c('cr', 'ds'),
                                        d = c(1, 2), k = c(20, 40)),
                          #s(lake, station, bs = 're'),
                          hu ~ Year),
                       family = hurdle_gamma(), # Gamma given that the lake froze
                       data = ice.eura,
                       chains = 4,
                       iter = 10000,
                       cores = 4,
                       control = list(adapt_delta = .999, max_treedepth = 20))
#saveRDS(hurdle.dur.eura, 'analysis/models/hurdle-dur-eura.rds')

hurdle.dur.na <- brm(bf(duration ~ t2(Year, lat, long, bs = c('cr', 'ds'),
                                      d = c(1, 2), k = c(20, 40)),
                        #s(lake, station, bs = 're'),
                        hu ~ Year),
                     family = hurdle_gamma(), # Gamma given that the lake froze
                     data = ice.na,
                     chains = 4,
                     iter = 10000,
                     cores = 4,
                     control = list(adapt_delta = .999, max_treedepth = 20))


# smooth spatio-temporal change in hurdle
hurdle.dur.na <- brm(bf(duration ~ t2(Year, lat, long, bs = c('cr', 'ds'),
                                      d = c(1, 2), k = c(20, 40)) +
                          #s(lake, station, bs = 're'),
                          hu ~ t2(Year, lat, long, bs = c('cr', 'ds'),
                                  d = c(1, 2), k = c(10, 20))),
                     family = hurdle_gamma(), # Gamma given that the lake froze
                     data = ice.na,
                     chains = 4,
                     iter = 10000,
                     cores = 4,
                     control = list(adapt_delta = .999))

