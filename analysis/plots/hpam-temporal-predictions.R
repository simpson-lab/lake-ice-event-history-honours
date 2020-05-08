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
source(here::here('functions/save.plt.R')) # to save plots easily
source(here::here('functions/post.ref.date.R')) # for n of days post June/Sept 30th

# palette for plots
pal <- c('#4477AA', '#ff8c00', '#66CCEE', '#009900',
         '#CCBB44', '#EE6677', '#AA3377', '#BBBBBB')
plot(1:length(pal), col = pal, cex = 5, pch = 15)

# labels for plots
june.lab <- expression(Days~after~June~30^{th})
sept.lab <- expression(Days~after~September~30^{th})

# for specific years
years <- c(1950, 1965, 1980, 1995, 2010) # years for predictions

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

# Import models ####
# PAMs for freeze and thaw
pam.freeze.na <- read_rds('analysis/models/pam-freeze-na4.rds')
pam.freeze.eura <- read_rds('analysis/models/pam-freeze-eura4.rds')
pam.thaw.na <- read_rds('analysis/models/pam-thaw-na4.rds')
pam.thaw.eura <- read_rds('analysis/models/pam-thaw-eura4.rds')

# Temporal predictions ####
# North America
pred.freeze.na <- 
  make_newdata(freeze.na,
               tend = unique(tend),
               Year = years) %>%
  group_by(Year) %>%
  add_cumu_hazard(pam.freeze.na,
                  se_mult = qnorm(0.945),
                  terms = c('s(tend)', 's(Year)', 's(tend,Year)')) %>%
  add_surv_prob(pam.freeze.na,
                se_mult = qnorm(0.945),
                terms = c('s(tend)', 's(Year)', 's(tend,Year)')) %>%
  mutate(p = 1 - surv_prob,
         p.lwr = 1 - surv_lower,
         p.upr = 1 - surv_upper) %>%
  select(tstart, tend, Year, cumu_hazard, cumu_lower, cumu_upper, p, p.lwr, p.upr)

pred.thaw.na <- 
  make_newdata(thaw.na,
               tend = unique(tend),
               Year = years) %>%
  group_by(Year) %>%
  add_cumu_hazard(pam.thaw.na,
                  se_mult = qnorm(0.945),
                  terms = c('s(tend)', 's(Year)', 's(tend,Year)')) %>%
  add_surv_prob(pam.thaw.na,
                se_mult = qnorm(0.945),
                terms = c('s(tend)', 's(Year)', 's(tend,Year)')) %>%
  mutate(p = 1 - surv_prob,
         p.lwr = 1 - surv_lower,
         p.upr = 1 - surv_upper) %>%
  select(tstart, tend, Year, cumu_hazard, cumu_lower, cumu_upper, p, p.lwr, p.upr)

# Eurasia
pred.freeze.eura <- 
  make_newdata(freeze.eura,
               tend = unique(tend),
               Year = years) %>%
  group_by(Year) %>%
  add_cumu_hazard(pam.freeze.eura,
                  se_mult = qnorm(0.945),
                  terms = c('s(tend)', 's(Year)', 's(tend,Year)')) %>%
  add_surv_prob(pam.freeze.eura,
                se_mult = qnorm(0.945),
                terms = c('s(tend)', 's(Year)', 's(tend,Year)')) %>%
  mutate(p = 1 - surv_prob,
         p.lwr = 1 - surv_lower,
         p.upr = 1 - surv_upper) %>%
  select(tstart, tend, Year, cumu_hazard, cumu_lower, cumu_upper, p, p.lwr, p.upr)

pred.thaw.eura <- 
  make_newdata(thaw.eura,
               tend = unique(tend),
               Year = years) %>%
  group_by(Year) %>%
  add_cumu_hazard(pam.thaw.eura,
                  se_mult = qnorm(0.945),
                  terms = c('s(tend)', 's(Year)', 's(tend,Year)')) %>%
  add_surv_prob(pam.thaw.eura,
                se_mult = qnorm(0.945),
                terms = c('s(tend)', 's(Year)', 's(tend,Year)')) %>%
  mutate(p = 1 - surv_prob,
         p.lwr = 1 - surv_lower,
         p.upr = 1 - surv_upper) %>%
  select(tstart, tend, Year, cumu_hazard, cumu_lower, cumu_upper, p, p.lwr, p.upr)

### save predictions
saveRDS(pred.freeze.eura, 'analysis/plots/predictions/pred-freeze-eura.rds')
saveRDS(pred.freeze.na, 'analysis/plots/predictions/pred-freeze-na.rds')
saveRDS(pred.thaw.eura, 'analysis/plots/predictions/pred-thaw-eura.rds')
saveRDS(pred.thaw.na, 'analysis/plots/predictions/pred-thaw-na.rds')
