# setup ####
# data accessing
library('here')      # for easier directory referencing, conflicts with lubridate::here
library('readr')     # to read in files as tibbles

# data processing
library('dplyr')     # for easier data wrangling
library('tidyr')     # for easier data wrangling
library('tibble')    # for fancy dataframes
library('lubridate') # makes working with dates smoother

# model fitting
library('pammtools') # tools for Piecewise-exponential Additive Mixed Models
library('mgcv')      # to fit GAMs

# graphics
library('ggplot2')   # for fancy plots
library('cowplot')   # for ggplots in grids
library('gratia')    # for pretty GAM plots

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

# Import models ####
pam.freeze.na <- read_rds(here::here('analysis/models/pam-freeze-na4.rds'))
pam.freeze.eura <- read_rds(here::here('analysis/models/pam-freeze-eura4.rds'))
pam.thaw.na <- read_rds(here::here('analysis/models/pam-thaw-na4.rds'))
pam.thaw.eura <- read_rds(here::here('analysis/models/pam-thaw-eura4.rds'))

# make new data for predictions ####
years <- seq(1950, 2010, by = 10)

# function to make predictions and save them
pred.save <- function(d, model, file.name = NULL, n.spatial = 30, dist = 0.2,
                      years = years) {
  # make new data
  pred <- make_newdata(d,
                       Year = years,
                       tend = unique(tend),
                       long = seq(min(d$long) - 5,
                                  max(d$long) + 5,
                                  length.out = n.spatial),
                       lat = seq(min(d$lat) - 5,
                                 max(d$lat) + 5,
                                 length.out = n.spatial),
                       lake = c('NA')) %>%
    
    # filter to avoid extrapolating too far
    filter(! exclude.too.far(long, lat, d$long, d$lat, dist))
  
  # group to estimate hazard for every year and at every location
  pred <- group_by(pred, Year, long, lat)
  
  # make predictions
  pred <- add_surv_prob(pred,
                        object = model,
                        ci = TRUE,
                        se_mult = qnorm(0.945), # 89% CIs
                        terms = c('s(tend)', 's(Year)', 's(long,lat)',
                                  'ti(tend,Year)', 'ti(tend,long,lat)',
                                  'ti(Year,long,lat)'))
  
  # filter to days with cumulative probability closest to 0.5 (i.e. expected occurrece date)
  pred <- mutate(pred,
                 p = 1 - surv_prob,
                 lwr = 1 - surv_lower,
                 upr = 1 - surv_upper) %>%
    select(Year, tend, long, lat, p, lwr, upr) %>%
    pivot_longer(c('p', 'lwr', 'upr'), names_to = 'stat', values_to = 'value') %>%
    group_by(Year, long, lat, stat) %>%
    filter(abs(value - .5) == min(abs(value - .5))) %>%
    filter(!duplicated(paste(Year, long, lat)))
  
  # save the predictions
  if(is.null(file.name)) {
    pred
  } else {
    saveRDS(pred, paste0(file.name, '-spatial-predictions.rds'))
  }
}

# predict and save
pred.na.f <- pred.save(freeze.na, pam.freeze.na, 'pam-freeze-na')
pred.eura.f <- pred.save(freeze.eura, pam.freeze.eura, 'pam-freeze-eura')
pred.na.t.1 <- pred.save(thaw.na, pam.thaw.na, years = years[1:3]) # too big for one run
pred.na.t.2 <- pred.save(thaw.na, pam.thaw.na, years = years[4:5])
pred.na.t.3 <- pred.save(thaw.na, pam.thaw.na, years = years[6:7])
pred.na.t <- saveRDS(rbind(pred.na.t.1, pred.na.t.2, pred.na.t.3),
                     'pam-thaw-na-spatial-predictions.rds')
pred.eura.t <- pred.save(thaw.eura, pam.thaw.eura, 'pam-thaw-eura')
