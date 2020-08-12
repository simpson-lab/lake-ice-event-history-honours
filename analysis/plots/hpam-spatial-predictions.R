# setup ####
# data accessing
library('here')      # for easier directory referencing
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
library('sp')        # for spatial objects
library('spData')    # for map data
library('sf')        # for simple features objects
library('raster')    # to clip predictions to land only

# Import data ####
# read in data
ice <- read_rds(here::here('data/lake-ice-data.rds')) %>%
  rename(Year = year) %>%
  filter(Year >= 1950) %>%
  mutate(continent = if_else(long > -30, 'Eurasia', 'North America'))
ice.na <- filter(ice, continent == 'North America')
ice.eura <- filter(ice, continent == 'Eurasia')

# change to PED format
freeze.na <-
  dplyr::select(ice.na, lake, station, Year, july.year, froze.bool, On.date, On.DOY, On.DOY.jul, long, lat) %>%
  as_ped(formula = Surv(time = On.DOY.jul,       # follow-up time
                        event = froze.bool) ~ .) # did the lake freeze? TRUE/FALSE
freeze.eura <-
  dplyr::select(ice.eura, lake, station, Year, july.year, froze.bool, On.date, On.DOY, On.DOY.jul, long, lat) %>%
  as_ped(formula = Surv(time = On.DOY.jul,       # follow-up time
                        event = froze.bool) ~ .) # did the lake freeze? TRUE/FALSE
thaw.na <-
  dplyr::select(ice.na, lake, station, Year, july.year, froze.bool, Off.date, Off.DOY, Off.DOY.oct, long, lat) %>%
  as_ped(formula = Surv(time = Off.DOY.oct,                            # follow-up time
                        event = froze.bool & !is.na(Off.DOY.oct)) ~ .) # lake thawed if it froze & thaw date not NA
thaw.eura <-
  dplyr::select(ice.eura, lake, station, Year, july.year, froze.bool, Off.date, Off.DOY, Off.DOY.oct, long, lat) %>%
  as_ped(formula = Surv(time = Off.DOY.oct,                            # follow-up time
                        event = froze.bool & !is.na(Off.DOY.oct)) ~ .) # lake thawed if it froze & thaw date not NA

# Import models ####
pam.freeze.na <- read_rds(here::here('analysis/models/pam-freeze-na4.rds'))
pam.freeze.eura <- read_rds(here::here('analysis/models/pam-freeze-eura4.rds'))
pam.thaw.na <- read_rds(here::here('analysis/models/pam-thaw-na4.rds'))
pam.thaw.eura <- read_rds(here::here('analysis/models/pam-thaw-eura4.rds'))

# make new data for predictions ####
YEARS <- c(1950, 1995, 2010)
YEARS <- seq(1950, 2010, by = 10)

# function to make predictions and save them
pred.save <- function(d, model, years, stepsize = 1, latlong.step.ratio = 0.5,
                      TENDS = NULL) {
  
  # filter lat/long newdata to land only
  longlat <-
    # spatial data
    expand_grid(x = seq(round(min(d$long)) - 5,
                        round(max(d$long)) + 5,
                        by = stepsize),
                # finer steps in y bc of distorsion in polar projection
                y = seq(round(min(d$lat)) - 5,
                        round(max(d$lat)) + 5,
                        by = stepsize * latlong.step.ratio)) %>%
    mutate(z = 0) %>%
    
    # convert to raster
    rasterFromXYZ() %>%
    
    # remove spatial points that are not on land
    mask(world) %>%
    
    # convert back to a tibble
    rasterToPoints() %>%
    as_tibble() %>%
    
    # remove unnecessary column
    dplyr::select(-z) %>%
    
    # add location column
    mutate(loc = paste(x, y))
  
  if(is.null(TENDS)) {
    # predict expected dates
    make_newdata(d,
                 lake = c('NA'), # to ensure that random effect is not included
                 Year = years,
                 tend = unique(tend),
                 long = seq(round(min(d$long)) - 5,
                            round(max(d$long)) + 5,
                            by = stepsize),
                 lat = seq(round(min(d$lat)) - 5,
                           round(max(d$lat)) + 5,
                           by = stepsize)) %>%
      mutate(loc = paste(long, lat)) %>%
      
      # filter to terrestrial locations only
      filter(loc %in% longlat$loc) %>%
      
      # group to estimate hazard for every year and at every location
      group_by(Year, long, lat) %>%
      
      # make predictions
      add_surv_prob(object = model,
                    ci = TRUE,
                    se_mult = qnorm(0.945), # 89% CIs
                    terms = c('s(tend)', 's(Year)', 's(long,lat)',
                              'ti(tend,Year)', 'ti(tend,long,lat)',
                              'ti(Year,long,lat)')) %>%
      
      # filter to days with cumulative probability closest to 0.5 (i.e. E(doy))
      mutate(p = 1 - surv_prob,
             lwr = 1 - surv_lower,
             upr = 1 - surv_upper) %>%
      dplyr::select(Year, tend, long, lat, p, lwr, upr) %>%
      pivot_longer(c('p', 'lwr', 'upr'), names_to = 'stat',
                   values_to = 'value') %>%
      group_by(Year, long, lat, stat) %>%
      filter(abs(value - .5) == min(abs(value - .5))) %>%
      filter(!duplicated(paste(Year, long, lat)))
  } else {
    # predict P for given days
    make_newdata(d,
                 lake = c('NA'), # to ensure that random effect is not included
                 Year = years,
                 tend = TENDS,
                 long = seq(round(min(d$long)) - 5,
                            round(max(d$long)) + 5,
                            by = stepsize),
                 lat = seq(round(min(d$lat)) - 5,
                           round(max(d$lat)) + 5,
                           by = stepsize)) %>%
      mutate(loc = paste(long, lat)) %>%
      
      # filter to terrestrial locations only
      filter(loc %in% longlat$loc) %>%
      
      # group to estimate hazard for every year and at every location
      group_by(Year, long, lat) %>%
      
      # make predictions
      add_surv_prob(object = model,
                    ci = TRUE,
                    se_mult = qnorm(0.945), # 89% CIs
                    terms = c('s(tend)', 's(Year)', 's(long,lat)',
                              'ti(tend,Year)', 'ti(tend,long,lat)',
                              'ti(Year,long,lat)')) %>%
      
      # filter to days with cumulative probability closest to 0.5 (i.e. E(doy))
      mutate(p = 1 - surv_prob,
             lwr = 1 - surv_lower,
             upr = 1 - surv_upper) %>%
      dplyr::select(Year, tend, long, lat, p, lwr, upr) %>%
      pivot_longer(c('p', 'lwr', 'upr'), names_to = 'stat', values_to = 'value')
  }

}

# predict and save (too big for a single run)
fun <- function(d, m) {
  # approximately 10GB per year of data
  purrr::map(YEARS,
             function(y) pred.save(d = d, model = m, y = y)) %>%
    bind_rows()
}

# predict average date (P = 0.5)
pred.na.f <- fun(freeze.na, pam.freeze.na)
saveRDS(pred.na.f, 'analysis/plots/predictions/pred-na-f.rds')

pred.eura.f <- fun(freeze.eura, pam.freeze.eura)
saveRDS(pred.eura.f, 'analysis/plots/predictions/pred-eura-f.rds')

pred.na.t <- fun(thaw.na, pam.thaw.na)
saveRDS(pred.na.t, 'analysis/plots/predictions/pred-na-t.rds')

pred.eura.t <- fun(thaw.eura, pam.thaw.eura)
saveRDS(pred.eura.t, 'analysis/plots/predictions/pred-eura-t.rds')

# predict P on given tends
TENDS <- c(75, 125, 175, 225)
pred.na.f.est <- pred.save(freeze.na, pam.freeze.na, YEARS, TENDS = TENDS)
saveRDS(pred.na.f.est, 'analysis/plots/predictions/pred-na-f-est-p.rds')

pred.eura.f.est <- pred.save(freeze.eura, pam.freeze.eura, YEARS, TENDS = TENDS)
saveRDS(pred.eura.f.est, 'analysis/plots/predictions/pred-eura-f-est-p.rds')

pred.na.t.est <- pred.save(thaw.na, pam.thaw.na, YEARS, TENDS = TENDS)
saveRDS(pred.na.t.est, 'analysis/plots/predictions/pred-na-t-est-p.rds')

pred.eura.t.est <- pred.save(thaw.eura, pam.thaw.eura, YEARS, TENDS = TENDS)
saveRDS(pred.eura.t.est, 'analysis/plots/predictions/pred-eura-t-est-p.rds')
