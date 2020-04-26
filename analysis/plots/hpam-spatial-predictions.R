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

# make new data for predictions (fine scale) ####
# find ranges of data
ranges <- expand.grid(event = c('freeze', 'thaw'), cont = c('na', 'eura')) %>%
  mutate(min.long = sapply(paste0(event, '.', cont),
                           function(x) floor(min(filter(get(x), Year == 2005)$long))),
         max.long = sapply(paste0(event, '.', cont),
                           function(x) ceiling(max(filter(get(x), Year == 2005,
                                                          long < 60)$long))),
         min.lat = sapply(paste0(event, '.', cont),
                          function(x) floor(min(filter(get(x), Year == 2005)$lat))),
         max.lat = sapply(paste0(event, '.', cont),
                          function(x) ceiling(max(filter(get(x), Year == 2005)$lat)))) %>%
  mutate(long = paste0('[', min.long, ', ', max.long, ']'),
         lat = paste0('[', min.lat, ', ', max.lat, ']'))

# spatial freezing data
mk.newd <- function(d, coord, dist = 0.15, max.y = 2005) {
  long0 <- c(coord$min.long, coord$max.long)
  lat0 <- c(coord$min.lat, coord$max.lat)
  
  newd <- make_newdata(d,
                       tend = unique(tend),
                       Year = c(1950, max.y),
                       long = seq_range(long0, by = 0.5),
                       lat = seq_range(lat0, by = 0.5))
  
  rbind(mutate(filter(newd, Year == 1950),
               too.far = exclude.too.far(filter(newd, Year == 1950)$long,
                                         filter(newd, Year == 1950)$lat,
                                         filter(d, Year == max.y)$long,
                                         filter(d, Year == max.y)$lat,
                                         dist)),
        mutate(filter(newd, Year == max.y),
               too.far = exclude.too.far(filter(newd, Year == max.y)$long,
                                         filter(newd, Year == max.y)$lat,
                                         filter(d, Year == max.y)$long,
                                         filter(d, Year == max.y)$lat,
                                         dist))) %>%
    filter(!too.far) %>%
    group_by(Year, long, lat)
}

pred <- function(newd, model) {
  add_surv_prob(newd, object = model, ci = TRUE,
                terms = c('s(tend)', 's(Year)', 's(long,lat)',
                          'ti(tend,Year)', 'ti(tend,long,lat)',
                          'ti(Year,long,lat)')) %>%
    mutate(p = 1 - surv_prob) %>%
    select(Year, tend, long, lat, p) %>%
    filter(abs(p - .5) == min(abs(p - .5))) %>%
    filter(!duplicated(paste(Year, long, lat)))
}

newd.na.f <- mk.newd(freeze.na, filter(ranges, cont == 'na', event == 'freeze'), dist = 1)
newd.eura.f <- mk.newd(freeze.eura, filter(ranges, cont == 'eura', event == 'freeze'), dist = 1)
newd.na.t <- mk.newd(thaw.na,  filter(ranges, cont == 'na', event == 'thaw'), dist = 1)
newd.eura.t <- mk.newd(thaw.eura,  filter(ranges, cont == 'eura', event == 'thaw'), dist = 1)

pred.na.f <- pred(newd.na.f, pam.freeze.na)
pred.eura.f <- pred(newd.eura.f, pam.freeze.eura)
pred.na.t <- pred(newd.na.t, pam.thaw.na)
pred.eura.t <- pred(newd.eura.t, pam.thaw.eura)

# save the predictions
saveRDS(list(pred.na.f = pred.na.f,
             pred.eura.f = pred.eura.f,
             pred.na.t = pred.na.t,
             pred.eura.t = pred.eura.t),
        'list-of-hpam-predictions.rds')
