# setup ####
# data accessing
library('here')      # for easier directory referencing, conflicts with lubridate::here

# data processing
library('dplyr')     # for easier data wrangling
library('tidyr')     # for easier data wrangling
library('tibble')    # for fancy dataframes

# model fitting
library('pammtools') # tools for Piecewise-exponential Additive Mixed Models
library('mgcv')      # to fit GAMs

# graphics
library('ggplot2')   # for fancy plots
library('cowplot')   # for ggplots in grids
library('gratia')    # for pretty GAM plots
source(here::here('functions/save.plt.R')) # to save plots easily

# set ggplot theme
theme_set(theme_bw())

# import models (only necessary if evaluating smooths) ####
pam.freeze.na <- read_rds(here::here('analysis/models/pam-freeze-na4.rds'))
pam.freeze.eura <- read_rds(here::here('analysis/models/pam-freeze-eura4.rds'))
pam.thaw.na <- read_rds(here::here('analysis/models/pam-thaw-na4.rds'))
pam.thaw.eura <- read_rds(here::here('analysis/models/pam-thaw-eura4.rds'))

# evaluate fs smooths
year.fs <-
  rbind(evaluate_smooth(pam.freeze.na, smooth = 's(Year,lake)', n = 100) %>%
          mutate(continent = 'North America', event = 'Freeze'),
        evaluate_smooth(pam.freeze.eura, smooth = 's(Year,lake)', n = 100) %>%
          mutate(continent = 'Eurasia', event = 'Freeze'),
        evaluate_smooth(pam.thaw.na, smooth = 's(Year,lake)', n = 100) %>%
          mutate(continent = 'North America', event = 'Thaw'),
        evaluate_smooth(pam.thaw.eura, smooth = 's(Year,lake)', n = 100) %>%
          mutate(continent = 'Eurasia', event = 'Thaw'))

tend.fs <-
  rbind(evaluate_smooth(pam.freeze.na, smooth = 's(tend,lake)', n = 100) %>%
          mutate(continent = 'North America', event = 'Freeze'),
        evaluate_smooth(pam.freeze.eura, smooth = 's(tend,lake)', n = 100) %>%
          mutate(continent = 'Eurasia', event = 'Freeze'),
        evaluate_smooth(pam.thaw.na, smooth = 's(tend,lake)', n = 100) %>%
          mutate(continent = 'North America', event = 'Thaw'),
        evaluate_smooth(pam.thaw.eura, smooth = 's(tend,lake)', n = 100) %>%
          mutate(continent = 'Eurasia', event = 'Thaw'))

# save evaluations
#saveRDS(year.fs, 'analysis/plots/predictions/year-fs-evaluation.rds')
#saveRDS(tend.fs, 'analysis/plots/predictions/tend-fs-evaluation.rds')

# read in evaluations
year.fs <- readRDS("analysis/plots/predictions/year-fs-evaluation.rds")
tend.fs <- readRDS("analysis/plots/predictions/tend-fs-evaluation.rds")

# little variation
plt.year.fs <- ggplot(year.fs, aes(Year, est, group = lake)) +
  facet_grid(event ~ continent) +
  geom_line(alpha = 0.2) +
  labs(x = element_blank(), y = expression(log~'['*widehat(lambda)(t)*']'))

# much more variation!
plt.tend.fs <- ggplot(tend.fs, aes(tend, est, group = lake)) +
  facet_grid(event ~ continent) +
  geom_line(alpha = 0.2) +
  labs(x = 'Number of days after reference date',
       y = expression(log~'['*widehat(lambda)(t)*']'))

plt.fs <- plot_grid(plt.year.fs, plt.tend.fs, nrow = 2, labels = c('a.', 'b.'))
save.plt(plt.fs, 'hpam-fs.pdf')
