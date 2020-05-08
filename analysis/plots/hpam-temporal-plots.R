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

# import models ####
# PAMs for freeze and thaw
pam.freeze.na <- read_rds('analysis/models/pam-freeze-na4.rds')
pam.freeze.eura <- read_rds('analysis/models/pam-freeze-eura4.rds')
pam.thaw.na <- read_rds('analysis/models/pam-thaw-na4.rds')
pam.thaw.eura <- read_rds('analysis/models/pam-thaw-eura4.rds')

# import predictions ####
pred.freeze.eura <- read_rds('analysis/plots/predictions/pred-freeze-eura.rds')
pred.freeze.na <- read_rds('analysis/plots/predictions/pred-freeze-na.rds')
pred.thaw.eura <- read_rds('analysis/plots/predictions/pred-thaw-eura.rds')
pred.thaw.na <- read_rds('analysis/plots/predictions/pred-thaw-na.rds')

# plots ####
# freezing ----
pred.freeze.year <- bind_rows(mutate(pred.freeze.na, continent = 'North America'),
                              mutate(pred.freeze.eura, continent = 'Eurasia')) %>%
  filter(tend > 100, tend < 260)

freeze.cumu <-
  ggplot(pred.freeze.year, aes(tend, cumu_hazard, col = factor(Year))) +
  facet_grid(continent ~ .) +
  geom_line() +
  scale_color_manual('Year', values = pal) +
  labs(x = june.lab, y = expression(widehat(Lambda)[freeze](t))) +
  theme(legend.position = 'top')

freeze.p <-
  ggplot(pred.freeze.year, aes(tend, p, col = factor(Year))) +
  facet_grid(continent ~ .) +
  geom_ribbon(aes(x = tend, ymin = p.lwr, ymax = p.upr, fill = factor(Year)),
              alpha = 0.1, inherit.aes = FALSE) +
  geom_line() +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal) +
  labs(x = june.lab, y = expression(widehat(F)[freeze](t))) +
  theme(legend.position = 'top')

# thawing ----
pred.thaw.year <- bind_rows(mutate(pred.thaw.na, continent = 'North America'),
                            mutate(pred.thaw.eura, continent = 'Eurasia')) %>%
  filter(tend > 175, tend < 260)

thaw.cumu <-
  ggplot(pred.thaw.year, aes(tend, cumu_hazard, col = factor(Year))) +
  facet_grid(continent ~ .) +
  geom_line() +
  scale_color_manual('Year', values = pal) +
  labs(x = sept.lab, y = expression(widehat(Lambda)[thaw](t))) +
  theme(legend.position = 'top')

thaw.p <-
  ggplot(pred.thaw.year, aes(tend, p, col = factor(Year))) +
  facet_grid(continent ~ .) +
  geom_ribbon(aes(x = tend, ymin = p.lwr, ymax = p.upr, fill = factor(Year)),
              alpha = 0.1, inherit.aes = FALSE) +
  geom_line() +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal) +
  labs(x = sept.lab, y = expression(widehat(F)[thaw](t))) +
  theme(legend.position = 'top')

leg.ft <- get_legend(thaw.p)

plt.ft <- plot_grid(freeze.cumu + theme(legend.position = 'none'),
                    thaw.cumu + theme(legend.position = 'none'),
                    freeze.p + theme(legend.position = 'none'),
                    thaw.p + theme(legend.position = 'none'),
                    labels = c('a.', 'b.', 'c.', 'd.'), ncol = 2)
plt.ft <- plot_grid(leg.ft, plt.ft, rel_heights = c(0.1, 1), ncol = 1); plt.ft
#save.plt(plt.ft, 'p-freeze-thaw.pdf', height = 6.5)

# factor smooths ####
# import predictions...
pred.fs <- readRDS('analysis/plots/predictions/fs-predictions.rds')

# ... or make new predictions (takes a while)
eval.fs <- function(model, d) {
  evaluate_smooth(model,
                  smooth = 's(Year,lake)',
                  newdata = make_newdata(d,
                                         Year = seq_range(Year, 100),
                                         lake = unique(lake)))
}

pred.fs.freeze.na <- eval.fs(pam.freeze.na, freeze.na)
pred.fs.freeze.eura <- eval.fs(pam.freeze.eura, freeze.eura)
pred.fs.thaw.na <- eval.fs(pam.thaw.na, thaw.na)
pred.fs.thaw.eura <- eval.fs(pam.thaw.eura, thaw.eura)

pred.fs <- bind_rows(mutate(pred.fs.freeze.na,
                            continent = 'North America',
                            event = 'Freeze'),
                     mutate(pred.fs.freeze.eura,
                            continent = 'Eurasia',
                            event = 'Freeze'),
                     mutate(pred.fs.thaw.na,
                            continent = 'North America',
                            event = 'Thaw'),
                     mutate(pred.fs.thaw.eura,
                            continent = 'Eurasia',
                            event = 'Thaw'))
#saveRDS(pred.fs, 'analysis/plots/predictions/fs-predictions.rds')

plt.fs <-
  ggplot(pred.fs, aes(Year, est, group = lake)) +
  facet_grid(event ~ continent) +
  geom_line(alpha = 0.2) +
  ylab('Deviation from the average hazard')
#save.plt(plt.fs, 'fs-plots.pdf')
