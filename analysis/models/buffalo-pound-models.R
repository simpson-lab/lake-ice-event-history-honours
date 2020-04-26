# setup ####
# install if necessary
#devtools::install_github('adibender/pammtools')

# data accessing
library('here')      # for easier directory referencing, conflicts w/ lubridate::here
library('readr')     # to read in files as tibbles

# data editing
library('dplyr')     # makes data editing easier
library('tibble')    # a tibble is a fancy data.frame
library('lubridate') # makes dealing with dates smoother
library('tidyr')     # easier to make new data for predictions using expand_grid()

# model fitting
library('pammtools') # tools for Piecewise-exponential Additive Mixed Models
library('mgcv')      # to fit GAMs
library('survival')  # survival analysis functions

# graphics
library('ggplot2')   # fancy plots
library('gganimate') # animated ggplot2 plots
library('cowplot')   # ggplot in grids
library('gratia')    # pretty GAM plots

# source functions
source(here::here('functions/post.ref.date.R')) # to change DOY to post Jun 30 / Sep 30
source(here::here('functions/save.plt.R'))      # to save plots easily

# change ggplot theme
theme_set(theme_bw() + theme(legend.position = 'top'))

# palette for plots
pal <- c('#4477AA', '#ff8c00', '#66CCEE', '#009900',
         '#CCBB44', '#EE6677', '#AA3377', '#BBBBBB')
plot(1:length(pal), col = pal, cex = 5, pch = 15)

# plot labels
june.lab <- expression(Days~after~June~30^{th})
sept.lab <- expression(Days~after~September~30^{th})

# Import and process data ####
# read in the buffalo pound data
################### assuming that dataset starts with lake already frozen ################################################
# weekly format
bp.ice.weekly <- read_csv(here::here('data', 'bp-weekly-data.csv'), guess_max = 1900) %>%
  transmute(Year = year,                                                   # rename
            week = week,                                                   # keep week column
            date = date(paste0(Year - 1, '-12-31')) + week * 7,
            frozen = if_else(week < iceonweek & week > iceoffweek, 0, 1), # was BP frozen?
            july.year = paste0(Year - 1, '-', Year))

# add a column for "number of days after June 30"
bp.ice.weekly$doy.jul <- post.ref.date('date', bp.ice.weekly, event = 'freeze')

# yearly format
bp.ice <- read_csv(here::here('data', 'bp-weekly-data.csv'), guess_max = 1900) %>%
  transmute(Year = year,
            on.date = as.Date(paste0(Year - 1, '-12-31')) + iceondoy,     # Dec 31 of previous year + DOY
            off.date = as.Date(paste0(Year - 1, '-12-31')) + iceoffdoy,   # Dec 31 of previous year + DOY
            season = paste0(Year, '-', Year + 1),                         # year starting in July, like academic year
            duration = as.numeric(off.date - lag(on.date)),               # period of ice-cover
            observed = 1L)                                                # always observed ("dead")

# visualize time-to-event nature of the data
bp.viz.jan <- ggplot(filter(bp.ice.weekly, Year %in% 1981:1989), aes(week, frozen)) +
  geom_ribbon(aes(week, ymin = 0, ymax = frozen), alpha = 0.25, fill = pal[1]) +
  geom_line(color = pal[1]) +
  facet_wrap(. ~ Year, ncol = 3, dir = 'h') +
  scale_y_continuous(breaks = 0:1, labels = c('No', 'Yes')) +
  labs(x = 'Week', y = 'Frozen'); bp.viz.jan

# visualize time-to-event nature of the data
bp.viz.jul <- ggplot(filter(bp.ice.weekly, Year %in% 1981:1989),
                     aes(doy.jul, frozen)) +
  geom_ribbon(aes(doy.jul, ymin = 0, ymax = frozen), alpha = 0.25, fill = pal[1]) +
  geom_line(color = pal[1]) +
  facet_wrap(. ~ july.year, ncol = 3, dir = 'h') +
  scale_y_continuous(breaks = 0:1, labels = c('No', 'Yes')) +
  labs(x = expression(Days~after~June~30^{th}), y = 'Frozen'); bp.viz.jul

bp.viz <- plot_grid(bp.viz.jan, bp.viz.jul, ncol = 1, labels = c('a.', 'b.'), hjust = 0)
#save.plt(bp.viz, dir = 'bp-viz.pdf')

# process data to piece-wise exponential data format
bp.ice <- bp.ice %>%
  dplyr::group_by(Year) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

# add new variables
bp.ice <- mutate(bp.ice,
                 on.doy.jul = post.ref.date('on.date', bp.ice, 'freeze'), # DOY of ice on post June 30
                 off.doy.oct = post.ref.date('off.date', bp.ice, 'thaw')) # DOY of ice off post Sept 30

# format to Piece-wise Exponential Data (PED)
bp.freeze <-
  select(bp.ice, -'off.doy.oct', -'off.date', -'duration') %>%
  as_ped(formula = Surv(time = on.doy.jul,     # follow-up time
                        event = observed) ~ .) # censored or not

bp.thaw <-
  select(bp.ice, -'on.doy.jul', -'on.date', -'duration') %>%
  as_ped(formula = Surv(time = off.doy.oct,    # follow-up time
                        event = observed) ~ .) # event

## column values:
# id        : ID for "patient" (i.e. Year)
# tstart    : first day of censored period
# tend      : end of censored period 
# interval  : censored period (can be more than one day)
# offset    : accounts for hazard within interval
# ped_status: 1/0 = dead/alive (1 if the event was observed)
# Year      : calendar year (e.g. 2001)
# season    : year starting in June (e.g. 2000-2001)

# Model fitting ####
# if fit terminates with step failure just re-fit: change seed => change paths of convergence

# model hazard of freezing
pam.bp.freeze <- gam(ped_status ~
                       s(tend, bs = 'cr', k = 10) +       # effect of DOY
                       s(Year, bs = 'cr', k = 10) +       # effect of Year
                       ti(tend, Year, bs = 'cr', k = 5),  # tensor interaction
                     data = bp.freeze,
                     family = poisson('log'),             # response of Y
                     offset = offset,                     # account for censoring
                     method = 'REML', # optimization of smoothness parameter via restricted marginal likelihood
                     control = gam.control(nthreads = 3)) # use multiple cores

# model hazard of thawing
pam.bp.thaw <- gam(ped_status ~
                     s(tend, bs = 'cr', k = 10) +       # smooth effect of DOY ("follow-up dates") on hazard
                     s(Year, bs = 'cr', k = 10) +       # smooth effect of Year on hazard
                     ti(tend, Year, bs = 'cr', k = 5),  # how the effect of tend changes over the years
                   data = bp.thaw,
                   family = poisson('log'),
                   offset = offset,
                   method = 'REML',
                   control = gam.control(nthreads = 3))

# model duration of ice cover
gam.bp.duration <- gam(duration ~ s(Year, bs = 'cr', k = 10), # smooth effect of Year
                       data = bp.ice,
                       family = Gamma('log'),
                       method = 'REML')

# model diagnostics
pam.bp.freeze # make sure edf is not close to k
pam.bp.thaw
layout(matrix(1:4, ncol = 2))
gam.check(gam.bp.duration)
layout(1)

# model summaries (p-values are aprroximate!)
summary(pam.bp.freeze)
summary(pam.bp.thaw)
summary(gam.bp.duration)

# Plots ####
## Freezing ----
# predictions grouped by year
pred.freeze.year <- 
  bp.freeze %>%
  make_newdata(tend = unique(tend),
               Year = round(seq_range(1980:2014, n = 5))) %>%
  group_by(Year) %>%
  add_hazard(pam.bp.freeze) %>%
  add_cumu_hazard(pam.bp.freeze) %>%
  add_surv_prob(pam.bp.freeze) %>%
  mutate(freeze.prob = 1 - surv_prob,
         f.p.lwr = 1 - surv_lower,
         f.p.upr = 1 - surv_upper)

# predictions grouped by tend
pred.freeze.tend <- 
  bp.freeze %>%
  make_newdata(tend = seq_range(tend, n = 5),
               Year = seq_range(Year, by = 1)) %>%
  group_by(tend) %>%
  add_hazard(pam.bp.freeze) %>%
  add_cumu_hazard(pam.bp.freeze) %>%
  add_surv_prob(pam.bp.freeze) %>%
  mutate(freeze.prob = 1 - surv_prob,
         f.p.lwr = 1 - surv_lower,
         f.p.upr = 1 - surv_upper)

# P(frozen) by year
plt.f.p.year <- ggplot(pred.freeze.year, aes(tend, freeze.prob, group = Year)) +
  geom_line(aes(color = factor(Year)), lwd = 1) +
  geom_ribbon(aes(ymin = f.p.lwr, ymax = f.p.upr, fill = factor(Year)), alpha = 0.1) +
  labs(x = june.lab, y = expression(widehat(F)[freeze](t))) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal); plt.f.p.year

# P(frozen) by tend
plt.f.p.tend <- ggplot(pred.freeze.tend, aes(Year, freeze.prob, group = tend)) +
  geom_ribbon(aes(ymin = f.p.lwr, ymax = f.p.upr, fill = factor(tend)), alpha = 0.15) +
  geom_line(aes(color = factor(tend)), lwd = 1) +
  labs(x = 'Year', y = expression(widehat(F)[freeze](t))) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal); plt.f.p.tend

# stepwise hazard
plt.f.step.year <- ggplot(pred.freeze.year, aes(x = tend, group = Year)) +
  geom_stephazard(aes(y = hazard, col = factor(Year)), lwd = 1) +
  xlim(c(120, 154)) +
  labs(x = june.lab, y = expression(widehat(lambda)[freeze](t))) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal); plt.f.step.year

# piecewise cumulative hazard
plt.f.haz.year <- ggplot(pred.freeze.year, aes(x = tend, group = Year)) +
  geom_hazard(aes(y = cumu_hazard, col = factor(Year)), lwd = 1) +
  xlim(c(120, 154)) +
  labs(x = june.lab, y = expression(widehat(Lambda)[freeze](t))) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal); plt.f.haz.year

# smooth term of year
f.year <-
  bp.freeze %>%
  make_newdata(tend = c(150), # December 1st
               Year = seq_range(Year, n = 100)) %>%
  group_by(Year) %>%
  add_hazard(pam.bp.freeze) %>%
  ggplot(aes(Year, hazard)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = pal[1], alpha = .2) +
  geom_line(col = pal[1]) +
  ylab(Hazard~of~freezing~on~December~1^{st}); f.year

## Thawing ----
# predictions grouped by year
pred.thaw.year <- 
  bp.thaw %>%
  make_newdata(tend = unique(tend), 
               Year = round(seq_range(1980:2014, n = 5))) %>%
  group_by(Year) %>%
  add_hazard(pam.bp.thaw) %>%
  add_cumu_hazard(pam.bp.thaw) %>%
  add_surv_prob(pam.bp.thaw) %>%
  mutate(thaw.prob = 1 - surv_prob,
         t.p.lwr = 1 - surv_lower,
         t.p.upr = 1 - surv_upper)

# predictions grouped by tend
pred.thaw.tend <- 
  bp.thaw %>%
  make_newdata(tend = seq_range(tend, n = 5),
               Year = seq_range(Year, by = 1)) %>%
  group_by(tend) %>%
  add_hazard(pam.bp.thaw) %>%
  add_cumu_hazard(pam.bp.thaw) %>%
  add_surv_prob(pam.bp.thaw) %>%
  mutate(thaw.prob = 1 - surv_prob,
         t.p.lwr = 1 - surv_lower,
         t.p.upr = 1 - surv_upper)

# P(ice-free) by year
plt.t.p.year <- ggplot(pred.thaw.year, aes(tend, thaw.prob, group = Year)) +
  geom_line(aes(color = factor(Year)), lwd = 1) +
  geom_ribbon(aes(ymin = t.p.lwr, ymax = t.p.upr, fill = factor(Year)), alpha = 0.1) +
  labs(x = sept.lab, y = expression(widehat(F)[thaw](t))) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal); plt.t.p.year

# P(ice-free) by tend
plt.t.p.tend <- ggplot(pred.thaw.tend, aes(Year, thaw.prob, group = tend)) +
  geom_ribbon(aes(ymin = t.p.lwr, ymax = t.p.upr, fill = factor(tend)), alpha = 0.15) +
  geom_line(aes(color = factor(tend)), lwd = 1) +
  labs(x = 'Year', y = expression(widehat(F)[thaw](t))) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal); plt.t.p.tend

# stepwise hazard
plt.t.step.year <- ggplot(pred.thaw.year, aes(x = tend, group = Year)) +
  geom_stephazard(aes(y = hazard, col = factor(Year)), lwd = 1) +
  xlim(c(187, 228)) +
  labs(x = june.lab, y = expression(widehat(lambda)[thaw](t))) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal); plt.t.step.year

# piecewise cumulative hazard
plt.t.haz.year <- ggplot(pred.thaw.year, aes(x = tend, group = Year)) +
  geom_hazard(aes(y = cumu_hazard, col = factor(Year)), lwd = 1) +
  xlim(c(187, 228)) +
  labs(x = june.lab, y = expression(widehat(Lambda)[thaw](t))) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal); plt.t.haz.year

## Duration ----
pred.dur <- tibble(Year = seq_range(bp.ice$Year, n = 400))
pred.dur <- bind_cols(pred.dur,
                      predict(gam.bp.duration, newdata = pred.dur, se = TRUE)) %>%
  mutate(mu = exp(fit),
         lwr = exp(fit - se.fit * 1.96),
         upr = exp(fit + se.fit * 1.96)) %>%
  as_tibble()

plt.dur <- ggplot() +
  geom_point(aes(Year, duration), bp.ice, alpha = 0.5) +
  geom_line(aes(Year, mu), pred.dur) +
  geom_ribbon(aes(Year, ymin = lwr, ymax = upr), pred.dur, alpha = 0.2) +
  labs(x = NULL, y = 'Days of ice cover'); plt.dur

## Group plots ----
plt.freeze <- plot_grid(plt.f.step.year + theme(legend.position = 'none') + xlab(NULL),
                        plt.f.haz.year + theme(legend.position = 'none') + xlab(NULL),
                        plt.f.p.year + theme(legend.position = 'none'),
                        rel_heights = c(1, 1, 1.1),
                        labels = c('a.', 'c.', 'e.'), vjust = 0,
                        ncol = 1); plt.freeze

plt.thaw <- plot_grid(plt.t.step.year + theme(legend.position = 'none') + xlab(NULL),
                      plt.t.haz.year + theme(legend.position = 'none') + xlab(NULL),
                      plt.t.p.year + theme(legend.position = 'none'),
                      rel_heights = c(1, 1, 1.1),
                      labels = c('b.', 'd.', 'f.'), vjust = 0,
                      ncol = 1); plt.thaw

plt.bp <- plot_grid(get_legend(plt.f.p.year),
                    plot_grid(plt.freeze, plt.thaw, ncol = 2),
                    plt.dur,
                    ncol = 1, rel_heights = c(.25, 3, 1), vjust = 0,
                    labels = c(NA, NA, 'g.'), scale = .95)
#save.plt(plt.bp, dir = 'bp.pdf', height = 10, width = 8)

# more functions from the pammtools package ####
int_info(bp.thaw) %>% slice(1:10)            # internal info on intervals from a ped object
sample_info(bp.thaw)                         # sample means for cont and mode for categorical variables
bp.thaw %>% group_by(Year) %>% sample_info() # same as above but grouped by year (useful if 2+ lakes)
ped_info(bp.thaw)                            # all unique intervals with covariates set to mean or mode
