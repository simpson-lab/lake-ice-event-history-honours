# setup ####
# install if necessary
#devtools::install_github('adibender/pammtools')
#install.packages('package-name')

# data accessing
library('here')      # for easier directory referencing, possible conflicts with lubridate::here
library('readr')     # to read in files as tibbles

# data editing
library('dplyr')     # makes data editing easier
library('tibble')    # a tibble is a fancy data.frame
library('lubridate') # makes dealing with dates smoother

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
source(here::here('functions/post.ref.date.R'))
source(here::here('functions/plot.pam.R'))

# change ggplot theme
theme_set(theme_bw())

# Import data ####
# read in the buffalo pound data
kallavesi <-
  read_rds(here::here('data/lake-ice-data.rds')) %>%
  filter(station == 'JK02')

# Data pre-processing ####
# format to Piece-wise Exponential Data (PED)
kallavesi.freeze <-
  select(kallavesi,
         'Year', 'july.year', 'On.date', 'On.DOY', 'On.DOY.jul', 'froze.bool') %>%
  as_ped(formula = Surv(time = On.DOY.jul,       # follow-up time
                        event = froze.bool) ~ .) # did the lake freeze?

kallavesi.thaw <-
  select(kallavesi,
         'Year', 'july.year', 'Off.date', 'Off.DOY', 'Off.DOY.oct', 'froze.bool') %>%
  as_ped(formula = Surv(time = Off.DOY.oct,                            # follow-up time
                        event = froze.bool & !is.na(Off.DOY.oct)) ~ .) # lake thawed if it froze and thaw date exists

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
# model the hazard of freezing
pam.freeze <- gam(ped_status ~
                    s(tend, bs = 'tp', k = 20) + # smooth effect of DOY ("follow-up dates") on hazard
                    s(Year, bs = 'tp', k = 20) + # smooth effect of year on hazard
                    ti(tend, Year, bs = 'tp', k = 10, np = FALSE), # tensor interaction between tend and Year
                  data = kallavesi.freeze,
                  family = poisson('log'),
                  offset = offset,
                  method = 'REML', # optimization of smoothness parameter via restricted marginal likelihood
                  control = gam.control(nthreads = 3))

# model the hazard of thawing
pam.thaw <- gam(ped_status ~
                  s(tend, bs = 'tp', k = 20) + # smooth effect of DOY ("follow-up dates") on hazard
                  s(Year, bs = 'tp', k = 20) + # smooth effect of year on hazard
                  ti(tend, Year, bs = 'tp', k = 10, np = FALSE), # tensor interaction between tend and Year
                data = kallavesi.thaw,
                family = poisson('log'),
                offset = offset,
                method = 'REML', # optimization of smoothness parameter via restricted marginal likelihood
                control = gam.control(nthreads = 3))

# model the duration of ice cover
gam.duration <- gam(duration ~ s(Year, bs = 'tp', k = 20), # smooth effect of year on duration
                    data = kallavesi,
                    family = Gamma('log'),
                    method = 'REML', # optimization of smoothness parameter via restricted marginal likelihood
                    control = gam.control(nthreads = 3))

# check if k is large enough
layout(matrix(1:4, ncol = 2))
gam.check(pam.freeze)
gam.check(pam.thaw)
gam.check(gam.duration)
layout(1)

# model summaries (aprroximate p-values)
summary(pam.freeze)
summary(pam.thaw)
summary(gam.duration)

# plots ####
### freezing
# predictions
pred.freeze <-
  kallavesi.freeze %>%
  make_newdata(tend = unique(tend), Year = seq_range(Year, n = 3)) %>%
  add_hazard(pam.freeze)

# stepwise hazard
ggplot(pred.freeze, aes(x = tend, group = Year)) +
  geom_stephazard(aes(y = hazard, col = Year)) +
  geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = Year), alpha = 0.25) +
  labs(x = expression(Days~after~June~30^{th}), y = expression(hat(lambda)(t)))

# smooth hazard
pred.freeze %>%
  group_by(Year) %>%
  add_cumu_hazard(pam.freeze) %>%
  ggplot(aes(x = tend, y = cumu_hazard, ymin = cumu_lower, ymax = cumu_upper, group = Year)) +
  geom_hazard(aes(col = Year)) +
  geom_ribbon(aes(fill = Year), alpha = 0.25) +
  labs(x = expression(Days~after~June~30^{th}), y = expression(hat(Lambda)(t)))

# P(freezing)
pred.freeze %>%
  group_by(Year) %>%
  add_surv_prob(pam.freeze) %>%
  mutate(freeze.prob = 1 - surv_prob,
         f.p.lwr = 1 - surv_lower,
         f.p.upr = 1 - surv_upper) %>%
  ggplot(aes(x = tend, y = freeze.prob, ymin = f.p.lwr, ymax = f.p.upr, group = Year)) +
  geom_line(aes(col = Year)) +
  geom_ribbon(aes(fill = Year), alpha = 0.2) +
  labs(x = expression(Days~after~June~30^{th}), y = 'Probability of freezing')

### thawing
pred.thaw <-
  kallavesi.thaw %>%
  make_newdata(tend = unique(tend), Year = seq_range(Year, n = 3)) %>%
  add_hazard(pam.thaw)

# stepwise hazard
ggplot(pred.thaw, aes(x = tend, group = Year)) +
  geom_stephazard(aes(y = hazard, col = Year)) +
  geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = Year), alpha = 0.25) +
  labs(x = expression(Days~after~September~30^{th}), y = expression(hat(lambda)(t)))

# smooth hazard
pred.thaw %>%
  group_by(Year) %>%
  add_cumu_hazard(pam.thaw) %>%
  ggplot(aes(x = tend, y = cumu_hazard, ymin = cumu_lower, ymax = cumu_upper, group = Year)) +
  geom_hazard(aes(col = Year)) +
  geom_ribbon(aes(fill = Year), alpha = 0.25) +
  labs(x = expression(Days~after~September~30^{th}), y = expression(hat(Lambda)(t)))

# P(thawing)
pred.thaw %>%
  group_by(Year) %>%
  add_surv_prob(pam.thaw) %>%
  mutate(thaw.prob = 1 - surv_prob,
         t.p.lwr = 1 - surv_lower,
         t.p.upr = 1 - surv_upper) %>%
  ggplot(aes(x = tend, y = thaw.prob, ymin = t.p.lwr, ymax = t.p.upr, group = Year)) +
  geom_line(aes(col = Year)) +
  geom_ribbon(aes(fill = Year), alpha = 0.2) +
  labs(x = expression(Days~after~September~30^{th}), y = 'Probability of thawing')

### effect of Year
plot_grid(plot.pam.year(pam.freeze, kallavesi.freeze, 'cornflowerblue', 0.05),
          plot.pam.year(pam.thaw, kallavesi.thaw, 'darkorange', 0.05),
          
          #plot.pam.doy(pam.freeze, kallavesi.freeze, 'cornflowerblue', 0.05),
          #plot.pam.doy(pam.thaw, kallavesi.thaw, 'darkorange', 0.05),
          
          draw(pam.freeze, select = 3),
          draw(pam.thaw, select = 3),
          ncol = 2)

### duration
pred.duration <- tibble(Year = seq_min_max(kallavesi$Year, 400))
pred.duration <- cbind(pred.duration,
                       as.data.frame(predict(gam.duration, pred.duration, se.fit = TRUE))) %>%
  mutate(mean = exp(fit),
         lwr = exp(fit + se.fit * qnorm(0.025)),
         upr = exp(fit + se.fit * qnorm(0.975))) %>%
  as_tibble()

ggplot(pred.duration) +
  geom_ribbon(aes(x = Year, ymin = lwr, ymax = upr), alpha = 0.2) + # confidence intervals
  geom_line(aes(Year, mean)) +                                      # mean line
  geom_point(aes(Year, duration), kallavesi, alpha = 0.5) +         # data points
  labs(x = 'Year C.E.', y = 'Duration')
