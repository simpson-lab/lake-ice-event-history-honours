# setup ####
# install if necessary
#devtools::install_github('adibender/pammtools')

# data accessing
library('here')      # for easier directory referencing, possible conflicts with lubridate::here
library('readr')     # to read in files as tibbles

# data editing
library('dplyr')     # makes data editing easier
library('tibble')    # a tibble is a fancy data.frame
library('lubridate') # makes dealing with dates smoother
library('tidyr')     # for easy tidying of data

# model fitting
library('pammtools') # tools for Piecewise-exponential Additive Mixed Models
library('mgcv')      # to fit GAMs
library('survival')  # survival analysis functions
library('segmented') # segmented regression

# graphics
library('ggplot2')   # fancy plots
library('gganimate') # animated ggplot2 plots
library('cowplot')   # ggplot in grids
library('gratia')    # pretty GAM plots
source(here::here('functions/save.plt.R')) # to save plots easily
source(here::here('functions/post.ref.date.R')) # days post reference DOY

# plot labels
june.lab <- expression(Days~after~June~30^{th})
sept.lab <- expression(Days~after~September~30^{th})

# palette for plots
pal <- c('#4477AA', '#ff8c00', '#66CCEE', '#009900',
         '#CCBB44', '#EE6677', '#AA3377', '#BBBBBB')
plot(1:length(pal), col = pal, cex = 5, pch = 15)

# change ggplot theme
theme_set(theme_bw())

# Import data ####
# read in the suwa data
suwa <-
  read_rds(here::here('data/lake-ice-data.rds')) %>%
  filter(lake == 'LAKE SUWA')

# Data pre-processing ####
# format to Piece-wise Exponential Data (PED)

## Freezing
# if the lake froze: and the freeze date is avaliable =>    keep the row
#                  : but the freeze date is NA =>           remove the row in the FREEZE data
# if the lake did NOT freeze =>                             keep the row with DOY = NA
suwa.freeze <-
  select(suwa, 'Year', 'july.year', 'On.date', 'On.DOY', 'On.DOY.jul', 'froze.bool') %>%
  filter(!is.na(froze.bool)) %>%                 # remove rows where freezing event is uncertain
  as_ped(formula = Surv(time = On.DOY.jul,       # follow-up time
                        event = froze.bool) ~ .) # did the lake freeze?

## Thawing
# if the lake did NOT freeze (and thus could not thaw) =>   remove the row from the THAW data
# if the lake froze and thawed: and the day is available => keep the row
#                               but the day is NA =>        remove the row from the THAW data
# if the lake froze but did not thaw =>                     keep the row with DOY as NA
suwa.thaw <-
  filter(suwa, froze.bool) %>% # lake must have frozen before thawing
  select('Year', 'july.year', 'Off.date', 'Off.DOY', 'Off.DOY.oct', 'froze.bool') %>%
  as_ped(formula = Surv(time = Off.DOY.oct,                            # follow-up time
                        event = froze.bool & !is.na(Off.DOY.oct)) ~ .) # event (froze and thaw date not NA)

## column values:
# id          : ID for "patient" (i.e. year)
# tstart      : first day of censored period
# tend        : end of censored period 
# interval    : censored period (can be more than one day)
# offset      : accounts for hazard within interval
# ped_status  : 1/0 = dead/alive (1 if the event was observed)
# Year        : calendar year (e.g. 2001)
# june.year   : year starting in June (e.g. 2000-2001)
# On.date     : freezing data
# Off.date    : thawing date
# On.DOY      : freezing day of year
# Off.DOY     : thawing day of year

# Model fitting ####
# segmented regression
plt.suwa.freeze <-
  ggplot() +
  geom_point(aes(Year, On.DOY.jul), suwa, alpha = .5, na.rm = TRUE) + # days post June 30
  labs(y = expression(Freezing~date~after~July~31^{st})); plt.suwa.freeze

lm.suwa.freeze <- lm(On.DOY.jul ~ Year, data = suwa) # fit a LM
segm.suwa.freeze <- segmented(lm.suwa.freeze, seg.Z = ~Year) # segment the LM

pred.segm.suwa.freeze <-
  predict(segm.suwa.freeze,
          tibble(Year = seq_range(suwa$Year, 400)),
          se.fit = TRUE, interval = 'predict') %>%
  as.data.frame() %>%
  transmute(Year = seq_range(suwa$Year, 400),
            mu = fit.fit,
            lwr = fit.lwr,
            upr = fit.upr)

plt.suwa.freeze +
  geom_ribbon(aes(x = Year, ymin = lwr, ymax = upr), pred.segm.suwa.freeze,
              inherit.aes = FALSE, alpha = 0.25, fill = pal[1]) +
  geom_line(aes(Year, mu), pred.segm.suwa.freeze, col = pal[1], lwd = 1)

summary(segm.suwa.freeze) # p-value of the second segment can't be calculated using ANOVA
davies.test(lm.suwa.freeze, ~ Year, k = 2) # test for change in slope
slope(segm.suwa.freeze) # use confidence intervals for the slope of individual segments

# model hazard of freezing
pam.freeze <- gam(ped_status ~
                    s(tend, bs = 'tp', k = 20) + # smooth effect of DOY ("follow-up dates") on hazard
                    s(Year, bs = 'tp', k = 20) + # smooth effect of year on hazard
                    ti(tend, Year, bs = 'tp', k = 10), # tensor interaction between tend and Year
                  data = suwa.freeze,
                  family = poisson('log'),
                  offset = offset,
                  method = 'REML', # optimization of smoothness parameter via restricted marginal likelihood
                  control = gam.control(nthreads = 3))

# model hazard of thawing (not used in thesis)
pam.thaw <- gam(ped_status ~
                  s(tend, bs = 'tp', k = 10) + # smooth effect of DOY ("follow-up dates") on hazard
                  s(Year, bs = 'tp', k = 20) + # smooth effect of year on hazard
                  ti(tend, Year, bs = 'tp', k = 10), # tensor interaction between tend and Year
                data = suwa.thaw,
                family = poisson('log'),
                offset = offset,
                method = 'REML', # optimization of smoothness parameter via restricted marginal likelihood
                control = gam.control(nthreads = 3))

# duration model is not very useful because of all of the zeros and the few non-zero values
ggplot(suwa, aes(Year, duration)) +
  geom_point(na.rm = TRUE)

# model diagnostics
layout(matrix(1:4, ncol = 2))
gam.check(pam.freeze)
gam.check(pam.thaw)
layout(1)

# model summaries (aprroximate p-values)
summary(pam.freeze)
summary(pam.thaw)

# Plots ####
pred.freeze <-
  suwa.freeze %>%
  make_newdata(tend = unique(tend), Year = seq_range(Year, n = 3)) %>%
  group_by(Year) %>%
  add_hazard(pam.freeze) %>%
  add_cumu_hazard(pam.freeze) %>%
  add_surv_prob(pam.freeze) %>%
  ungroup() %>%
  mutate(freeze.prob = 1 - surv_prob,
         f.p.lwr = 1 - surv_lower,
         f.p.upr = 1 - surv_upper)
  

# to find the average freezing day
pred.freeze.full <-
  suwa.freeze %>%
  make_newdata(tend = unique(tend), Year = seq_range(Year, n = 400)) %>%
  group_by(Year) %>%
  add_hazard(pam.freeze) %>%
  add_cumu_hazard(pam.freeze) %>%
  add_surv_prob(pam.freeze) %>%
  ungroup() %>%
  mutate(p = 1 - surv_prob)

### test #####
# expected freezing date
pred.freeze.mean <-
  filter(pred.freeze.full, round(p, 1) == 0.5)%>%
  group_by(Year, p) %>%
  mutate(tend = mean(tend)) %>%
  ungroup()

# upper CI limit
pred.freeze.upr <-
  filter(pred.freeze.full, round(p, 2) == 0.98) %>%
  group_by(Year, p) %>%
  mutate(tend = mean(tend)) %>%
  ungroup()

# lower CI limit
pred.freeze.lwr <-
  filter(pred.freeze.full, round(p, 2) == 0.03) %>%
  group_by(Year, p) %>%
  mutate(tend = mean(tend)) %>%
  ungroup()

sort(c(mean = nrow(pred.freeze.mean),
       upr = nrow(pred.freeze.upr),
       lwr = nrow(pred.freeze.lwr)), decreasing = TRUE)

pred.freeze.lines <-
  left_join(select(pred.freeze.mean, Year, tend),
            select(pred.freeze.upr, Year, tend), 'Year', suffix = c('.mu', '.upr')) %>%
  left_join(rename(select(pred.freeze.lwr, Year, tend), tend.lwr = tend), 'Year') %>%
  transmute(Year = Year, mu = tend.mu, lwr = tend.lwr, upr = tend.upr)

# cumulative density function and mean (orange)
ggplot() +
  geom_tile(aes(Year, tend, fill = p), pred.freeze.full) +
  geom_line(aes(Year, tend), pred.freeze.mean, color = pal[2], lwd = 1) +
  scale_fill_distiller('Probability of being frozen', type = 'div', palette = 5, direction = 1) +
  ylab(june.lab) +
  theme(legend.position = 'top')

# plot pam and segmented regression
preds <- bind_rows(mutate(pred.segm.suwa.freeze, Model = 'CSR'),
                   mutate(pred.freeze.lines, Model = 'PAM')) %>%
  as_tibble()

plt.csr.pam.freeze <- 
  plt.suwa.freeze +
  
  geom_ribbon(aes(Year, ymin = lwr, ymax = upr, fill = Model), preds, alpha = 0.25, na.rm = TRUE) +
  geom_line(aes(Year, mu, color = Model), filter(preds, Model == 'CSR'), lwd = 1) +
  geom_smooth(aes(Year, mu, color = Model), filter(preds, Model == 'PAM'), method = 'gam',
              formula = y ~ s(x), na.rm = TRUE) +
  scale_color_manual('Model', values = pal) +
  scale_fill_manual('Model', values = pal) +
  theme(legend.position = 'bottom'); plt.csr.pam.freeze
#save.plt(plt.csr.pam.freeze, 'suwa-csr-pam-freeze.pdf', height = 6, width = 6)

# other plots ----
# stepwise hazard
ggplot(pred.freeze, aes(x = tend, group = Year)) +
  geom_stephazard(aes(y = hazard, col = factor(Year))) +
  geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = factor(Year)), alpha = 0.25) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal) +
  labs(x = june.lab, y = expression(hat(lambda)(t)))

# smooth hazard
ggplot(pred.freeze, aes(x = tend, y = cumu_hazard, ymin = cumu_lower,
                        ymax = cumu_upper, group = Year)) +
  geom_hazard(aes(col = factor(Year))) +
  geom_ribbon(aes(fill = factor(Year)), alpha = 0.25) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal) +
  labs(x = june.lab, y = expression(hat(Lambda)(t)))

# P(frozen)
ggplot(pred.freeze, aes(x = tend, y = freeze.prob, ymin = f.p.lwr, ymax = f.p.upr, group = Year)) +
  geom_line(aes(col = factor(Year))) +
  geom_ribbon(aes(fill = factor(Year)), alpha = 0.2) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal) +
  labs(x = june.lab, y = 'Probability of being frozen')

### thawing
pred.thaw <-
  suwa.thaw %>%
  make_newdata(tend = unique(tend), Year = seq_range(Year, n = 3)) %>%
  add_hazard(pam.thaw)

# stepwise hazard
ggplot(pred.thaw, aes(x = tend, group = Year)) +
  geom_stephazard(aes(y = hazard, col = Year)) +
  geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = Year), alpha = 0.25) +
  labs(x = sept.lab, y = expression(hat(lambda)(t)))

# smooth hazard
pred.thaw %>%
  group_by(Year) %>%
  add_cumu_hazard(pam.thaw) %>%
  ggplot(aes(x = tend, y = cumu_hazard, ymin = cumu_lower, ymax = cumu_upper, group = Year)) +
  geom_hazard(aes(col = Year)) +
  geom_ribbon(aes(fill = Year), alpha = 0.25) +
  labs(x = sept.lab, y = expression(hat(Lambda)(t)))

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
  labs(x = sept.lab, y = 'Probability of being ice-free')
