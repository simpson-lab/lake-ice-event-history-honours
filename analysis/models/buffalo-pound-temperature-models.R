# setup ####
# data accessing
library('here')      # for easier directory referencing
library('readr')     # to read in files as tibbles

# data editing
library('dplyr')     # makes data editing easier
library('tibble')    # a tibble is a fancy data.frame
library('lubridate') # makes dealing with dates smoother

# model fitting
library('survival')  # for Surv function
library('pammtools') # tools for Piecewise-exponential Additive Mixed Models
library('mgcv')      # to fit GAMs

# graphics
library('ggplot2')   # fancy plots
library('cowplot')   # ggplot in grids

# source functions
source(here::here('functions/post.ref.date.R')) # DOY => post Jun 30 / Sep 30
source(here::here('functions/save.plt.R'))      # to save plots easily

# change ggplot theme
theme_set(theme_bw() + theme(legend.position = 'top'))

# palette for plots
pal <- c('#4477AA', '#ff8c00', '#66CCEE', '#009900',
         '#CCBB44', '#EE6677', '#AA3377', '#BBBBBB')

# plot labels
june.lab <- expression(Days~after~June~30^{th})
sept.lab <- expression(Days~after~September~30^{th})
temp.lab <- expression(Air~temperature~(degree*C))

# Import and process data ####
# read in the buffalo pound data
bp.ice <-
  read_csv(here::here('data', 'bp-weekly-data.csv'),
                          guess_max = 1900) %>%
  select(year, week, temp, iceonweek, iceoffweek, iceondoy, iceoffdoy) %>%
  transmute(Year = year,
            water.temp = temp,
            week = week,
            on.week = iceonweek,
            off.week = iceoffweek,
            on.date = as.Date(paste0(Year - 1, '-12-31')) + iceondoy,   # Dec 31 of previous year + DOY
            off.date = as.Date(paste0(Year - 1, '-12-31')) + iceoffdoy, # Dec 31 of previous year + DOY
            froze = 1L,
            date = date(paste0(Year - 1, '-12-31')) + week * 7,
            frozen = if_else(week < iceonweek &
                               week > iceoffweek, 0, 1)) # was BP frozen?
bp.ice <-
  mutate(bp.ice,
         # DOY of ice on post June 30
         on.doy.jul = post.ref.date('on.date', bp.ice, 'freeze'),
         # DOY of ice off post Sept 30
         off.doy.oct = post.ref.date('off.date', bp.ice, 'thaw'),
         # DOY post June 30
         doy.jul = post.ref.date('date', bp.ice, 'freeze'),
         # DOY post Sept 30
         doy.oct = post.ref.date('date', bp.ice, 'thaw'))

# add air temperature data
bp.temp <- read_rds('data/daily-weather/SK_data_dmtemp.rds') %>%
  filter(Location == 'DAVIDSON') %>%
  transmute(air.temp  = Value,
            Year = year(Date),
            date = Date) %>%
  filter(year(date) %in% unique(bp.ice$Year))
bp.temp <- mutate(bp.temp,
                  doy.jul = post.ref.date('date', bp.temp, 'freeze'),
                  doy.oct = post.ref.date('date', bp.temp, 'thaw'))

bp.ice <- filter(bp.ice, Year %in% year(bp.temp$date))

# format to Piecewise Exponential Data (PED) ####
# (https://adibender.github.io/pammtools//articles/data-transformation.html)

################################################################################
############ how to add better air.temp values for the predictions? ############
################################################################################

# freeze dates
bp.freeze <-
  as_ped(data = list(bp.ice %>%
                       filter(Year < 2008) %>%
                       select(on.doy.jul, froze, Year) %>%
                       group_by(Year) %>%
                       slice(1) %>%
                       ungroup(),
                     bp.ice %>%
                       select(on.doy.jul, froze, Year, date) %>%
                       left_join(select(bp.temp, -doy.oct, -Year), by = 'date') %>%
                       select(Year, doy.jul, air.temp)),
         formula = Surv(on.doy.jul, froze) ~ . +
           concurrent(air.temp, tz_var = 'doy.jul'),
         id = 'Year') %>%
  filter(! is.na(air.temp))

# thaw dates
bp.thaw <-
  as_ped(data = list(bp.ice %>%
                       filter(Year < 2008) %>%
                       select(off.doy.oct, froze, Year) %>%
                       group_by(Year) %>%
                       slice(1) %>%
                       ungroup(),
                     bp.ice %>%
                       select(off.doy.oct, Year, date) %>%
                       left_join(select(bp.temp, -doy.jul, -Year), by = 'date') %>%
                       select(Year, doy.oct, air.temp)),
         formula = Surv(off.doy.oct, froze) ~ . + # always thawed (if frozen)
           concurrent(air.temp, tz_var = 'doy.oct'),
         id = 'Year') %>%
  filter(! is.na(air.temp))

## column values:
# tstart    : first day of period
# tend      : end of period
# interval  : period (can be more than one day)
# offset    : log(interval), accounts for hazard within interval
# ped_status: 1/0 = frozen/not freezen
# Year      : calendar year (e.g. 2001)
# air.temp  : air temperature on the day

# Model fitting ####
# strong correlation in the effects of temp and doy
plot_grid(ggplot(bp.temp, aes(yday(date), air.temp, )) +
            geom_point(alpha = 0.1) +
            labs(x = 'Day of year', y = temp.lab, title = 'Full data'),
          plot_grid(ggplot(bp.freeze, aes(tend, air.temp)) +
                      geom_point() +
                      labs(x = june.lab, y = temp.lab, title = 'pam.freeze data'),
                    ggplot(bp.thaw, aes(tend, air.temp)) +
                      geom_point() +
                      labs(x = sept.lab, y = temp.lab, title = 'pam.thaw data'),
                    ncol = 2),
          ncol = 1)

# model hazard of freezing
pam.bp.freeze <- pamm(ped_status ~
                        s(air.temp, bs = 'tp', k = 5) +
                        s(tend, bs = 'tp', k = 5) +
                        s(Year, bs = 'tp', k = 5) +
                        ti(Year, tend, bs = 'tp', k = 5),
                      data = bp.freeze) # method = 'REML' by default

# model hazard of thawing
pam.bp.thaw <- pamm(ped_status ~
                      s(air.temp, bs = 'tp', k = 5) +
                      s(tend, bs = 'tp', k = 5) +
                      s(Year, bs = 'tp', k = 5) +
                      ti(tend, Year, bs = 'tp', k = 5),
                    data = bp.thaw)

# Plots ####
# F(t) = 1 - exp(-Lambda(t)) = 0.999 if Lambda(t) = -log(0.001) = 6.907755
## Freezing ----
# predictions grouped by year
pred.freeze.year <-
  bp.freeze %>%
  make_newdata(tend = unique(tend),
               Year = round(seq_range(Year, n = 5))) %>%
  group_by(Year) %>%
  add_hazard(pam.bp.freeze, se_mult = qnorm(0.945)) %>%
  add_cumu_hazard(pam.bp.freeze, se_mult = qnorm(0.945)) %>%
  add_surv_prob(pam.bp.freeze, se_mult = qnorm(0.945)) %>%
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
bp.freeze %>%
  make_newdata(tend = c(150), # December 1st
               Year = seq_range(Year, n = 100)) %>%
  group_by(Year) %>%
  add_hazard(pam.bp.freeze) %>%
  ggplot(aes(Year, hazard)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = pal[1], alpha = .2) +
  geom_line(col = pal[1]) +
  ylab(Hazard~of~freezing~on~December~1^{st})

# smooth term of air temperature
bp.freeze %>%
  make_newdata(air.temp = seq_range(air.temp, 400)) %>%
  group_by(air.temp) %>%
  add_hazard(pam.bp.freeze) %>%
  ggplot(aes(air.temp, hazard)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = pal[1], alpha = .2) +
  geom_line(col = pal[1]) +
  labs(x = expression(Air~temperaure~(degree*C)), y = 'Hazard of freezing')

## Thawing ----
# predictions grouped by year
pred.thaw.year <-
  bp.thaw %>%
  make_newdata(tend = unique(tend),
               Year = round(seq_range(Year, n = 5))) %>%
  group_by(Year) %>%
  add_hazard(pam.bp.thaw, se_mult = qnorm(0.945)) %>%
  add_cumu_hazard(pam.bp.thaw, se_mult = qnorm(0.945)) %>%
  add_surv_prob(pam.bp.thaw, se_mult = qnorm(0.945)) %>%
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

# stepwise hazard
plt.t.step.year <- ggplot(pred.thaw.year, aes(x = tend, group = Year)) +
  geom_stephazard(aes(y = hazard, col = factor(Year)), lwd = 1) +
  xlim(c(187, 228)) +
  labs(x = sept.lab, y = expression(widehat(lambda)[thaw](t))) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal) +
  ylim(c(0, 3.5)); plt.t.step.year

# piecewise cumulative hazard
plt.t.haz.year <- ggplot(pred.thaw.year, aes(x = tend, group = Year)) +
  geom_hazard(aes(y = cumu_hazard, col = factor(Year)), lwd = 1) +
  xlim(c(187, 228)) +
  labs(x = sept.lab, y = expression(widehat(Lambda)[thaw](t))) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal) +
  ylim(c(0, 10)); plt.t.haz.year

## Group plots ----
plt.freeze <- plot_grid(plt.f.step.year + theme(legend.position = 'none') + xlab(NULL),
                        plt.f.haz.year + theme(legend.position = 'none') + xlab(NULL),
                        plt.f.p.year + theme(legend.position = 'none'),
                        rel_heights = c(1, 1, 1.1),
                        labels = c('a.', 'c.', 'e.'), vjust = 0,
                        ncol = 1)

plt.thaw <- plot_grid(plt.t.step.year + theme(legend.position = 'none') + xlab(NULL),
                      plt.t.haz.year + theme(legend.position = 'none') + xlab(NULL),
                      plt.t.p.year + theme(legend.position = 'none'),
                      rel_heights = c(1, 1, 1.1),
                      labels = c('b.', 'd.', 'f.'), vjust = 0,
                      ncol = 1)

plt.bp <- plot_grid(get_legend(plt.f.p.year),
                    plot_grid(plt.freeze, plt.thaw, ncol = 2),
                    ncol = 1, rel_heights = c(.25, 3), vjust = 0,
                    labels = c(NA, NA, 'g.'), scale = .95)
