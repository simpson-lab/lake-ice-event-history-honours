# setup ####
# data accessing
library('readr')     # to read in files as tibbles

# data wrangling
library('dplyr')     # makes data wrangling easier
library('lubridate') # makes dealing with dates smoother

# model fitting
library('pammtools') # tools for Piecewise-exponential Additive Mixed Models
library('mgcv')      # to fit GAMs

# graphics
library('ggplot2')   # fancy plots
library("patchwork") # ggplot in grids
library("gridExtra") # for tables in ggplots
library("gratia")

# source functions
source('functions/post.ref.date.R') # DOY => post Jun 30 / Sep 30
source('functions/save.plt.R')      # to save plots easily

# change ggplot theme
theme_set(theme_bw() + theme(legend.position = 'top'))

# palette for plots
pal <- c('#4477AA', '#ff8c00', '#66CCEE', '#009900',
         '#CCBB44', '#EE6677', '#AA3377', '#BBBBBB')
plot(1:length(pal), col = pal, cex = 5, pch = 15)

# plot labels
june.lab <- expression(Days~after~June~30^{th})

# Import and process data ####
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-ntl.33.37
mendota <-
  read_csv('data/ntl33_v10.csv', col_types = 'ccdcdcdd') %>%
  filter(lakeid == 'ME', year4 >= 1950) %>%
  transmute(year = year4,
            july.year = season,
            on.date = date(ice_on),
            off.date = date(ice_off),
            on.doy = yday(on.date),
            off.doy = yday(off.date),
            froze = !is.na(on.doy),
            thawed = !is.na(on.doy) & !is.na(off.doy),
            duration = ice_duration)
mendota <- mutate(mendota,
                  on.doy.jul = post.ref.date('on.date', mendota, 'freeze'),
                  off.doy.oct = post.ref.date('off.date', mendota, 'thaw'))

# lake mendota always froze after 1950
filter(mendota, is.na(froze + thawed) | froze + thawed != 2)

freeze <-
  select(mendota, year, july.year, froze, on.date, on.doy, on.doy.jul) %>%
  as_ped(formula = Surv(time = on.doy.jul,       # follow-up time
    event = froze) ~ ., # did the lake freeze? TRUE/FALSE
  cut = 150:210) # split into 1-day intervals

## column values:
#' id        : ID for "patient" (i.e. year)
#' tstart    : first day of period
#' tend      : end of period
#' interval  : period (can be more than one day)
#' offset    : accounts for increase in hazard in intervals longer than one day
#' ped_status: 1 if the event was observed, 0 otherwise
#' year      : calendar year (e.g. 2001)
#' season    : year starting in June (e.g. 2000-2001)

# basic scatterplot
mendota_ts_plt <-
  ggplot(mendota) +
  geom_line(aes(year, on.doy.jul), alpha = 0.25) +
  geom_point(aes(year, on.doy.jul), alpha = 1) +
  labs(x = "Year", y = june.lab)
mendota_ts_plt

# visualize time-to-event nature of the data
viz.f <-
  freeze %>%
  bind_rows(tibble(year = unique(mendota$year), tend = 150, ped_status = 0)) %>%
  bind_rows(tibble(year = unique(mendota$year), tend = 210, ped_status = 1)) %>%
  mutate(status = if_else(ped_status == 1, "thawed", "frozen")) %>%
  ggplot() +
  geom_line(aes(tend, year, color = status, group = year)) +
  geom_point(aes(tend, year, color = status, group = year), alpha = 0.5) +
  scale_color_manual("Froze", values = pal[2:1], labels = c("No", "Yes")) +
  labs(x = june.lab, y = "Year") +
  theme(legend.position = "none")
viz.f

# Model fitting ####
#' by default, `pamm()` uses:
#' family = poisson(link = 'log') for the likelihood of Y
#' offset = offset to account for change in hazard in  periods longer than 1 day
#' method = 'REML' for smoothness parameter optimization via REML

# model hazard of freezing
pam.f <-
  pamm(ped_status ~
         s(tend, bs = 'tp', k = 10) +      # effect of DOY
         s(year, bs = 'tp', k = 10) +      # effect of year
         ti(tend, year, bs = 'tp', k = 5), # change in effect of DOY over years
       data = freeze)

# Plots ####
## Freezing ----
# predictions grouped by year
# GLS: ignore warnings; bug in pammtools, see https://github.com/adibender/pammtools/issues/235
pred.freeze <-
  freeze %>%
  make_newdata(tend = unique(tend),
               year = unique(year)) %>%
  group_by(year) %>% # calculate values for each year
  add_hazard(pam.f, se_mult = qnorm(0.945)) %>%
  add_cumu_hazard(pam.f, se_mult = qnorm(0.945)) %>%
  add_surv_prob(pam.f, se_mult = qnorm(0.945)) %>%
  mutate(freeze.prob = 1 - surv_prob,
         f.p.lwr = 1 - surv_lower,
         f.p.upr = 1 - surv_upper)

# tile plot of P(F) by year and ti
ggplot(pred.freeze, aes(year, tend, fill = freeze.prob)) +
  geom_tile() +
  scale_fill_distiller(expression(widehat(F)[freeze](t)), palette = 1,
                       direction = 1) +
  scale_x_continuous(NULL, expand = c(0, 0)) +
  scale_y_continuous(june.lab, expand = c(0, 0))

# P(frozen) by year
plt.f.p.year <-
  filter(pred.freeze, year %in% c(1950, 1980, 2010, 2020)) %>%
  mutate(year = factor(year)) %>%
  ggplot(aes(tend, freeze.prob, group = year)) +
  geom_line(aes(color = year), lwd = 1) +
  geom_ribbon(aes(ymin = f.p.lwr, ymax = f.p.upr, fill = year), alpha = 0.1) +
  labs(x = june.lab, y = expression(widehat(F)[freeze](t))) +
  scale_color_manual("Year", values = pal[c(1:3, 7)],
    aesthetics = c("color", "fill")) +
  theme(legend.key.width = unit(1, "cm"), legend.position = "top")
plt.f.p.year

# P(frozen) by tend
plt.f.p.tend <-
  filter(pred.freeze, tend %% 10 == 0, tend >= 170) %>%
  mutate(tend = factor(tend)) %>%
  ggplot(aes(year, freeze.prob, group = tend)) +
  geom_ribbon(aes(ymin = f.p.lwr, ymax = f.p.upr, fill = tend), alpha = 0.1) +
  geom_line(aes(color = tend), lwd = 1) +
  labs(x = "Year", y = expression(widehat(F)[freeze](t))) +
  scale_color_manual(june.lab, values = pal[c(1:3, 7)], aesthetics = c("color", "fill")) +
  theme(legend.key.width = unit(1, "cm"), legend.position = "top")
plt.f.p.tend

# stepwise hazard
plt.f.step.year <-
  filter(pred.freeze, year %% 20 == 0) %>%
  mutate(year = factor(year)) %>%
  ggplot(aes(x = tend, group = year)) +
  geom_stephazard(aes(y = hazard, col = year), lwd = 1) +
  xlim(c(155, NA)) +
  labs(x = june.lab, y = expression(widehat(lambda)[freeze](t))) +
  scale_color_manual("Year", values = pal[c(1:3, 7)],
    aesthetics = c("color", "fill"))
plt.f.step.year

# piecewise cumulative hazard
plt.f.haz.year <-
  filter(pred.freeze, year %% 20 == 0) %>%
  mutate(year = factor(year)) %>%
  ggplot(aes(x = tend, group = year)) +
  geom_hazard(aes(y = cumu_hazard, col = year), lwd = 1) +
  xlim(c(155, NA)) +
  labs(x = june.lab, y = expression(widehat(Lambda)[freeze](t))) +
  scale_color_manual("Year", values = pal[c(1:3, 7)],
    aesthetics = c("color", "fill"))
plt.f.haz.year

# ti term
gratia::draw(pam.f, select = 3, dist = 1, rug = FALSE) # quite flat

# smooth term of year
freeze %>%
  make_newdata(tend = c(185), # January 1st
               year = seq_range(year, n = 100)) %>%
  group_by(year) %>%
  add_surv_prob(pam.f) %>%
  mutate(p = 1 - surv_prob, lwr = 1 - surv_lower, upr = 1 - surv_upper) %>%
  ggplot(aes(year, p)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = pal[1], alpha = .3) +
  geom_line(col = pal[1]) +
  ylab(Probability~of~beng~frozen~on~January~1^{st})

# estimated posterior mean and 89% CI for the freeze date
est <- 
  pred.freeze %>%
  select(year, tend, freeze.prob, f.p.lwr, f.p.upr) %>%
  group_by(year) %>%
  summarize(doy = tend[which.min(abs(0.5 - freeze.prob))],
            lwr = tend[which.min(abs(0.5 - f.p.lwr))],
            upr = tend[which.min(abs(0.5 - f.p.upr))])

# add estimated day and 95% CIs
scatter.f.est <-
  scatter.f +
  geom_ribbon(aes(year, ymin = lwr, ymax = upr), est, fill = "steelblue",
    alpha = 0.3) +
  geom_step(aes(year, doy), est, color = "steelblue", lwd = 1)
scatter.f.est

tbl_grob <- freeze %>%
  filter(year == 1976) %>%
  select(-july.year) %>%
  relocate(year, .after = last_col()) %>%
  tableGrob(theme =
    ttheme_default(padding = unit(c(3, 3), "mm")))

plt <- scatter.f.est + tbl_grob + plt.f.p.year + plt.f.p.tend +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a", tag_suffix = ".")
plt

ggsave("plots/mendota-freeze-viz.pdf", height = 3.5, width = 6.86, scale = 2)
# save.plt(plt, 'mendota-freeze-viz.pdf', height = 3.5, width = 6.86, scale = 2)
