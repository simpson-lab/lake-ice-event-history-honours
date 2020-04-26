# setup ####
# data accessing
library('here')      # for easier directory referencing
library('readr')     # to read in files as tibbles

# data editing
library('dplyr')     # makes data editing easier
library('tibble')    # a tibble is a fancy data.frame

# model fitting
library('pammtools') # tools for Piecewise-exponential Additive Mixed Models
library('survival')  # survival analysis functions
library('mgcv')      # to fit Tweedie GAM
library('brms')      # to fit censored Gamma GAM

# graphics
library('ggplot2')   # fancy plots
library('cowplot')   # ggplot in grids
library('gratia')    # pretty GAM plots
library('mgcViz')    # pretty GAM diasgnostic plots
source(here::here('functions/save.plt.R')) # to save plots easily

# change ggplot theme
theme_set(theme_bw())

# palette for plots
pal <- c('#4477AA', '#ff8c00', '#66CCEE', '#009900',
         '#CCBB44', '#EE6677', '#AA3377', '#BBBBBB')
plot(1:length(pal), col = pal, cex = 5, pch = 15)

# plot labels
june.lab <- expression(Days~after~June~30^{th})
sept.lab <- expression(Days~after~September~30^{th})

stechlinsee <- read_rds(here::here('data/lake-ice-data.rds')) %>%
  filter(lake == 'STECHLINSEE')

ggplot(stechlinsee, aes(Year, duration, col = froze.bool)) +
  geom_point() +
  scale_color_manual('Froze', values = 2:1)

# freeze and thaw models ####
stech.freeze <-
  select(stechlinsee, Year, july.year, froze.bool, On.date, On.DOY, On.DOY.jul) %>%
  as_ped(formula = Surv(time = On.DOY.jul,       # follow-up time
                        event = froze.bool) ~ .) # did the lake freeze? TRUE/FALSE
stech.thaw <-
  select(stechlinsee, Year, july.year, froze.bool, Off.date, Off.DOY, Off.DOY.oct) %>%
  as_ped(formula = Surv(time = Off.DOY.oct,         # follow-up time
                        event = froze.bool &
                          !is.na(Off.DOY.oct)) ~ .) # thawed if froze & date not NA

pam.freeze.stech <-
  gam(ped_status ~                              # frozen? y/n
        s(tend, bs = 'cr', k = 10) +            # smooth of DOY ("follow-up dates")
        s(Year, bs = 'cr', k = 10) +            # smooth of year
        ti(tend, Year, bs = 'cr', k = c(5, 5)), # change in DOY over the years
      data = stech.freeze,
      family = poisson('log'),
      offset = offset,
      method = 'REML') # "fast REML"

pam.thaw.stech <-
  gam(ped_status ~                              # frozen? y/n
        s(tend, bs = 'cr', k = 10) +            # smooth of DOY ("follow-up dates")
        s(Year, bs = 'cr', k = 10) +            # smooth of year
        ti(tend, Year, bs = 'cr', k = c(5, 5)), # change in DOY over the years
      data = stech.thaw,
      family = poisson('log'),
      offset = offset,
      method = 'REML') # "fast REML"

# duration models ####
# a regular Tweedie gives a really bad fit:
# (zeros are not allowed for Gamma family)
stech.tw <- gam(duration ~ s(Year, bs = 'cr', k = 10),
                data = stechlinsee,
                family = tw(link = 'log'),
                method = 'REML')

sum(stechlinsee$duration == 0) / nrow(stechlinsee) # half of response values are 0

# qqplot
stech.tw.qq <-
  qq_plot(stech.tw, method = 'simulate', n_simulate = 1e4, level = 0.95) +
  labs(title = NULL, subtitle = NULL)

# normal density for reference
stech.tw.resid <- tibble(x = seq(-10, 10, length.out = 400),
                         d = dnorm(x, 0, sd(resid(stech.tw))))

stech.tw.hist <-
  ggplot(mapping = aes(resid(stech.tw))) +
  geom_density(alpha = .2, fill = 'black') +
  geom_rug() +
  geom_line(aes(x = x, y = d), stech.tw.resid, col = 'red',
            inherit.aes = FALSE) +
  labs(x = 'Residuals', y = 'Density')

plt.stech.tw.resid <- plot_grid(stech.tw.qq, stech.tw.hist,
                                ncol = 1, labels = c('a.', 'b.')); plt.stech.tw.resid
#save.plt(plt.stech.tw.resid, 'stechlinsee-tw-diagnostics.pdf', width = 4)

# effect of year given that the lake froze
stech.hurdle <- readRDS('analysis/models/stech-hurdle-gamma.rds')
stech.hurdle <- brm(bf(duration ~ s(Year, bs = 'cr', k = 10),
                       hu ~ 1),
                    family = hurdle_gamma(),
                    data = stechlinsee,
                    chains = 4,
                    iter = 4000,
                    cores = 4,
                    control = list(adapt_delta = 0.99999, max_treedepth = 20),
                    seed = 1)
#saveRDS(stech.hurdle, 'analysis/models/stech-hurdle-gamma.rds')

summary(stech.hurdle)

# model predictions ####
# pamms ----
## Freezing ----
# predictions grouped by year
pred.freeze.year <- 
  stech.freeze %>%
  make_newdata(tend = unique(tend),
               Year = round(seq_range(Year, n = 5))) %>%
  group_by(Year) %>%
  add_hazard(pam.freeze.stech) %>%
  add_cumu_hazard(pam.freeze.stech) %>%
  add_surv_prob(pam.freeze.stech) %>%
  mutate(freeze.prob = 1 - surv_prob,
         f.p.lwr = 1 - surv_lower,
         f.p.upr = 1 - surv_upper)

# P(frozen) by year
plt.f.p.year <- ggplot(pred.freeze.year, aes(tend, freeze.prob, group = Year)) +
  geom_line(aes(color = factor(Year)), lwd = 1) +
  geom_ribbon(aes(ymin = f.p.lwr, ymax = f.p.upr, fill = factor(Year)), alpha = 0.1) +
  labs(x = june.lab, y = 'Probability being frozen') +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal); plt.f.p.year

# stepwise hazard
plt.f.step.year <- ggplot(pred.freeze.year, aes(x = tend, group = Year)) +
  geom_stephazard(aes(y = hazard, col = factor(Year)), lwd = 1) +
  xlim(c(180, 260)) +
  labs(x = june.lab, y = expression(hat(lambda)(t))) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal); plt.f.step.year

# piecewise cumulative hazard
plt.f.haz.year <- ggplot(pred.freeze.year, aes(x = tend, group = Year)) +
  geom_hazard(aes(y = cumu_hazard, col = factor(Year)), lwd = 1) +
  xlim(c(180, 260)) +
  labs(x = june.lab, y = expression(hat(Lambda)(t))) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal); plt.f.haz.year

## Thawing ----
# predictions grouped by year
pred.thaw.year <- 
  stech.thaw %>%
  make_newdata(tend = unique(tend), 
               Year = round(seq_range(Year, n = 5))) %>%
  group_by(Year) %>%
  add_hazard(pam.thaw.stech) %>%
  add_cumu_hazard(pam.thaw.stech) %>%
  add_surv_prob(pam.thaw.stech) %>%
  mutate(thaw.prob = 1 - surv_prob,
         t.p.lwr = 1 - surv_lower,
         t.p.upr = 1 - surv_upper)

# P(ice-free) by year
plt.t.p.year <- ggplot(pred.thaw.year, aes(tend, thaw.prob, group = Year)) +
  geom_line(aes(color = factor(Year)), lwd = 1) +
  geom_ribbon(aes(ymin = t.p.lwr, ymax = t.p.upr, fill = factor(Year)), alpha = 0.1) +
  labs(x = sept.lab, y = 'Probability of being ice-free') +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal); plt.t.p.year

# stepwise hazard
plt.t.step.year <- ggplot(pred.thaw.year, aes(x = tend, group = Year)) +
  geom_stephazard(aes(y = hazard, col = factor(Year)), lwd = 1) +
  xlim(c(160, NA)) +
  labs(x = june.lab, y = expression(hat(lambda)(t))) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal); plt.t.step.year

# piecewise cumulative hazard
plt.t.haz.year <- ggplot(pred.thaw.year, aes(x = tend, group = Year)) +
  geom_hazard(aes(y = cumu_hazard, col = factor(Year)), lwd = 1) +
  xlim(c(160, NA)) +
  labs(x = june.lab, y = expression(hat(Lambda)(t))) +
  scale_color_manual('Year', values = pal) +
  scale_fill_manual('Year', values = pal); plt.t.haz.year

# duration ----
# tweedie
pred.tw <- tibble(Year = seq_min_max(stechlinsee$Year, 400))
pred.tw <- 
  cbind(pred.tw,
        predict(stech.tw, pred.tw, se.fit = TRUE)) %>%
  mutate(mu = exp(fit), upr = exp(fit + se.fit * 1.96), lwr = exp(fit - se.fit * 1.96))

# hurdle gamma
pred.hg <- tibble(Year = seq_min_max(stechlinsee$Year, 400))
pred.hg <- cbind(pred.hg, fitted(stech.hurdle, pred.hg)) %>% as_tibble()

# plot
pred <- rbind(select(pred.tw, Year, mu, lwr, upr) %>% mutate(Model = 'Tweedie'),
              transmute(pred.hg, Year = Year, mu = Estimate, lwr = Q2.5, upr = Q97.5) %>%
                mutate(Model = 'Hurdle gamma'))

plt.stech.fit <- 
  ggplot() +
  geom_point(aes(Year, duration), stechlinsee, alpha = 0.5) +
  geom_ribbon(aes(Year, ymin = lwr, ymax = upr, fill = Model), pred, alpha = 0.2) +
  geom_line(aes(Year, mu, color = Model), pred, lwd = 1) +
  ylab('Duration of ice cover') +
  scale_color_manual('Model', values = pal) +
  scale_fill_manual('Model', values = pal) +
  theme(legend.position = 'bottom', text = element_text(size = 15)); plt.stech.fit
#save.plt(plt.stech.fit, 'stechlinsee-models.pdf')
