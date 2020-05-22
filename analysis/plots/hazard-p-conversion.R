library('ggplot2') # for plots
library('cowplot') # for plot grids
library('dplyr')   # for data processing
library('tidyr')   # for data processing
source('functions/save.plt.R') # for saving plots easily

theme_set(theme_bw())

# create table of values
d <- tibble(log.hazard = seq(-3, 1.5, length.out = 5e4),
            hazard = exp(log.hazard),
            cumu_hazard = .5 * hazard^2 - 0, # integral (T) dT from 0 to t
            surv = exp(0 - cumu_hazard),
            p = 1 - surv)

# values when p = 0.5
filter(d, round(p, 4) == .5)

# values when p = 0.99
filter(d, round(p, 5) == .99)

# plotting function
plt.fun <- function(x) {
  d$x <- d[[x]]
  
  xlab <- case_when(x == 'log.hazard' ~ expression(log(lambda(t))),
                    x == 'hazard' ~ expression(lambda(t)),
                    x == 'cumu_hazard' ~ expression(Lambda(t)),
                    x == 'surv' ~ expression(S(t)))
  
  ggplot(d, aes(x, p)) +
    geom_hline(yintercept = 0.5, lty = 'dashed', color = 'grey') +
    geom_line() +
    labs(x = xlab, y = expression(P(T<=t)))
}

# plot
plot_grid(plt.fun('log.hazard'),
          plt.fun('hazard'),
          plt.fun('cumu_hazard'),
          plt.fun('surv'),
          labels = c('a.', 'b.', 'c.', 'd.'))
#save.plt(last_plot(), 'hazard-viz.pdf')
