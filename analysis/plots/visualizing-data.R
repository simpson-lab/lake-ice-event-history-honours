# setup ####
# data accessing
library('here')      # for easier directory referencing, possible conflicts with lubridate::here
library('readr')     # to read in files as tibbles

# data editing
library('dplyr')     # makes data editing easier
library('tibble')    # a tibble is a fancy data.frame
library('pammtools') # change data to piece-wise esponential data
library('lubridate') # deal with dates smoothly
library('tidyr')     # for expand_grid()

# graphics
library('ggplot2')   # fancy plots
library('cowplot')   # ggplot in grids
library('ggExtra')   # marginal histograms
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

# Import data ####
# read in data
ice <- read_rds(here::here('data/lake-ice-data.rds')) %>%
  filter(Year >= 1800) %>%
  mutate(continent = if_else(long > -30, 'Eurasia', 'North America'))

# map
WorldData <-
  map_data('world') %>%
  filter(lat > 0) %>%
  fortify() %>%
  as_tibble()

plt.map <- ggplot(filter(ice, !duplicated(station)), aes(long, lat)) +
  geom_map(map = WorldData, aes(group = group, map_id = region), data = WorldData,
           inherit.aes = FALSE, col = 'black', fill = 'grey90', size = 0.5) +
  coord_map('azequidistant', xlim = c(-180, 180), ylim = c(30, 90)) +
  geom_point(alpha = .2, col = pal[1], pch = 19) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = c(45, 62), labels = c('', '')) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid = element_line(linetype = 'dashed', color = 'black')); plt.map
#save.plt(plt.map, 'lake-map.pdf', height = 8, width = 8)

# lenght of time series
ice.length <-
  ice %>%
  group_by(station) %>%
  summarise(length = n(),
            Year = min(Year),
            continent = unique(continent)) %>%
  ungroup()

plt.lake.length.na <-
  ggplot(filter(ice.length, continent == 'North America'), aes(Year, length)) +
  geom_point(size = 1e-4, color = 'white') + # invisible points for marginal density plot
  geom_hex(col = 'grey') +
  labs(x = 'Start of the time series', y = 'Number of observations in the record') +
  theme(legend.position = 'none', text = element_text(size = 15)) +
  scale_fill_viridis_c('Number of lakes', option = 'B', direction = -1, limits = c(0, 20))
plt.lake.length.legend <- # plot legend
  get_legend(plt.lake.length.na +
               theme(legend.position = 'bottom', legend.key.width = unit(.4, 'in')))
plt.lake.length.na <- ggMarginal(plt.lake.length.na, type = 'histogram', fill = 'grey75')

plt.lake.length.eura <-
  ggplot(filter(ice.length, continent == 'Eurasia'), aes(Year, length)) +
  geom_point(size = 1e-4, color = 'white') + # invisible points for marginal density plot
  geom_hex(col = 'grey') +
  labs(x = 'Start of the time series', y = 'Number of observations in the record') +
  theme(legend.position = 'none', text = element_text(size = 15)) +
  scale_fill_viridis_c('Number of lakes', option = 'B', direction = -1, limits = c(0, 20))
plt.lake.length.eura <- ggMarginal(plt.lake.length.eura, type = 'histogram', fill = 'grey75')

plt.lake.length.marg <- plot_grid(NULL, plt.lake.length.na, NULL,  # NULLs add two cols
                                  NULL, plt.lake.length.eura, NULL,# of white space
                                  ncol = 3, rel_widths = c(0.3, 1, 0.3),
                                  labels = c(NA, 'a.', NA, NA, 'b.', NA))
plt.lake.length.leg <- plot_grid(plt.lake.length.marg,
                                 plt.lake.length.legend,
                                 ncol = 1, rel_heights = c(1, .05))
#save.plt(plt.lake.length.leg, 'time-series-length.pdf', width = 7, height = 10)

# sampling and location bias ####
plt.obs.dens.na <- ggplot(filter(ice, continent == 'North America'), aes(Year, lat)) +
  geom_point(size = 1e-4, color = 'white') + # invisible points for marginal density plot
  geom_hex(col = 'grey') +
  theme(legend.position = 'none', text = element_text(size = 15)) +
  ylab('Latitude (degrees North)') +
  scale_fill_viridis_c('Number of lakes', option = 'B', direction = -1, limits = c(0, 500))

plt.obs.dens.eura <- ggplot(filter(ice, continent == 'Eurasia'), aes(Year, lat)) +
  geom_point(size = 1e-4, color = 'white') + # invisible points for marginal density plot
  geom_hex(col = 'grey') +
  theme(legend.position = 'none', text = element_text(size = 15)) +
  ylab('Latitude (degrees North)') +
  scale_fill_viridis_c('Number of lakes', option = 'B', direction = -1, limits = c(0, 500))
plt.obs.dens.legend <- # plot legend
  get_legend(plt.obs.dens.eura +
               theme(legend.position = 'bottom', legend.key.width = unit(.4, 'in')))

plt.obs.dens.marg <-
  plot_grid(NULL,
            ggMarginal(plt.obs.dens.na, type = 'histogram', fill = 'grey75'),
            NULL, NULL,
            ggMarginal(plt.obs.dens.eura, type = 'histogram', fill = 'grey75'),
            NULL, ncol = 3, labels = c(NA, 'a.', NA, NA, 'b.', NA),
            rel_widths = c(.3, 1))
plt.obs.dens.marg <- plot_grid(plt.obs.dens.marg,
                               plt.obs.dens.legend,
                               ncol = 1, rel_heights = c(1, .1))
#save.plt(plt.obs.dens.marg, 'obs-density.pdf', width = 7, height = 10)

# filter record to start in 1950
ice <- filter(ice, Year >= 1950)

# Date shift ####
# plot dates
dates.fun <- function(x, ylab) {
  x +
    geom_hex(na.rm = TRUE, col = 'grey') +
    scale_fill_viridis_c('Number of events', option = 'B', direction = -1,
                         limits = c(0, 250)) +
    ylab(ylab) +
    theme(legend.position = 'none')
}

plt.dates.leg <- get_legend(dates.fun(ggplot(ice, aes(Year, On.DOY)), 'Day of year') +
                              theme(legend.position = 'bottom',
                                    legend.key.width = unit(.4, 'in')))

plt.dates <-
  plot_grid(dates.fun(ggplot(ice, aes(Year, On.DOY)), 'Day of year'),
            dates.fun(ggplot(ice, aes(Year, On.DOY.jul)), june.lab),
            dates.fun(ggplot(ice, aes(Year, Off.DOY)), 'Day of year'),
            dates.fun(ggplot(ice, aes(Year, Off.DOY.oct)), sept.lab),
            labels = c('a.', 'b.', 'c.', 'd.'), ncol = 2)
plt.dates <- plot_grid(plt.dates, plt.dates.leg, rel_heights = c(1, 0.1), ncol = 1)
#save.plt(plt.dates, 'dates.pdf', height = 8, width = 8)

# example: filter to only specific lakes
eg <- filter(ice,                           # subset to only three lakes:
             'WASCANA LAKE' == lake |       # Wascana Lake, Canada
               'MOHANSIC' == lake |         # Mohansic Lake, United States
               station == 'ARAI1',          # Lake Suwa, Japan
             Year %in% 1970:1979)           # subset to a single decade

ggplot(eg, aes(Year, On.DOY.jul)) +
  facet_grid(. ~ lake) +
  geom_point()

# change to piece-wise esponential data format
# freezing
eg.freeze <-
  ice %>%
  select(lake, station, Year, july.year, froze.bool, On.date, On.DOY,
         On.DOY.jul) %>%
  as_ped(formula = Surv(time = On.DOY.jul,           # follow-up time
                        event = froze.bool) ~ .) %>% # did the lake freeze? TRUE/FALSE
  filter(station %in% unique(eg$station) & Year %in% 1970:1979)

eg.freeze <- bind_rows(eg.freeze,
                       expand_grid(id = NA, tstart = NA, tend = max(eg.freeze$tend) + 10,
                                   interval = NA, offset = NA, ped_status = 1,
                                   lake = unique(eg$lake), Year = unique(eg$Year),
                                   july.year = NA, On.date = NA, On.DOY = NA)) %>%
  mutate(ped_status = case_when(lake == 'LAKE SUWA' & Year %in% c(1971, 1978) ~ 0,
                                lake == 'WASCANA LAKE' & Year == 1974 ~ NA_real_,
                                TRUE ~ ped_status)) %>%
  mutate(censored = lake == 'WASCANA LAKE' & Year == 1974) # only one censored observation

# thawing
eg.thaw <-
  mutate(ice, froze = froze.bool) %>%
  select(lake, station, Year, july.year, froze.bool, Off.date,
         Off.DOY, Off.DOY.oct) %>%
  as_ped(formula = Surv(time = Off.DOY.oct,             # follow-up time
                        event = froze.bool &
                          !is.na(Off.DOY.oct)) ~ .) %>% # lake thawed if it froze & date not NA
  filter(station %in% unique(eg$station) & Year %in% 1970:1979)

eg.thaw <- bind_rows(eg.thaw,
                     expand_grid(id = NA, tstart = NA, tend = max(eg.thaw$tend) + 10,
                                 nterval = NA, offset = NA, ped_status = 1,
                                 lake = unique(eg$lake), Year = unique(eg$Year), july.year = NA,
                                 Off.date = NA, Off.DOY = NA))

plt.eg.freeze <- ggplot(eg.freeze, aes(tend, ped_status)) +
  facet_grid(Year ~ lake) +
  geom_line(color = pal[1]) +
  geom_ribbon(aes(ymin = 0, ymax = ped_status), fill = pal[1], alpha = 0.5) +
  geom_hline(aes(col = censored, yintercept = 0), lwd = 1, show.legend = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = 0:1, labels = c('No', 'Yes')) +
  scale_color_manual(values = c(pal[1], pal[4])) +
  labs(x = june.lab, y = 'Frozen')

plt.eg.thaw <- ggplot(eg.thaw, aes(tend, ped_status)) +
  facet_grid(Year ~ lake) +
  geom_line(color = pal[2]) +
  geom_ribbon(aes(ymin = 0, ymax = ped_status), fill = pal[2], alpha = .5) +
  geom_hline(yintercept = 0, lwd = 1, col = pal[2]) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = 0:1, labels = c('No', 'Yes')) +
  labs(x = sept.lab, y = 'Thawed')

plt.eg <- plot_grid(plt.eg.freeze, plt.eg.thaw, nrow = 1, labels = c('a.', 'b.'), hjust = 0)
#save.plt(plt.eg, 'three-lake-example.pdf', width = 8, height = 6)
