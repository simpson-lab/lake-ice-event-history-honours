# setup ####
# data accessing
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
library('gganimate') # for animated ggplots
source('functions/save.plt.R') # to save plots easily

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
  mutate(continent = if_else(long > -30, 'Eurasia', 'North America'))

# static map ####
WorldData <-
  map_data('world') %>%
  filter(lat > 0) %>%
  fortify() %>%
  as_tibble()

lakes <-
  group_by(ice, lake) %>%
  mutate(last = max(year),
         post1995 = case_when(last > 1995 ~ 'Yes',
                              last > 1950 ~ 'No',
                              TRUE ~ 'No data after 1950 (removed)') %>%
           factor(levels = c('Yes', 'No', 'No data after 1950 (removed)'))) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(post1995)

plt.map <-
  ggplot(lakes, aes(long, lat)) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           data = WorldData, inherit.aes = FALSE, col = 'grey30',
           fill = 'grey90', size = 0.5) +
  coord_map('azequidistant', xlim = c(-180, 180), ylim = c(30, 90)) +
  geom_point(aes(color = post1995, shape = post1995, alpha = post1995,
                 size = post1995)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = c(45, 62), labels = c('', '')) +
  scale_color_manual('Data after 1995', values = pal[c(1, 3, 2)]) +
  scale_shape_manual('Data after 1995', values = c(19, 19, 17)) +
  scale_alpha_manual('Data after 1995', values = c(0.75, 0.75, 1)) +
  scale_size_manual('Data after 1995', values = c(0.75, 0.75, 1.5)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = 'top') +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1)))
plt.map
#save.plt(plt.map, 'lake-map.pdf', height = 4, width = 4, scale = 2)

# only lakes with no data after 1950 
group_by(ice, lake) %>%
  mutate(last = max(year)) %>%
  filter(last < 1950, !duplicated(lake)) %>%
  select(continent, lake, last, long, lat) %>%
  ggplot(aes(long, lat)) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           data = WorldData, inherit.aes = FALSE, col = 'black',
           fill = 'grey90', size = 0.5) +
  coord_map('azequidistant', xlim = c(-180, 180), ylim = c(30, 90)) +
  geom_point(color = 'red', size = 3, alpha = 0.5) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank())

# animated map ####
# number of frames = number of years + 10 for end_pause
nframes <- length(1950:max(ice$year)) + 10

plt.world <-
  ggplot(WorldData) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'grey85', col = 'black', size = 0.5) +
  
  # add points
  geom_point(aes(long, lat, col = froze.bool, group = year),
             filter(ice, year >= 1950)) +
  
  # froze: y = blue, n = orange
  scale_color_manual('Froze', values = c('#ff8000', '#3366cc'),
                     labels = c('No', 'Yes')) +
  
  # polar projection
  coord_map('azequidistant', xlim = c(-180, 180), ylim = c(30, 90)) +
  
  # set scales
  scale_y_continuous(breaks = 3:8 * 10, labels = NULL) + # 30 to 80 every 10
  scale_x_continuous(breaks = NULL) +                    # no long
  labs(title = 'Year: {frame_time}', x = element_blank(), y = element_blank()) +
  transition_time(year) + # adds the animation layer
  
  # change style of parallels and other small things
  theme_minimal() +
  theme(panel.grid.major = element_line(colour = "black", linetype = 'dashed'),
        legend.position = 'top',
        axis.ticks = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        plot.title = element_text(size = 24)) +
  guides(color = guide_legend(override.aes = list(size = 5)))

# animate the ggplot
animate(plt.world, duration = 30, nframes = nframes, end_pause = 10,
        width = 800, height = 800)

# save the animation using the here function
anim_save('plots/animated-map-1950.gif')

# lenght of time series ####
ice.length <-
  ice %>%
  group_by(station) %>%
  summarise(length = n(),
            year = min(year),
            continent = unique(continent)) %>%
  ungroup()

plt.lake.length.na <-
  ggplot(filter(ice.length, continent == 'North America'), aes(year, length)) +
  geom_point(color = 'transparent') + # invisible for marginal density plot
  geom_hex(col = 'grey') +
  labs(x = 'Start of the time series',
       y = 'Number of observations in the record') +
  theme(legend.position = 'none', text = element_text(size = 15)) +
  scale_fill_viridis_c('Number of lakes', option = 'B', direction = -1,
                       limits = c(0, 20))
plt.lake.length.legend <- # plot legend
  get_legend(plt.lake.length.na +
               theme(legend.position = 'bottom',
                     legend.key.width = unit(.4, 'in')))
plt.lake.length.na <- ggMarginal(plt.lake.length.na, type = 'histogram',
                                 fill = 'grey75')
plt.lake.length.eura <-
  ggplot(filter(ice.length, continent == 'Eurasia', year >=1800),
         aes(year, length)) +
  geom_point(color = 'transparent') + # invisible for marginal density plot
  geom_hex(col = 'grey') +
  labs(x = 'Start of the time series',
       y = 'Number of observations in the record') +
  theme(legend.position = 'none', text = element_text(size = 15)) +
  scale_fill_viridis_c('Number of lakes', option = 'B', direction = -1,
                       limits = c(0, 20))
plt.lake.length.eura <- ggMarginal(plt.lake.length.eura, type = 'histogram',
                                   fill = 'grey75')
plt.lake.length.marg <-
  plot_grid(NULL, plt.lake.length.na, NULL,  # NULLs add two cols
            NULL, plt.lake.length.eura, NULL,# of white space
            ncol = 3, rel_widths = c(0.1, 1, 0.1),
            labels = c(NA, 'a.', NA, NA, 'b.', NA), label_x = - 0.05)
# plt.lake.length.marg <-
#   plot_grid(plt.lake.length.na, plt.lake.length.eura,
#             ncol = 1, labels = c('a.', 'b.'))
plt.lake.length.leg <- plot_grid(plt.lake.length.marg,
                                 plt.lake.length.legend,
                                 ncol = 1, rel_heights = c(1, 0.1))
#save.plt(plt.lake.length.leg, 'time-series-length.pdf', width = 8, height = 10)

# sampling and location bias ####
plt.obs.dens.na <-
  ggplot(filter(ice, continent == 'North America'), aes(year, lat)) +
  geom_point(color = 'transparent') + # invisible for marginal density plot
  geom_hex(col = 'grey') +
  theme(legend.position = 'none', text = element_text(size = 15)) +
  labs(x = 'Year', y = 'Latitude (degrees North)') +
  ylim(c(35, 85)) +
  scale_fill_viridis_c('Number of lakes', option = 'B', direction = -1,
                       limits = c(0, 500))

plt.obs.dens.eura <- ggplot(filter(ice, continent == 'Eurasia', year >= 1800),
                            aes(year, lat)) +
  geom_point(color = 'transparent') + # invisible for marginal density plot
  geom_hex(col = 'grey') +
  theme(legend.position = 'none', text = element_text(size = 15)) +
  labs(x = 'Year', y = 'Latitude (degrees North)') +
  ylim(c(35, 85)) +
  scale_fill_viridis_c('Number of lakes', option = 'B', direction = -1,
                       limits = c(0, 500))
plt.obs.dens.legend <- # plot legend
  get_legend(plt.obs.dens.eura +
               theme(legend.position = 'bottom',
                     legend.key.width = unit(.4, 'in')))

plt.obs.dens.marg <-
  plot_grid(NULL,
            ggMarginal(plt.obs.dens.na, type = 'histogram', fill = 'grey75'),
            NULL, NULL,
            ggMarginal(plt.obs.dens.eura, type = 'histogram', fill = 'grey75'),
            NULL,
            ncol = 3, labels = c(NA, 'c.', NA, NA, 'd.', NA),
            rel_widths = c(0.1, 1, 0.1), label_x = - 0.05)
plt.obs.dens.marg <- plot_grid(plt.obs.dens.marg, plt.obs.dens.legend,
                               ncol = 1, rel_heights = c(1, 0.1))
#save.plt(plt.obs.dens.marg, 'obs-density.pdf', width = 8, height = 10)

plt.hex <-
  plot_grid(
    NULL, # white space
    # time series length
    plot_grid(plt.lake.length.na, plt.lake.length.eura,
              labels = c('a.', 'b.')) %>%
      plot_grid(plt.lake.length.legend, ncol = 1, rel_heights = c(1, 0.1)),
    NULL, # white space
    # observation density
    plot_grid(ggMarginal(plt.obs.dens.na, type = 'histogram', fill = 'grey75'),
              ggMarginal(plt.obs.dens.eura, type = 'histogram', fill = 'grey75'),
              labels = c('c.', 'd.')) %>%
      plot_grid(plt.obs.dens.legend, ncol = 1, rel_heights = c(1, 0.1)),
    NULL, # white space
    ncol = 1, rel_heights = c(0.02, 1, 0.05, 1, 0.02))
plt.hex
#save.plt(plt.hex, 'ts-hex.pdf', width = 5, height = 5, scale = 2)

# filter record to start in 1950
ice <- filter(ice, year >= 1950)

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

plt.dates.leg <- get_legend(dates.fun(ggplot(ice, aes(year, On.DOY)),
                                      'Day of year') +
                              theme(legend.position = 'bottom',
                                    legend.key.width = unit(.4, 'in')))

plt.dates <-
  plot_grid(dates.fun(ggplot(ice, aes(year, On.DOY)), 'Day of year'),
            dates.fun(ggplot(ice, aes(year, On.DOY.jul)), june.lab),
            dates.fun(ggplot(ice, aes(year, Off.DOY)), 'Day of year'),
            dates.fun(ggplot(ice, aes(year, Off.DOY.oct)), sept.lab),
            labels = c('a.', 'b.', 'c.', 'd.'), ncol = 2)
plt.dates <- plot_grid(plt.dates, plt.dates.leg, rel_heights = c(1, 0.1),
                       ncol = 1)
#save.plt(plt.dates, 'dates.pdf', height = 8, width = 8)

# example: filter to only specific lakes
eg <- filter(ice,                           # subset to only three lakes:
             'WASCANA LAKE' == lake |       # Wascana Lake, Canada
               'MOHANSIC' == lake |         # Mohansic Lake, United States
               station == 'ARAI1',          # Lake Suwa, Japan
             year %in% 1970:1974)           # subset to a single decade

ggplot(eg, aes(year, On.DOY.jul)) +
  facet_grid(. ~ lake) +
  geom_point()

# change to piece-wise esponential data format
# freezing
eg.freeze <-
  ice %>%
  select(lake, station, year, july.year, froze.bool, On.date, On.DOY,
         On.DOY.jul) %>%
  as_ped(formula = Surv(time = On.DOY.jul,           # follow-up time
                        event = froze.bool) ~ .) %>% # did lake freeze? T/F
  filter(station %in% unique(eg$station) & year %in% 1970:1974)

eg.freeze <-
  bind_rows(eg.freeze,
            expand_grid(id = NA, tstart = NA, tend = max(eg.freeze$tend) + 10,
                        interval = NA, offset = NA, ped_status = 1,
                        lake = unique(eg$lake), year = unique(eg$year),
                        july.year = NA, On.date = NA, On.DOY = NA)) %>%
  mutate(ped_status = case_when(lake == 'LAKE SUWA' & year %in% c(1971, 1978)~0,
                                lake == 'WASCANA LAKE' & year == 1974~ NA_real_,
                                TRUE ~ ped_status)) %>%
  mutate(censored = lake == 'WASCANA LAKE' & year == 1974) # only 1 censored obs

# thawing
eg.thaw <-
  mutate(ice, froze = froze.bool) %>%
  select(lake, station, year, july.year, froze.bool, Off.date,
         Off.DOY, Off.DOY.oct) %>%
  as_ped(formula = Surv(time = Off.DOY.oct,             # follow-up time
                        event = froze.bool & # lake thawed if froze & date ! NA
                          !is.na(Off.DOY.oct)) ~ .) %>% 
  filter(station %in% unique(eg$station) & year %in% 1970:1974)

eg.thaw <- bind_rows(eg.thaw,
                     expand_grid(id = NA, tstart = NA,
                                 tend = max(eg.thaw$tend) + 10,
                                 nterval = NA, offset = NA, ped_status = 1,
                                 lake = unique(eg$lake), year = unique(eg$year),
                                 july.year = NA, Off.date = NA, Off.DOY = NA))

plt.eg.freeze <- ggplot(eg.freeze, aes(tend, ped_status)) +
  facet_grid(year ~ lake) +
  geom_line(color = pal[1]) +
  geom_ribbon(aes(ymin = 0, ymax = ped_status), fill = pal[1], alpha = 0.5) +
  geom_hline(aes(col = censored, yintercept = 0), lwd = 1, show.legend = FALSE)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = 0:1, labels = c('No', 'Yes')) +
  scale_color_manual(values = c(pal[1], pal[4])) +
  labs(x = june.lab, y = 'Frozen')

plt.eg.thaw <- ggplot(eg.thaw, aes(tend, ped_status)) +
  facet_grid(year ~ lake) +
  geom_line(color = pal[2]) +
  geom_ribbon(aes(ymin = 0, ymax = ped_status), fill = pal[2], alpha = .5) +
  geom_hline(yintercept = 0, lwd = 1, col = pal[2]) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = 0:1, labels = c('No', 'Yes')) +
  labs(x = sept.lab, y = 'Thawed')

plt.eg <- plot_grid(plt.eg.freeze, plt.eg.thaw, nrow = 1, # warnings are ok
                    labels = c('a.', 'b.'), hjust = 0)
#save.plt(plt.eg, 'three-lake-example.pdf', width = 8, height = 4)
