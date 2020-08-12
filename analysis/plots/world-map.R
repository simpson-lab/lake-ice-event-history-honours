library('here')      # for easy directory referencing
library('dplyr')     # for easier data wrangling
library('ggplot2')    # for fancy plots
library('gganimate') # for animated ggplots

theme_set(theme_bw()) # change ggplot theme

# read in ice phenology data
ice <- readRDS(here('data', 'lake-ice-data.rds')) %>%
  filter(year > 1950 & !is.na(froze.bool)) %>%
  group_by(lake, long, lat) %>%
  summarise(post95 = max(year)) %>%
  mutate(post95 = post95 > 1995)

# read in map data
WorldData <- map_data('world') %>%
  filter(lat > 0) %>%
  fortify() %>%
  as_tibble()

# static plot ####
plt.world.static <-
  ggplot(WorldData) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'grey85', col = 'black', size = 0.5) +
  
  # add points
  geom_point(aes(long, lat, col = post95), ice, alpha = 0.5) +

  # polar projection
  coord_map('azequidistant') +
  
  # set scales
  scale_y_continuous(breaks = 3:8 * 10, labels = NULL) + # 30 to 80 every 10
  scale_x_continuous(breaks = NULL) +                    # no long
  labs(x = element_blank(), y = element_blank()) +
  
  # change style of parallels and other small things
  scale_color_manual('Data after 1995', values = c('#ff8000', '#3366cc')) +
  cowplot::theme_map() +
  theme(panel.grid.major = element_line(colour = "black", linetype = 'dashed'),
        legend.position = 'top', axis.ticks = element_blank())
plt.world.static
ggsave('plots/lake-map.pdf', plt.world.static, width = 10, height = 10)

# animated plot ####
# number of frames = number of years + 10 for end_pause
nframes <- length(1820:max(ice$year)) + 10 

plt.world <-
  ggplot(WorldData) +
  geom_map(map = WorldData, aes(group = group, map_id = region),
           fill = 'grey85', col = 'black', size = 0.5) +
  
  # add points
  geom_point(aes(long, lat, col = froze.bool, group = year), ice) +
  
  # froze: y = blue, n = orange
  scale_color_manual('Froze', values = c('#ff8000', '#3366cc'),
                     labels = c('No', 'Yes')) +
  
  # polar projection
  coord_map('azequidistant', xlim = c(-180, 180), ylim = c(30, 90)) +
  
  # set scales
  scale_y_continuous(breaks = 3:8 * 10, labels = NULL) + # 30 to 80 every 10
  scale_x_continuous(breaks = NULL) +                    # no long
  labs(title = 'Year: {frame_time}', x = element_blank(), y = element_blank()) +
  transition_time(year) +
  
  # change style of parallels and other small things
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
anim_save(here('plots', 'animated-map.gif'))
