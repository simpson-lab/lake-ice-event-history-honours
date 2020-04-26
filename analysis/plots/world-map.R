library(here) # for easy directory referencing
library(tidyverse) # for ggplot2, dplyr, etc.
library(gganimate) # for animated ggplots

theme_set(theme_bw()) # change ggplot theme

# read in ice phenology data
ice <- readRDS(here('data', 'lake-ice-data.rds')) %>% filter(Year %in% 1820:2010 & !is.na(froze.bool))

# number of frames = number of years + 10 for end_pause
nframes <- length(1820:max(ice$Year)) + 10 

# read in map data
WorldData <- map_data('world') %>% filter(lat > 0) %>% fortify() %>% as_tibble()

# create plot
plt.world <- ggplot(WorldData) +
  geom_map(map = WorldData, aes(group = group, map_id = region), fill = 'grey75', col = 'black', size = 0.5) +
  
  # add points
  geom_point(aes(longitude, latitude, col = froze.bool, group = Year), ice) +              # add points for each lake
  scale_color_manual('Froze', values = c('#ff8000', '#3366cc'), labels = c('No', 'Yes')) + # froze: y = blue, n = orange
  coord_map('azequidistant', xlim = c(-180, 180), ylim = c(30, 90)) +                      # polar projection
  
  # set scales
  scale_y_continuous(breaks = 3:8 * 10, labels = NULL) +                                   # from 30 to 80 every 10
  scale_x_continuous(breaks = NULL) +                                                      # no longitude
  labs(title = 'Year: {frame_time}', x = element_blank(), y = element_blank()) +
  transition_time(Year) +
  
  # change style of parallels and other small things
  theme(panel.grid.major = element_line(colour = "black", linetype = 'dashed'),
        legend.position = 'top',
        axis.ticks = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        plot.title = element_text(size = 24)) +
  guides(color = guide_legend(override.aes = list(size = 5)))

# animate the ggplot
animate(plt.world, duration = 30, nframes = nframes, end_pause = 10, width = 800, height = 800)

# save the animation using the here function
anim_save(here('plots', 'animated-map.gif'))
