library('here')      # for easy directory referencing
library('tidyverse') # for ggplot2, dplyr, etc.
library('gganimate') # for animated ggplots
library('ggridges')  # for ridgeline plots

theme_set(theme_bw()) # change ggplot theme

# read in data
ice <- readRDS(here('data', 'lake-ice-data.rds')) %>%
  mutate(Year = year) %>%
  filter(Year >= 1820 & Year %% 5 == 0)

# number of frames = number of years + 10 for end_pause
nframes <- length(1820:max(ice$Year)) + 10 

## create plots
# density of freeze events
plt.dens.freeze <- ggplot(ice, aes(x = On.DOY.jul)) +
  geom_density_line(fill = '#67a9cf') +
  geom_point(aes(y = 0), alpha = 0.25) +
  labs(title = 'Year: {frame_time}', x = 'Days of year after June 20',
       y = 'Density of freeze events') +
  scale_y_continuous(breaks = NULL) +
  transition_time(Year) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        plot.title = element_text(size = 24))

# density of thaw events
plt.dens.thaw <- ggplot(ice, aes(x = Off.DOY.oct)) +
  geom_density_line(fill = '#ef8a62') +
  geom_point(aes(y = 0), alpha = 0.25) +
  labs(title = 'Year: {frame_time}', x = 'Days after January 1',
       y = 'Density of thaw events') +
  scale_y_continuous(breaks = NULL) +
  transition_time(Year) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        plot.title = element_text(size = 24))

## animate plots
anim.dens.freeze <- animate(plt.dens.freeze, duration = 30, nframes = nframes,
                            end_pause = 10, width = 800, height = 350)
anim.dens.thaw <- animate(plt.dens.thaw, duration = 30, nframes = nframes,
                          end_pause = 10, width = 800, height = 350)

# save plots
anim_save(here('plots', 'animated-density-freeze.gif'),
          animation = anim.dens.freeze)
anim_save(here('plots', 'animated-density-thaw.gif'),
          animation = anim.dens.thaw)

#### Note: #################
# this plot works
ggplot(subset(ice, Year %in% 1820:2010), aes(x = Off.DOY)) +
  geom_density_line(fill = '#ef8a62') +
  geom_point(aes(y = 0), alpha = 0.25) +
  labs(title = 'Year: {frame_time}', x = 'Days after January 1', y = NULL) +
  transition_time(Year) +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        plot.title = element_text(size = 24))

# this plot crops density plots at the end of the animation (Year >= 1990)
ggplot(subset(ice, Year %in% 1800:2010), aes(x = Off.DOY)) +
  geom_density_line(fill = '#ef8a62') +
  geom_point(aes(y = 0), alpha = 0.25) +
  labs(title = 'Year: {frame_time}', x = 'Days after January 1', y = NULL) +
  transition_time(Year) +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        plot.title = element_text(size = 24))
