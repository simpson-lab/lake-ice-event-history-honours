WorldData <-
  map_data('world') %>%
  filter(lat > 30) %>%
  fortify() %>%
  as_tibble()

plt.spatial <- function(lake.na, lake.eura, event = c('freeze', 'thaw'),
                        continent = c('na', 'ea')) {
  
  if(length(event) > 1 | !(event %in% c('freeze', 'thaw'))) {
    stop('Please choose an event, either "freeze" or "thaw"')
  }
  
  na <-
    make_newdata(get(paste0(event, '.na')),
                 lake = lake.na,
                 long = seq_range(long, 40),
                 lat = seq_range(lat, 40)) %>%
    add_surv_prob(get(paste0('pam.', event, '.na'))) %>%
    mutate(p = 1 - surv_prob) %>%
    select(p, long, lat)
  
  eura <- make_newdata(get(paste0(event, '.eura')),
                       lake = lake.eura,
                       long = seq_range(long, 40),
                       lat = seq_range(lat, 40)) %>%
    add_surv_prob(get(paste0('pam.', event, '.eura'))) %>%
    mutate(p = 1 - surv_prob) %>%
    select(p, long, lat)
  
  ggplot() +
    geom_tile(aes(long, lat, fill = p), na) +
    geom_tile(aes(long, lat, fill = p), eura) +
    geom_polygon(aes(long, lat, group = group),
                 WorldData, fill = 'transparent', color = 'grey') +
    coord_map('azequidistant', xlim = c(-180, 180), ylim = c(30, 90)) +
    labs(x = NULL, y = NULL) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    scale_fill_viridis_c(option = 'B', direction = -1)
}

x <- plt.spatial(lake.na, lake.eura, 'freeze')
