# setup ####
# data accessing
library('here')      # for easier directory referencing

# data editing
library('dplyr')     # makes data editing easier
library('tibble')    # a tibble is a fancy data.frame

# graphics
library('ggplot2')   # fancy plots
library('ggdag')     # ggplot DAGs

# Directed acyclical graphs ####
pam.dag <-
  dagify(ice ~ wind + air.temp + stratification + depth + salinity +Z, # frozen?
         
         wind ~ air.temp,                                              # effects
         air.temp ~ radiation + continentality + altitude + Year,
         stratification ~ air.temp + depth + radiation + wind + salinity,
         continentality ~ longitude + latitude,
         radiation ~ latitude + altitude + DOY,
         Z ~ Year,
         latent = c('stratification', 'depth', 'radiation', 'continentality',
                    'salinity', 'Z'),
         exposure = c('wind', 'air.temp', 'stratification', 'depth',
                      'altitude', 'radiation', 'continentality', 'longitude',
                      'latitude', 'DOY', 'Year'),
         outcome = c('ice'),
         coords = rbind(c(name = 'DOY',            x = 6, y = 5),
                        c(name = 'Year',           x = 2, y = 0),
                        c(name = 'radiation',      x = 5, y = 8),
                        c(name = 'longitude',      x = 3, y = 0),
                        c(name = 'latitude',       x = 6, y = 0),
                        c(name = 'continentality', x = 4, y = 1.5),
                        c(name = 'air.temp',       x = 3, y = 3),
                        c(name = 'wind',           x = 0, y = 3),
                        c(name = 'altitude',       x = 5, y = 4),
                        c(name = 'depth',          x = 0, y = 8),
                        c(name = 'stratification', x = 3, y = 8),
                        c(name = 'salinity',       x = 1, y = 10),
                        c(name = 'Z',              x = 1, y = 0),
                        c(name = 'ice',            x = 0, y = 6)) %>%
           as_tibble()) %>%
  tidy_dagitty() %>%
  node_collider() %>%
  mutate(observed = if_else(name %in% c('DOY', 'Year', 'longitude', 'latitude',
                                        'air.temp', 'wind', 'altitude','ice'),
                            'Yes', 'No'),
         observed = factor(observed, levels = c('Yes', 'No')))

# DAG for single-station model
filter(pam.dag,
       ! name %in% c('longitude', 'latitude', 'depth', 'altitude',
                     'continentality') &
         ! to %in% c('longitude', 'latitude', 'depth') ) %>%
  ggdag(text_col = 'grey25', text_size = 2.5,node = FALSE, text = FALSE) +
  geom_dag_node(aes(shape = observed, col = colliders), size = 20) +
  geom_dag_text(color = 'grey25', size = 2.5) +
  geom_dag_edges_link() +
  scale_color_brewer('Collider', type = 'qual', palette = 4, direction = -1,
                     labels = c('No', 'Yes')) +
  scale_shape_manual('Observed', values = c(16, 17)) +
  theme_dag_blank() +
  theme(legend.position = 'bottom')

## DAG for hierarchical models
ggdag(pam.dag, text_col = 'grey25', text_size = 2.5,
               text = FALSE, node = FALSE) +
  geom_dag_node(aes(shape = observed, color = colliders), size = 20) +
  geom_dag_text(color = 'grey25', size = 2.5) +
  geom_dag_edges_link() +
  scale_color_brewer('Collider', type = 'qual', palette = 4, direction = -1,
                     labels = c('No', 'Yes')) +
  scale_shape_manual('Observed', values = c(16, 17)) +
  theme_dag_blank() +
  theme(legend.position = 'bottom')

# cannot condition on air temperature because it is a collider and thus
# conditioning would break the causal paths and cause estimation issues.
# However, conditioning on long, lat, altitude, year, and DOY would effectively
# account for all causes of air temperature. If we do add air.temperature,
# then the effects of year, tend, long, lat and altitude are only their effect
# after accounting for the effect of air temperature, which does not leave much.
# Thus, we could fit a model with air.temp and without the other variables, but
# that model (air.temp causes changes in hazard) is not interesting or useful,
# as it cannot measure the change over time and space (edited