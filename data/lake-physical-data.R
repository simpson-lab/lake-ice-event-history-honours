# attach packages
library(readr)     # read files as tibbles
library(dplyr)     # easier data editing
library(here)      # easier directory referencing

# source external function
source(here::here('functions/post.ref.date.R'))

# data is available at https://nsidc.org/data/G01377?qt-data_set_tabs=1#qt-data_set_tabs
physical <- read_csv(here('data', 'liag_physical_character_table.csv'))
physical[physical == -999] <- NA # NAs are indicated as -999
summary(physical)
glimpse(physical)
