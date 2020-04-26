# attach packages
library(here)      # easier directory referencing
library(readr)     # read files as tibbles
library(dplyr)     # easier data editing

# data is available at: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.267.2
NCC.data <- read_csv(here('data', 'LakeIceIncidencewithCharacteristics_Final.csv'))
glimpse(NCC.data)
