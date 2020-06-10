# setup ####
# attach packages
library('tibble')    # fancy data frames
library('readr')     # read files as tibbles
library('dplyr')     # easier data editing
library('lubridate') # smoother date converting
library('here')      # easier directory referencing; possible conflict with lubridate::here
library('sp')        # spatial features in maps
library('rgdal')     # to write kml files
library('ggplot2')   # fancy plots
library('cowplot')   # ggplots in grids
library('ggmap')     # maps using ggplot
library('geonames')  # to convert coordinates to altitude
# set username using options(geonamesUsername = '<myusername>')
# (note that online services must be activated first in the geonames account)

theme_set(theme_bw())
VIEW.COMMENTS <- FALSE # to prevent view() functions if sourcing file

# source function to convert DOY to DOY after June 30 or September 30 (see note below) ####
source(here::here('functions/post.ref.date.R'))

# import and set up data ####
# Freezing Dates (Global Lake and River Ice Phenology Database)
# data is available at https://nsidc.org/data/G01377?qt-data_set_tabs=1#qt-data_set_tabs
ice <- read_csv(here('data', 'liag_freeze_thaw_table.csv')) %>%
  filter(lakeorriver == 'L') # lakes only (no rivers)

ice[ice == -999] <- NA # Change -999 values to NA

# update locations using data from https://doi.org/10.1038/s41558-018-0393-5
# (available at: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.267.2)
NCC.data <- read_csv(here('data', 'LakeIceIncidencewithCharacteristics_Final.csv'))

# updated coordinates from NCC.data (Sharma et al. 2019)
ice <- left_join(ice,
                 NCC.data %>%
                   select(lakecode, lakename, Latitude_dd, Longitude_dd,
                          Elevation_m)) %>%
  # if coord from NCC.data are NA take values from ice
  mutate(long = case_when(is.na(Longitude_dd) ~ longitude,
                          TRUE ~ Longitude_dd),
         lat = case_when(is.na(Latitude_dd) ~ latitude,
                         TRUE ~ Latitude_dd))

# remove rows with no coordinates
ice <- filter(ice, !is.na(long) | !is.na(lat)) # drop lakes with unknown location

# rename for easier referencing
ice <- rename(ice,
              lake = lakename,
              station = lakecode) # observation station

# lake thaws before freezing: freezing year is wrong
if(VIEW.COMMENTS) view(filter(ice, iceoff_year == 1939, lake == 'SIMCOE'))
ice <- mutate(ice,
              iceon_year = if_else(iceoff_year == 1939 & ice$lake == 'SIMCOE',
                                  true = 1938, false = iceon_year))

# add new variables ####
ice <- mutate(ice,
              # Year at the beginning of the season
              Year = as.numeric(substr(season, 1, case_when(nchar(as.character(season)) == 7 ~ 4,
                                                            nchar(as.character(season)) == 5 ~ 3))),
              
              # Year starting in July, ~ like academic year
              july.year = paste0(Year, '-', Year + 1),
              
              # first day of ice cover
              On.date = date(case_when(is.na(iceon_day) | is.na(iceon_month) | is.na(iceon_year) ~ NA_character_,
                                       TRUE ~ paste(iceon_year, iceon_month, iceon_day, sep = '-'))),
              
              # last day of ice cover + 1 = first day without ice
              Off.date = date(case_when(is.na(iceoff_day) | is.na(iceoff_month) | is.na(iceoff_year) ~ NA_character_,
                                        TRUE ~ paste(iceoff_year, iceoff_month, iceoff_day, sep = '-'))),
              
              # DOY of ice on
              On.DOY = yday(On.date),
              
              # DOY of ice off
              Off.DOY = yday(Off.date))

# must seprate into two mutate rounds because of references in post.ref.date
ice <- mutate(ice,
              
              # DOY of ice on post June 30 (year starts in July)
              On.DOY.jul = post.ref.date(column = 'On.date', d = ice, event = 'freeze'),
              
              # DOY of ice off post Sept 30 (year starts in October)
              Off.DOY.oct = post.ref.date(column = 'Off.date', d = ice, event = 'thaw'),
              
              # period of ice cover
              Off.On = as.numeric(Off.date - On.date),
              
              # difference between period and number of days with ice 
              Off.On.Duration = Off.On - duration,
              
              # boolean version of froze column
              froze.bool = case_when(froze == 'Y' ~ TRUE,
                                     froze == 'N' ~ FALSE,
                                     froze == 'U' ~ NA),
              
              # lakes that did not freeze in a given year get a DOY of -1
              On.DOY.jul.plot = case_when(is.na(froze.bool) ~ On.DOY.jul,
                                          froze.bool == TRUE ~ On.DOY.jul,
                                          froze.bool == FALSE ~ -1),
              
              # lakes that did not thaw in a given year get a DOY of -1
              Off.DOY.oct.plot = case_when(is.na(froze.bool) ~ Off.DOY.oct,
                                           froze.bool == TRUE ~ Off.DOY.oct,
                                           froze.bool == FALSE ~ -1))

# day of year is 'days after December 31 (of the previous year)' so for freeze and thaw events we can use
# days after June 30 for freeze (i.e. July 1 is 1 and June 30 is 365 or 366), and
# days after September 30 for thaw (i.e. October 1 is 1 and September 30 is 365 or 366)
plt.ref.dates <- 
  plot_grid(ggplot(ice, aes(On.DOY, 1)) +
              geom_point(alpha = .1) +
              geom_vline(aes(xintercept = yday('2000-06-30')),
                         col = 'cornflowerblue', lwd = 1) +
              labs(title = 'Freeze dates', y = element_blank(), x = 'DOY') +
              scale_y_continuous(breaks = NULL),
            ggplot(ice, aes(Off.DOY, 1)) +
              geom_point(alpha = .1) +
              geom_vline(aes(xintercept = yday('2000-09-30')),
                         col = 'darkorange', lwd = 1) +
              labs(title = 'Thaw dates', y = element_blank(), x = 'DOY') +
              scale_y_continuous(breaks = NULL),
            ncol = 1)

# check locations using GoogleMaps ####
stations <-
  select(ice, station, lake, long, lat, Elevation_m) %>%
  filter(!duplicated(station)) # one row per station

# add missing elevations
# update stations and add altitude
stations <-
  mutate(stations,
         altitude = purrr::map2_dbl(lat, long,
                                    function(a, b) {
                                      # GNgtopo30 is more accurate than GNsrtm3? 
                                      x <- GNgtopo30(lat = a, lng = b)
                                      as.numeric(x$gtopo30[1])
                                    }),
         altitude = if_else(is.na(Elevation_m), altitude, Elevation_m))
stations <- select(stations, -Elevation_m)
#write_csv(stations, 'data/stations.csv')

# some locations have very inaccurate altitudes
filter(stations, altitude < 0)

# change problematic altitudes (assuming lakes are coastal)
stations <- mutate(stations,
                   altitude = if_else(altitude > 0, altitude, 1))

# add altitude to individual observations
ice <- left_join(ice, select(stations, station, altitude), by = 'station')

# Fix locations ####
# did not change location if: - the pin was near the coast,
#                             - the pin was in the lake,
#                             - the pin was 0.01 degrees near the lake,
#                             - the lake could not be found nearby.
# Changed lake names for pins with different names but for the same lake.
# Did not change lake names or coordinates for pins in a different lake.
# lake map available at https://drive.google.com/open?id=1ZxMmIzqDVm8jBMm655EHgcMCR4KE43QJ&usp=sharing
# green = unchanged; red = changed location or name; yellow = problematic; purple = misc.

# issues:
#       - the lake labelled as "Salosjarvi" actually has coordinates for Pettama. Salosjarvi is further north.
#       - "LITTLE LAKE (BRADLEY)" renamed to "BRADLEY LAKE" (near Valmy, WI, USA)
#       - "CLOVER POND", CT was left unchanged because the point is only 1 km away from the (very small) pond
#       - "PRIOR" was assumend to be Upper Prior Lake since "LOWER PRIOR" is also present
#       - "KEITELE SW" and "KEITELE SE" considered as separate lakes due to distance and branching of the lake
#       - "LITTLE HAWK" changed to "BIG HAWK" because of location
#       - "DEASE LAKE" left unchanged because could not find a lake with that name (only a nearby town)

# not found: "CRYSTAL BOG", "CROOKED LAKE" (Minnesota), "DORSET LAKES MEAN (MAIN BASIN)", "LAKE SUUREEN-KUEL",
#            "BIG KANDI", "LAKE TCHUKCHAGIRSKOYE", "LAKE ULAKHAN-KUEL", "ADIES LAKE", "LAKE EMANDGA", "LAKE TCHALPAN",
#            "LAKE BOLSHIE CHANY", "POLICE LAKE", "MISSION LAKE", "DOG POND LAKE", "MODULE LAKE", "WATER LAKE",
#            "BAGNELL LAKE", "ATTILA LAKE", "BACK BAY" (no lake found, only location and bay), "CARMI LAKE" (in BC),
#            "FENOR LAKE", "LAKE ECORCES" (not nearby), "LITTLE GILLAM LAKE", "MCKENZIE LAKE" (none nearby),
#            "PILING LAKE", "SMALL LAKE", "TARDIF LAKE", "WACOUNO LAKE", "WINTER HOUSE POND"
#            (on other side of Newfoundland?).

ice <- mutate(ice,
              # if the lake name contains a string on the left, change its name the the name after the "~"
              lake = case_when(grepl('SUWA', lake) ~ 'LAKE SUWA',               
                               grepl('LAKE SUPERIOR', lake) |
                                 lake == 'PORT ARTHUR HARBOUR' ~ 'LAKE SUPERIOR',
                               lake == 'CHEQUAMEGON BAY' ~ 'LAKE SUPERIOR',
                               lake == 'THUNDER BAY' ~ 'LAKE SUPERIOR',
                               grepl('TROUT LAKE', lake) &
                                 !grepl('BIG', lake) ~ 'TROUT LAKE',
                               grepl('KENOGAMISIS', lake) ~ 'KENOGAMISIS LAKE',
                               grepl('UPPER TWIN LAKE AND LOWER TWIN LAKE', lake) ~ 'UPPER TWIN LAKE',
                               grepl('LAKE ONTARIO', lake) | grepl('TORONTO', lake) ~ 'LAKE ONTARIO',
                               grepl('LANGELMAVESI', lake) ~ 'LANGELMAVESI',
                               grepl('OULUJARVI', lake) ~ 'OULUJARVI',
                               grepl('JOUTSJARVI', lake) ~ 'JOUTSJARVI',
                               grepl('VATIAJARVI', lake) ~ 'VATIANJARVI', # spelling mistake: vatiajarvi -> vatiaNjarvi
                               grepl('VATIANJARVI', lake) ~ 'VATIANJARVI',
                               grepl('SALOSJARVI', lake) ~ 'PETTAMA',     # wrong lake name?
                               grepl('JOUTSJARVI', lake) ~ 'JOUTSJARVI',
                               grepl('OULUJARVI', lake) ~ 'OULUJARVI',
                               grepl('KIVIJARVI', lake) ~ 'KIVIJARVI',
                               station == 'WRS395' ~ 'TROUT LAKE, ON',         # in Ontario, Canada
                               station %in% c('JJM7', 'KMS22') ~ 'TROUT LAKE', # in Wisconsin, USA
                               grepl('BRADLEY', lake) ~ 'BRADLEY LAKE',        # originally "LITTLE LAKE (BRADLEY)"
                               grepl('VESIJARVI', lake) ~ 'VESIJARVI',
                               grepl('LANGELMAVESI', lake) ~ 'LANGELMAVESI',
                               grepl('PAIJANNE', lake) ~ 'PAIJANNE',
                               grepl('CHAUTAUQUA', lake) ~ 'CHAUTAUQUA',
                               grepl('VATTERN', lake) ~ 'LAKE VATTERN',
                               grepl('INARI', lake) ~ 'INARI',
                               grepl('BAIKAL', lake) ~ 'LAKE BAIKAL',
                               grepl('BAYKAL', lake) ~ 'LAKE BAIKAL',
                               grepl('ONEGA', lake) ~ 'LAKE ONEGA',
                               grepl('QINGHAI LAKE', lake) ~ 'QINGHAI LAKE',
                               grepl('LAKE OF BAYS', lake) ~ 'LAKE OF BAYS',
                               grepl('LITTLE HAWK', lake) ~ 'BIG HAWK',
                               grepl('ERIE', lake) |
                                 lake == 'BUFFALO HARBOR TO LONG POINT' |
                                 lake == 'INNER BAY OF LONG POINT BAY (LAKE E' ~ 'LAKE EIRIE',
                               lake == 'HAMILTON HARBOUR' |
                                 lake == 'TORONTO HARBOUR' |
                                 lake == 'PICTON BAY' |
                                 lake == 'LAKE ONTARIO - HORSEY BAY' ~ 'LAKE ONTARIO',
                               lake == 'GREEN BAY AT MENOMINEE' |
                                 lake == 'GRAND TRAVERSE BAY' ~ 'LAKE MICHIGAN',
                               lake == 'COLPOYS BAY' |
                                 lake == 'BEAUSOLEIL BAY' |
                                 lake == 'GEORGIAN BAY' |
                                 lake == 'SOUTH BAY (LAKE HURON)' ~ 'LAKE HURON',
                               lake == 'NISUTLIN BAY' ~ 'TESLIN LAKE',
                               grepl('GREAT SLAVE LAKE', lake) ~ 'GREAT SLAVE LAKE',
                               grepl('UNNAMED', lake) ~ station,
                               grepl('SMALL LAKE NEAR STATION ', lake) ~ station,
                               grepl('LAKE ADJACENT TO STATION', lake) ~ station,
                               grepl('KENOGAMISIS', lake) ~ 'KENOGAMISIS LAKE',
                               grepl('KOIVUJARVI', lake) ~ 'LAKE KOIVUJARVI',
                               grepl('CACHE', lake) ~ 'LAC CACHE',
                               station == 'WRS321' ~ 'LONG LAKE (ONTARIO)',
                               
                               # otherwise keep it the same
                               TRUE ~ lake),
              
              # fix latitude
              lat = case_when(grepl('SUWA', lake) ~ 36.04871,
                                   grepl('VATIANJARVI', lake) ~ 62.48,
                                   grepl('SALOSJARVI', lake) ~ 62.25836,
                                   grepl('REHJA', lake) ~ 64.20207,
                                   grepl('NUASJARVI', lake) ~ 64.15193,
                                   lake == 'ANDERSON LAKE' ~ 46.17072,
                                   lake == 'FISHER LAKE' ~ 45.91931,
                                   lake == 'SPRUCE LAKE' ~ 46.05276,
                                   lake == 'MYSTERY LAKE' ~ 46.05646,
                                   lake == 'LAKE NIPISSING' ~ 46.24391,
                                   lake == '(BIG) LONG LAKE' ~ 45.75577,
                                   lake == 'BLACK OTTER LAKE' ~ 44.333,
                                   grepl('LAKE NASIJARVI', lake) ~ 61.69841,
                                   grepl('KITUSJARVI (3548)', lake) ~ 62.27939,
                                   lake == 'JAASJARVI - HARTOLA (1457)' ~ 61.61585,
                                   lake == 'POROVESI' ~ 63.54326,
                                   lake == 'PALOVESI' ~ 61.89657,
                                   lake == 'KUIVAJARVI' ~ 60.78028,
                                   lake == 'DEEP' ~ 43.02253,
                                   lake == 'CHAUTAUQUA SOUTH' ~ 42.12263,
                                   lake == 'DUCK LAKE' ~ 42.38561,
                                   lake == 'FOUNTAIN' ~ 43.65612,
                                   lake == 'DETROIT LAKE' ~ 46.7952,
                                   lake == 'LOWER PRIOR' ~ 44.73465,
                                   lake == 'PRIOR' ~ 44.71361,
                                   lake == 'LEECH' ~ 47.14375,
                                   lake == 'NIINIVESI W' ~ 62.78378,
                                   lake == 'KYNSIVESI W' ~ 62.40342,
                                   lake == 'LAKE KUBENSKOYE (PESKI)' ~ 59.53351,
                                   lake == 'LAKE B.LEPRINDO' ~ 56.61743,
                                   lake == 'LAKE NICHATKA' ~ 57.62684,
                                   lake == 'LAKE ASLI-KUL' ~ 54.31316,
                                   lake == 'LAKE EMANDGA' ~ 65.26671,
                                   lake == 'LAKE SYABERO' ~ 58.81276,
                                   lake == 'LAKE ASLI-KUL' ~ 54.30517,
                                   lake == 'MINNOW LAKE' ~ 46.49246,
                                   lake == 'MARSH LAKE' ~ 60.51123,
                                   lake == 'ATLIN LAKE' ~ 60.72,
                                   lake == 'ANAHIM LAKE' ~ 52.50532,
                                   lake == 'ATLIN LAKE' ~ 59.56296,
                                   grepl('KENOGAMISIS', lake) ~ 49.70625,
                                   station == 'WRS206' ~ 49.68991,
                                   station == 'RB19' ~ 63.50447,
                                   lake == 'BEAR LAKE' ~ 55.22172,
                                   lake == 'LAC CACHE' ~ 49.83469,
                                   lake == 'CANOE LAKE' ~ 55.1545,
                                   lake == 'CHARLIE LAKE' ~ 56.27729,
                                   lake == 'LITTLE KETTLE LAKE' ~ 56.128077,
                                   lake == 'GILMAN LAKE' ~ 49.9114,
                                   lake == 'GLENMORE RESERVOIR ELBOW RIVER' ~ 50.98058,
                                   lake == 'GRAND LAKE' ~ 46.16428,
                                   lake == 'LAC MALFAIT' ~ 48.57695,
                                   lake == 'LAKE ALEX' ~ 49.265102,
                                   lake == 'LAKE DORE' ~ 49.86725,
                                   lake == 'LAKE ONATCHIWAY' ~ 49.0385,
                                   lake == 'LANDING LAKE' ~ 55.28424,
                                   station == 'WRS321' ~ 46.71514,
                                   lake == 'LYNN LAKE' ~ 56.83653,
                                   lake == 'MC LEOD LAKE' ~ 54.99179,
                                   lake == 'MEADOW LAKE' ~ 54.11209,
                                   lake == 'NUKKO LAKE' ~ 54.0722,
                                   lake == 'OBA LAKE' ~ 48.63629,
                                   lake == 'PLAYGREEN LAKE' ~ 53.89415,
                                   lake == 'PORTAGE BAY' ~ 51.67219,
                                   lake == 'UPPER TWIN LAKE AND LOWER TWIN LAKE' ~ 50.15144,
                                   lake == 'WILLIAMS LAKE' ~ 52.11761,
                                   TRUE ~ lat),
              
              # fix longitude
              long = case_when(grepl('SUWA', lake) ~ 138.0837,
                                    grepl('VATIANJARVI', lake) ~ 25.89999,
                                    grepl('SALOSJARVI', lake) ~ 25.03221,
                                    grepl('REHJA', lake) ~ 27.88733,
                                    grepl('NUASJARVI', lake) ~ 28.08605,
                                    lake == 'ANDERSON LAKE' ~ -89.34409,
                                    lake == 'FISHER LAKE' ~ -88.24577,
                                    lake == 'SPRUCE LAKE' ~ -89.56701,
                                    lake == 'MYSTERY LAKE' ~ -89.5655,
                                    lake == 'LAKE NIPISSING' ~ -79.76244,
                                    lake == '(BIG) LONG LAKE' ~ -91.65017,
                                    lake == 'BLACK OTTER LAKE' ~ -88.637,
                                    lake == 'LAKE NASIJARVI (3568)' ~ 23.75618,
                                    grepl('KITUSJARVI (3548)', lake) ~ 24.04439,
                                    lake == 'JAASJARVI - HARTOLA (1457)' ~ 26.15,
                                    lake == 'POROVESI' ~ 27.12763,
                                    lake == 'PALOVESI' ~ 23.93334,
                                    lake == 'KUIVAJARVI' ~ 23.86147,
                                    lake == 'DEEP' ~ -77.57017,
                                    lake == 'CHAUTAUQUA SOUTH' ~ -79.36177,
                                    lake == 'DUCK LAKE' ~ -84.78612,
                                    lake == 'FOUNTAIN' ~ -93.37752,
                                    lake == 'DETROIT LAKE' ~ -95.836,
                                    lake == 'LOWER PRIOR' ~ -93.41108,
                                    lake == 'PRIOR' ~ -93.44584,
                                    lake == 'LEECH' ~ -94.37896,
                                    lake == 'NIINIVESI W' ~ 26.75474,
                                    lake == 'KYNSIVESI W' ~ 26.25936,
                                    lake == 'LAKE KUBENSKOYE (PESKI)' ~ 39.53982,
                                    lake == 'LAKE B.LEPRINDO' ~ 117.50898,
                                    lake == 'LAKE NICHATKA' ~ 117.53266,
                                    lake == 'LAKE ASLI-KUL' ~ 54.57873,
                                    lake == 'LAKE EMANDGA' ~ 135.737,
                                    lake == 'LAKE SYABERO' ~ 29.15864,
                                    lake == 'LAKE ASLI-KUL' ~ 54.56775,
                                    lake == 'MINNOW LAKE' ~ -80.95655,
                                    lake == 'MARSH LAKE' ~ -134.33395,
                                    lake == 'ATLIN LAKE' ~ -135.07,
                                    lake == 'ANAHIM LAKE' ~ -125.34947,
                                    lake == 'ATLIN LAKE' ~ -133.75358,
                                    grepl('KENOGAMISIS', lake) ~ -86.86118,
                                    station == 'WRS206' ~ -86.92038,
                                    station == 'RB19' ~ 26.26149,
                                    lake == 'BEAR LAKE' ~ -118.93826,
                                    lake == 'LAC CACHE' ~ -74.41248,
                                    lake == 'CANOE LAKE' ~ -108.26149,
                                    lake == 'CHARLIE LAKE' ~ -120.95487,
                                    lake == 'LITTLE KETTLE LAKE' ~ -95.325759,
                                    lake == 'GILMAN LAKE' ~ -74.34697,
                                    lake == 'GLENMORE RESERVOIR ELBOW RIVER' ~ -114.12024,
                                    lake == 'GRAND LAKE' ~ -60.12862,
                                    lake == 'LAC MALFAIT' ~ -67.79,
                                    lake == 'LAKE ALEX' ~ -71.450989,
                                    lake == 'LAKE DORE' ~ -74.3763,
                                    lake == 'LAKE ONATCHIWAY'~ -71.05773,
                                    lake == 'LANDING LAKE' ~ -97.36118,
                                    station == 'WRS321' ~ -80.90402,
                                    lake == 'LYNN LAKE' ~ -101.04453,
                                    lake == 'MC LEOD LAKE' ~ -123.03985,
                                    lake == 'MEADOW LAKE' ~ -108.34765,
                                    lake == 'NUKKO LAKE' ~ -123.00959,
                                    lake == 'OBA LAKE' ~ -84.29181,
                                    lake == 'PLAYGREEN LAKE' ~ -98.13644,
                                    lake == 'PORTAGE BAY' ~ -98.79889,
                                    lake == 'UPPER TWIN LAKE AND LOWER TWIN LAKE' ~ -86.58052,
                                    lake == 'WILLIAMS LAKE' ~ -122.07128,
                                    TRUE ~ long))

# update stations
stations <- select(ice, station, lake, long, lat, altitude) %>%
  filter(!duplicated(station))

# check if there are any distinct lakes with the same name
if(VIEW.COMMENTS) {
  filter(stations,
         duplicated(lake) |
           duplicated(lake, fromLast = TRUE)) %>%
    arrange(lake) %>%
    View()
}

# distinguish between different lakes
ice <- mutate(ice,
              lake = case_when(station == 'WRS228' ~ 'CROOKED LAKE (BRITISH COLUMBIA)',
                               station == 'GAH6' ~ 'EAGLE LAKE (MAINE)',
                               station == 'MINN34' ~ 'GREEN LAKE (MINNESOTA)',
                               station == 'MICH03' ~ 'GULL LAKE (MICHIGAN)',
                               station == 'WRS388' ~ 'SWAN LAKE (BRITISH COLUMBIA)',
                               station == 'WRS395' ~ 'TROUT LAKE (ONTARIO)',
                               TRUE ~ lake))

# check comments ####
if(VIEW.COMMENTS) {
  # date changed by author
  View(filter(ice, comments == 'calendar correction for ice_on: -30 days of original data; applied on 1 February 2012'))
  
  # ignored the following comments because both stations were kept in the data
  View(filter(ice, comments %in% c('preferred is date from weather station: 5 January 1953 (J.Magnuson 1 Feb 2012)',
                                   'preferred is date from weather station: 30 December 1976 (J.Magnuson 1 Feb 2012)')))
  
  # discontinuous ice cover
  View(filter(ice, comments == 'open Jan 4 and close Jan 10'))
  View(filter(ice, comments == 'open Dec 15 and close Dec 17'))
  View(filter(ice, comments == '0 Interim Freeze Thaw Cycles'))    # (continuous ice cover)
  View(filter(ice, comments %in% c('1 Interim Freeze Thaw Cycles', # duration is less than Off.DOY - On.DOY, so no issues
                                   '2 Interim Freeze Thaw Cycles',
                                   '3 Interim Freeze Thaw Cycles',
                                   '4 Interim Freeze Thaw Cycles')))
  
  # discrepancies
  View(filter(ice, grepl('hardcopies at the library of Paul Smith College', comments))) # not an issue
  View(filter(ice, comments == 'note some discrepancies with dataset provided by Curt Stager CS1')) # not an issue: included both stations
  
  # missing date causes right-censored duration
  View(filter(ice, comments %in% c('duration over 36 days',
                                   'duration over 40 days',
                                   'duration over 49 days',
                                   'duration over 60 days',
                                   'duration over 72 days',
                                   'duration over 90 days')))
  
  # no ice
  View(filter(ice, comments == 'no complete freezing')) 
}

# other issues ####
# error: thawed before freezing (comment is NA)
if(VIEW.COMMENTS) view(ice[which(ice$Off.date < ice$On.date), ])

# check for other errors or problems
sum(ice$Off.On < ice$duration, na.rm = TRUE) # nrows w duration > ice-on period
sum(ice$Off.On < ice$duration, na.rm = TRUE) / nrow(ice) # percent

hist(ice$Off.On.Duration[ice$Off.On < ice$duration], breaks = 60,
     xlab = 'Off - On - duration',
     main = 'Number of Events with Duration > off.DOY - On.DOY')

# save data ####
# keep only necessary columns
ice <- select(ice,
              lake, station, Year, july.year, froze.bool, On.date, Off.date,
              On.DOY, Off.DOY, On.DOY.jul, Off.DOY.oct, duration, country, long,
              lat, altitude) %>%
  mutate(lake = factor(lake),
         station = factor(station))

# fixed: no unamed lakes
sum(grepl('UNNAMED', ice$lake))

saveRDS(ice, here('data/lake-ice-data.rds'))
