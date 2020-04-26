# find n days after the reference date for an entire dataset
post.ref.date <- function(column, d, event = c('freeze', 'thaw')) {
  
  if(length(event) > 1 | (event[1] != 'freeze' & event[1] != 'thaw'))
    stop('Please choose an event, either "freeze" or "thaw".')
  
  ref.date <- case_when(event == 'freeze' ~ '-06-30',
                        event == 'thaw' ~ '-09-30')
  
  # case_when and if_else cannot deal with NAs
  d[is.na(d[[column]]), column] <- as.Date('0001-01-01')
  
  d <- mutate(d,
              Date = d[[column]],
              Year = year(Date),
              DOY = yday(Date),
              NA.DOY = if_else(Date == as.Date('0001-01-01'), TRUE, FALSE),
              
              # leap years have 366 days, but if year starts after Februart then the years *before* have 366 days
              n.days.prev.year = yday(paste0(Year - 1, '-12-31')),
              ref.date.prev.year = yday(paste0(year(Date) - 1, ref.date)),
              new.DOY = if_else(yday(Date) < yday(paste0(Year, ref.date)),        # if before ref.date
                                DOY + n.days.prev.year - ref.date.prev.year,      # DOY + (365 or 366) - ref.date
                                yday(Date) - yday(paste0(year(Date), ref.date)))) # ref.date: subtract yday(ref.date)
  
  # change '0001-01-01' back to NA
  d$new.DOY <- case_when(d$NA.DOY == TRUE ~ NA_real_,
                         TRUE ~ d$new.DOY)
  
  d$new.DOY
}
