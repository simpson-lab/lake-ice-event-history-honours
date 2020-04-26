save.plt <- function(plt = last_plot(), dir, width = 8, height = 6, dpi = 300) {
  ggsave(filename = paste0('plots/', dir), plot = plt,
         width = width, height = height, units = 'in', dpi = 300)
}
