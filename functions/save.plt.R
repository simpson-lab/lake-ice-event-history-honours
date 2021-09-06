save.plt <- function(plt = ggplot2::last_plot(), dir, width = 8, height = 6,
                     dpi = 300, scale = 1) {
  ggsave(filename = paste0('plots/', dir), plot = plt,
         width = width, height = height, units = 'in', dpi = 300, scale = scale)
}
