#' A custom theme for the plots of the project
#'
#'
#' @title themeFX
#' @return ggplot object
#' @import ggplot2
#' @export
#' @author Firas Sadiyah


themeFX <- function() {
    theme_classic() +
        theme(
              axis.text.x           = element_text(angle  = 00, vjust = 0.5, hjust = 0.5, family = "Helvetica", size = 12),
              axis.title.x.bottom   = element_text(family = "Helvetica", size = 12),
              axis.text.y           = element_text(angle  = 00, vjust = 0.5, hjust = 0.5, family = "Helvetica", size = 12),
              axis.title.x          = element_text(colour = "Black", family = "Helvetica", size = 12, margin = margin(t = 20, r = 0 , b = 0 , l = 0)),
              axis.title.y          = element_text(colour = "Black", family = "Helvetica", size = 12, margin = margin(t =  0, r = 20, b = 0 , l = 0)),
              plot.title            = element_text(colour = "Black", family = "Helvetica", size = 12, margin = margin(t = 20, r = 0 , b = 20, l = 0)),
              legend.text           = element_text(family = "Helvetica", size = 12),
              legend.title          = element_text(family = "Helvetica", size = 12),
              #axis.line             = element_blank(),
              strip.text            = element_text(family = "Helvetica", size = 12),
              strip.background      = element_blank(),
              panel.background      = element_rect(fill   = "transparent", colour = NA), # bg of the panel
              plot.background       = element_rect(fill   = "transparent", colour = NA), # bg of the plot
              panel.grid.major      = element_blank(),                                   # no major grid
              panel.grid.minor      = element_blank(),                                   # no minor grid
              legend.background     = element_rect(fill   = "transparent", colour = NA), # no legend bg
              legend.box.background = element_rect(fill   = "transparent", colour = NA)  # no legend panel bg
              )
}
