#' Draw a venn plot of the differential regions
#'
#'
#' @title plotVolcano
#' @param x a list; the list of overlapped regions
#' @param cols a list; the list of colours to be used in the venn plot
#' @param fileName a string; the name of the exported plot
#' @return ggplot2 object
#' @import ggplot2
#' @import eulerr
#' @export
#' @author Firas Sadiyah


vennPlot <- function(x, cols, fileName) {
   plot(euler(x),
   fills = list(fill = cols, alpha = 0.50),
           labels = list(col = "black", font = 4),
           legend = list(col = "black", font = 4),
           main = "",
           quantities = TRUE,
           shape = "ellipse",
           lty = 1) -> gg

    print(gg)
    dev.copy(pdf, file = here::here("Plots", paste0("vennPlot_", fileName, ".pdf")))
    dev.off()
}
