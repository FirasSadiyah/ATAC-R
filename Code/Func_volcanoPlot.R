#' Draw a volcano plot for each contrast and label specifc set of genes
#'
#'
#' @title plotVolcano
#' @param DT a \code{data.table} object; the table of ATAC-seq peaks data
#' @param a a string; the name of contrast to include in the plot
#' @param opened_col a string; the value of the colour of opened regions
#' @param closed_col a string; the value of the colour of closed regions
#' @param ns_col a string; the value of the colour of non-significant regions
#' @param ymax an integer; the upper limit of the y axis
#' @param labList a list; the list of genes to label on the plot
#' @return ggplot2 object
#' @import data.table
#' @import ggplot2
#' @export
#' @author Firas Sadiyah


plotVolcano <- function(DT, a, opened_col, closed_col, ns_col, ymax=50, labList = NULL) {

    paste0(a, "_FDR") -> a_FDR
    paste0(a, "_LFC") -> a_LFC

          keyvals  <- rep(ns_col, nrow(DT))
    names(keyvals) <- rep("NS", nrow(DT))
           keyvals[which(DT[[a_FDR]] < parList$FDR & DT[[a_LFC]] >  parList$LFC)] <- opened_col
    names(keyvals)[which(DT[[a_FDR]] < parList$FDR & DT[[a_LFC]] >  parList$LFC)] <- "Opened"
           keyvals[which(DT[[a_FDR]] < parList$FDR & DT[[a_LFC]] < -parList$LFC)] <- closed_col
    names(keyvals)[which(DT[[a_FDR]] < parList$FDR & DT[[a_LFC]] < -parList$LFC)] <- "Closed"

    EnhancedVolcano(DT,
                    lab = NA, x = a_LFC, y = a_FDR,
                    #selectLab = rownames(DT)[which(names(keyvals) %in% c("Opened","Closed"))][1:10],
                    selectLab = labList,
                    pCutoff = parList$FDR, FCcutoff = parList$LFC, pointSize = 2.0, labSize = 3.0, colAlpha = 0.8,
                    drawConnectors = FALSE, widthConnectors = 0.75, colConnectors = "grey40",
                    xlab = bquote(~Log[2]~ "Fold Change"),
                    ylab = bquote(~-Log[10]~adjusted~italic(P)),
                    title = "",
                    subtitle = a,
                    gridlines.major = FALSE, gridlines.minor = FALSE,
                    border = "full",
                    borderWidth = 1.0, borderColour = "black", xlim = c(-4, 4), ylim = c(0, ymax),
                    legendPosition = "bottom", legendLabSize = 9, legendIconSize = 3,
                    colCustom = keyvals,
                    raster = TRUE
    ) -> gg

    print(gg)
    dev.copy(pdf, file = here::here("Plots", paste0("volcanoPlot_", a, ".pdf")))
    dev.off()

}
