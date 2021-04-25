#' Draw a line plot for each contrast with facet for each genomic annotation
#'
#'
#' @title linesPlot
#' @param DT a \code{data.table} object; the table of ATAC-seq peaks data
#' @param cl a list; list of contrast to include in the plot
#' @param opened_col a string; a value of the colour of opened regions
#' @param closed_col a string; a value of the colour of closed regions
#' @param x an integer; a value representing the index of the contrast in \code{cl}
#' @return ggplot2 object
#' @import data.table
#' @import ggplot2
#' @export
#' @author Firas Sadiyah


linesPlot <- function(DT, cl, opend_col, closed_col, x) {

    cn <- list()
    for (i in cl) {
        cn[i] <- list(i = c(FDR = paste0(i, "_FDR"), LFC = paste0(i, "_LFC"), DAR = paste0(i, "_DAR")))
    }

    sel <- function(x, y) {cn[[x]][[y]]}

    subset <- DT[DT[[sel(x, "DAR")]] %in% c("Opened", "Closed")]

    lines <- melt(subset,
                 id.vars       = c("peaks"),
                 measure.vars  = c(sel(1, "LFC"),
                                   sel(2, "LFC"),
                                   sel(3, "LFC")),
                 variable.name = "Condition_FC",
                 value.name    = "log2FC"
    )

    lines[, condition := gsub("_LFC", "", Condition_FC)]
    lines[, condition := factor(condition, levels = c(cl[[1]], cl[[2]], cl[[3]]))]
    data.table::setkey(lines, peaks)

    # add in change info for DA status
    DT.DAR = DT[, list(peaks, DT[[sel(1, "DAR")]], DT[[sel(2, "DAR")]], DT[[sel(3, "DAR")]])]
    names(DT.DAR) <- c("peaks", sel(1, "DAR"), sel(2, "DAR"), sel(3, "DAR"))
    data.table::setkey(DT.DAR, peaks)

    lines <- DT.DAR[lines]
    lines_anno <- atacPeaks$annt[lines]

    ggplot(lines_anno, aes(x = condition, y = log2FC)) +
        ggrastr::rasterise(geom_line(inherit.aes = FALSE,
                                     aes(x = condition, y = log2FC, group = peaks),
                                     alpha = 0.1, size = .15, linetype = "dashed"), dpi = 300) +
        ggrastr::rasterise(geom_point(aes(colour = lines_anno[[sel(x, "DAR")]])), dpi = 300) +
        ggtitle(cn[[x]]["DAR"]) +
        xlab("Condition") +
        ylab(bquote(~log[2]~"FC")) +
        ylim(c(-10, 10)) +
        scale_colour_manual(name = "",
                            values = c("Closed" = closed_col, "Opened" = opend_col),
                            breaks = c("Opened", "Closed")
                            ) +
        facet_wrap(~ annotation, ncol = 3) +
        theme(panel.background = element_blank()) +
        themeFX() -> gg

    print(gg)
    dev.copy(pdf, file = here::here("Plots", paste0("linesPlot_", cn[[x]]["DAR"], ".pdf")), height = 12, width = 10)
   dev.off()
}
