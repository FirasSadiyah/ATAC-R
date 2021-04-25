#' Draw a scatter plot of the differential regions
#'
#'
#' @title scatterPlot
#' @param DT a \code{data.table} object; the table of ATAC-seq peaks data
#' @param a a string; the name of contrast 1
#' @param b a string; the name of contrast 2
#' @param sig_col a string; the value of the colour of significant regions
#' @param ns_col a string; the value of the colour of non-significant regions
#' @return ggplot2 object
#' @import ggplot2
#' @export
#' @author Firas Sadiyah


scatterPlot <- function(DT, a, b, sig_col, ns_col) {

    paste0(a, "_FDR") -> a_FDR
    paste0(a, "_LFC") -> a_LFC
    paste0(b, "_FDR") -> b_FDR
    paste0(b, "_LFC") -> b_LFC

    format(round(rsq(DT[[a_LFC]], DT[[b_LFC]]), 2), nsmall = 2) -> rsquared
    round(cor(DT[[a_LFC]], DT[[b_LFC]], use = "pairwise.complete.obs", method = "spearman"), 2) -> rho

    ggplot(DT, aes(x = DT[[a_LFC]], y = DT[[b_LFC]])) +
        ggrastr::rasterise(geom_point(alpha = 0.1,
                                      colour = ifelse(DT[[a_FDR]] < parList$FDR
                                                    & DT[[b_FDR]] < parList$FDR, sig_col, ns_col)), dpi = 300) +
        geom_vline(xintercept =  1, linetype = 3, alpha = 0.6) +
        geom_vline(xintercept = -1, linetype = 3, alpha = 0.6) +
        geom_hline(yintercept =  1, linetype = 3, alpha = 0.6) +
        geom_hline(yintercept = -1, linetype = 3, alpha = 0.6) +
        ylim(c(-10, 10)) +
        xlim(c(-8,  8))  +
        ggtitle("") +
        xlab(a_LFC) +
        ylab(b_LFC) +
        ggplot2::annotate(geom = "text", label = bquote(rho == .(rho)),
                          x = Inf, y = max(DT[[a_LFC]]),
                          color = "black", size = 4, fontface = "plain",
                          hjust = 1, vjust = 1) +
        themeFX() -> gg

        print(gg)
        dev.copy(pdf, file = here::here("Plots", paste0("scatterPlot_", a, "_vs_", b, ".pdf")))
        dev.off()
}
