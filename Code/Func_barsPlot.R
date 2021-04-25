#' Plot the number of peaks in opened or closed regions
#'
#'
#' @title barsPlot
#' @param DT a \code{data.table} object; the table of ATAC-seq peaks data
#' @param contrastList a list; the list of contrasts to include
#' @param opened_col a string; the value of the colour of opened regions
#' @param closed_col a string; the value of the colour of closed regions
#' @param fileName a string; the name of the exported plot
#' @param n an inter; the number of contrasts to compare
#' @return ggplot2 object
#' @import data.table
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_text
#' @export
#' @author Firas Sadiyah


barsPlot <- function(DT, contrastList, opened_col, closed_col, fileName, n = 3) {

    DT_UP <- list()
    DT_DN <- list()

    for (contrast in contrastList) {
       paste0(contrast, "_FDR") -> contrast_FDR
       paste0(contrast, "_LFC") -> contrast_LFC
       DT_UP[contrast] <- DT[DT[[contrast_FDR]] < parList$FDR & DT[[contrast_LFC]] >  parList$LFC, .N]
       DT_DN[contrast] <- DT[DT[[contrast_FDR]] < parList$FDR & DT[[contrast_LFC]] < -parList$LFC, .N]
    }

    print(DT_UP)
    print(DT_DN)

    if (n == 3) {
    conditions <- factor(rep(c("6", "18", "72"), 2), levels = c(6, 18, 72))
    amount <- c(unlist(DT_UP), unlist(DT_DN))
    change <- c(rep("Opened", 3), rep("Closed", 3))
    DT_bars <- data.frame(conditions, change, amount)
    } else if (n == 4) {
    conditions = factor(rep(c("0", "6", "18", "72"), 2), levels = c(0, 6, 18, 72))
    amount = c(unlist(DT_UP), unlist(DT_DN))
    change = c(rep("Opened", 4), rep("Closed", 4))
    DT_bars = data.frame(conditions, change, amount)
    }

    colnames(DT_bars) <- c("condition", "change", "amount")
    print(DT_bars)
    ggplot(DT_bars, aes(x = condition, fill = change)) +
        geom_bar(data = subset(DT_bars, change == "Opened"),
                 aes(y = amount),
                 position = "stack",
                 stat = "identity") +
        geom_text(data = subset(DT_bars, change == "Opened"),
                  aes(y = amount, label = amount),
                  position = position_dodge(width = 1),
                  vjust = -0.5) +
        geom_bar(data = subset(DT_bars, change == "Closed"),
                 aes(y = -amount),
                 position = "stack",
                 stat = "identity") +
        geom_text(data = subset(DT_bars, change == "Closed"),
                  aes(y = -amount, label = -amount),
                  position = position_dodge(width = 1),
                  vjust = +1.5) +
        ggtitle("", subtitle = "") +
        ylab("Peaks") +
        xlab("Time points") +
        ylim(c(-5500, 5500)) +
        scale_fill_manual(name = "Change",
                          values = c("Opened" = opened_col, "Closed" = closed_col),
                          breaks = c("Opened", "Closed")) +
        themeFX() -> gg

        print(gg)
        dev.copy(pdf, file = here::here("Plots", paste0("barsPlot_", fileName, ".pdf")))
        dev.off()
}
