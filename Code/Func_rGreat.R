#' Perform Gene Ontology enrichment analysis using \code{rGREAT}
#'
#'
#' @title doGreat
#' @param contrast a string; the contrast to analyse
#' @param peakSet a string; the set of ATAC-seq peaks
#' @return a list of rGREAT results for different ontologies
#' @import rGREAT
#' @import data.table
#' @export
#' @author Firas Sadiyah


doGreat <- function(contrast, peakSet) {

    Great <- list()

    message("Getting enrichment for opened regions")
	Great$jobUP <- suppressMessages(submitGreatJob(makeGRangesFromDataFrame(contrast[["DTs"]][["DT"]][peakSet$UP,
	                                                                    	c("chrom", "start", "end", "peaks")]), species = "mm9"))
    Great$goUP$MF <- suppressMessages(getEnrichmentTables(Great$jobUP, category = "GO")[["GO Molecular Function"]])
    Great$goUP$BP <- suppressMessages(getEnrichmentTables(Great$jobUP, category = "GO")[["GO Biological Process"]])
    Great$goUP$CC <- suppressMessages(getEnrichmentTables(Great$jobUP, category = "GO")[["GO Cellular Component"]])

    for (i in c("MF", "BP", "CC")) {
        x <- as.data.table(Great$goUP[[i]])
        xSig <- x[x$Binom_Adjp_BH < parList$FDR, ]
        xSig <- xSig[order(rank(Binom_Adjp_BH))]
        xSig$name <- factor(xSig$name, levels = rev(unique(xSig$name)))
        i_Sig <- paste0(i, ".Sig")
        Great$goUP[i_Sig] <- list(xSig)
    }

    message("Getting enrichment for closed regions")
	Great$jobDN <- suppressMessages(submitGreatJob(makeGRangesFromDataFrame(contrast[["DTs"]][["DT"]][peakSet$DN,
		                                                                    c("chrom", "start", "end", "peaks")]), species = "mm9"))
    Great$goDN$MF <- suppressMessages(getEnrichmentTables(Great$jobDN, category = "GO")[["GO Molecular Function"]])
    Great$goDN$BP <- suppressMessages(getEnrichmentTables(Great$jobDN, category = "GO")[["GO Biological Process"]])
    Great$goDN$CC <- suppressMessages(getEnrichmentTables(Great$jobDN, category = "GO")[["GO Cellular Component"]])

    for (i in c("MF", "BP", "CC")) {
        x <- as.data.table(Great$goDN[[i]])
        xSig <- x[x$Binom_Adjp_BH < parList$FDR, ]
        xSig <- xSig[order(rank(Binom_Adjp_BH))]
        xSig$name <- factor(xSig$name, levels = rev(unique(xSig$name)))
        i_Sig <- paste0(i, ".Sig")
        Great$goDN[i_Sig] <- list(xSig)
    }

    return(Great)

}

#' Draw a bar plot for the top entries of \code{rGREAT} analysis
#'
#'
#' @title plotrGreat
#' @param data a \code{data.table} object; a table of ATAC-seq peaks data
#' @param ylabel a string; the title of the y axis
#' @param n an integer; the number of entries to include in the plot
#' @return ggplot2 object
#' @import ggplot2
#' @export
#' @author Firas Sadiyah


plotGreat <-  function(data, ylabel, n=10) {
    ggplot(data = data[1:n, ], aes(x = name, y = (-log10(Binom_Adjp_BH)))) +
       geom_bar(stat = "identity", fill = "steelblue") +
       theme_bw() +
       theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) +
       coord_flip() + labs(y = bquote(-log[10] ~ "adj" ~ italic(P) ~ "values"), x = ylabel) +
       geom_hline(yintercept = -log10(parList$FDR), linetype = "dashed") +
    themeFX() +
    theme(axis.text.y = element_text(hjust = 0.95, vjust = 0.2)) -> gg

    name <- gsub("\\$", "_", deparse(substitute(data)))
    name <- gsub("_*.Great", "", name)

    print(gg)
    dev.copy(pdf, file = here::here("Plots", paste0("rGreat_", name, ".pdf")))
    dev.off()
}
