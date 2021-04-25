#' Write a genomic ranges of specific regions as .bed file
#'
#'
#' @title writeFasta
#' @param contrast a string; the contrast to use as a source for the genomic ranges
#' @param regions a string; the regions to include
#' @param direction a string; the direction of change, e.g. UP, DN, UD
#' @return .fa file
#' @import data.table
#' @import GenomicRanges
#' @import Biostrings
#' @import BSgenome.Mmusculus.UCSC.mm9
#' @export
#' @author Firas Sadiyah


writeBed <- function(contrast, regions, direction = NULL) {
    regions <- deparse(substitute(regions))

    if (!missing(direction)) {
        direction <- deparse(substitute(direction))
        write.table(
            contrast$DTs$DT[contrast$peaks[[regions]][[direction]], c("chrom", "start", "end", "peaks")],
            here::here(
                "Out", "bed",
                paste0(deparse(substitute(contrast)), "_", regions, "_", direction, ".bed")
            ),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
        )
    } else {
        if (is.null(contrast$Heatmap[[regions]]$drawObj)) {
            message("Order of clusters not found, please run corresponding heatmap first.")
        } else {
            clusters <- ComplexHeatmap::row_order(contrast$Heatmap[[regions]]$drawObj)
            for (cl in seq_along(clusters)) {
                write.table(
                    contrast$DTs$DT[clusters[[cl]], c("chrom", "start", "end", "peaks")],
                    here::here(
                        "Out", "bed",
                        paste0(deparse(substitute(contrast)), "_", regions, "_C", cl, ".bed")
                    ),
                    quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
                )
            }
        }
    }
}
