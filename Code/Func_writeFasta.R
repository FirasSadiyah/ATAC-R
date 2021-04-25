#' Extract fasta sequence of genomic ranges and write to a .fa file
#'
#'
#' @title writeFasta
#' @param contrast a string; the contrast to use as a source for the genomic ranges
#' @param regions a string; the regions to include
#' @param direction a string; the direction of change, e.g. UP, DN, UD
#' @param byLength an integer; the length of extension to apply to each genomic range
#' @return .fa file
#' @import data.table
#' @import GenomicRanges
#' @import Biostrings
#' @import BSgenome.Mmusculus.UCSC.mm9
#' @export
#' @author Firas Sadiyah


writeFasta <- function(contrast, regions, direction = NULL, byLength = 250) {

    regions <- deparse(substitute(regions))

    if (!missing(direction)) {
        direction <- deparse(substitute(direction))
        Biostrings::getSeq(
            BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9,
            GenomicRanges::resize(GenomicRanges::makeGRangesFromDataFrame(
                contrast$DTs$DT[contrast$peaks[[regions]][[direction]], c("chrom", "start", "end", "peaks")]
            ),
            byLength,
            fix = "center"
            )
        ) -> x

        names(x) <- paste0("peak_", contrast$peaks[[regions]][[direction]])

        Biostrings::writeXStringSet(x, here::here(
            "Out", "fasta",
            paste0(deparse(substitute(contrast)), "_", regions, "_", direction, "_", byLength, "bp", ".fa")
        ))
    } else {
        if (is.null(contrast$Heatmap[[regions]]$drawObj)) {
            message("Order of clusters not found, please run corresponding heatmap first.")
        } else {
            clusters <- ComplexHeatmap::row_order(contrast$Heatmap[[regions]]$drawObj)
            for (cl in seq_along(clusters)) {
                Biostrings::getSeq(
                    BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9,
                    GenomicRanges::resize(GenomicRanges::makeGRangesFromDataFrame(
                        contrast$DTs$DT[clusters[[cl]], c("chrom", "start", "end", "peaks")]
                    ),
                    byLength,
                    fix = "center"
                    )
                ) -> x

                names(x) <- paste0("peak_", clusters[[cl]])

                Biostrings::writeXStringSet(x, here::here(
                    "Out", "fasta",
                    paste0(deparse(substitute(contrast)), "_", regions, "_CL", cl, "_", byLength, "bp", ".fa")
                ))
            }
        }
    }
}
