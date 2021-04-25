#' Extract different sets from a venn diagram
#'
#'
#' @title extractVennSets
#' @param venn a \code{VennDiagram} object
#' @return a list of different sets in the input \code{VennDiagram} object
#' @import VennDiagram
#' @export
#' @author Firas Sadiyah


extractVennSets <- function(venn) {

       VennDiagram::get.venn.partitions(venn, keep.elements = TRUE) -> venn.sets

       get.venn.partitions(venn)$..set.. -> setsNames

       setsList <- list()

       for (set in seq_along(setsNames)) {
              i_set  <- venn.sets$..values..[set]
              setsList$set <- c(i_set, atacPeaks$seq[unlist(i_set)])
              names(setsList)[set] <- setsNames[set]

       }

       return(setsList)

}
