#' Install necessary packages for the project
#'
#'
#' @title installPackages
#' @return
#' @export
#' @author Firas Sadiyah


installPackages <- function() {

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Update all installed CRAN packages
update.packages(ask = FALSE)

# Update all installed Bioconductor packages
BiocManager::install(ask = FALSE)

# Install CRAN package if not already installed
installCran <- function(package) {

  # First test if the package is already installed
  if (!package %in% rownames(installed.packages())) {
    install.packages(package)
  }else {
    print(paste0(package, " is already installed"))
  }

}

# Install Bioconductor package if not already installed
installBioconductor <- function(package) {

  # First test if the package is already installed
  if (!package %in% rownames(installed.packages())) {
    BiocManager::install(package, update = FALSE)
  }else {
    print(paste0(package, " is already installed"))
  }

}

cranPackages <- c(
                  "data.table",
                  "devtools",
                  "eulerr",
                  "extrafont",
                  "ggridges",
                  "ggseqlogo",
                  "gt",
                  "here",
                  "languageserver",
                  "microbenchmark",
                  "RColorBrewer",
                  "remotes",
                  "scales",
                  "tidyverse",
                  "VennDiagram",
                  "statmod",
                  "magick",
                  "bookdown",
                  "miniUI",
                  "RSQLite"
                  )

lapply(cranPackages, installCran)

biocPackages <- c(
                  "BSgenome.Mmusculus.UCSC.mm9",
                  "ChIPseeker",
                  "clusterProfiler",
                  "csaw",
                  "edgeR",
                  "EnhancedVolcano",
                  "EnrichedHeatmap",
                  "MotifDb",
                  "org.Mm.eg.db",
                  "PCAtools",
                  "rGREAT",
                  "rtracklayer",
                  "TxDb.Mmusculus.UCSC.mm9.knownGene",
                  "regioneR",
                  "karyoploteR"
                  )

lapply(biocPackages, installBioconductor)

}
