#' Construct a \code{ComplexHeatmap} object from count matrix and BigWig tracks
#'
#'
#' @title heatmapPlot
#' @param DT a \code{data.table} object; the table of ATAC-seq peaks data
#' @param hm a list; list of paramters of the heatmap, e.g. annotation, unit, width, order
#' @param rhm a list; list of paramters of the specific region, e.g. number of cluster
#' @param enriched a Boolean; whether to include enriched heatmaps
#' @param dend a Boolean; whether to include a dendrogram
#' @return a \{ComplexHeatmap} object
#' @import data.table
#' @import ComplexHeatmap
#' @import EnrichedHeatmap
#' @import GenomicRanges
#' @import dendextend
#' @export
#' @author Firas Sadiyah


heatmapPlot <- function(DT, hm, rhm, enriched = TRUE, dend = FALSE) {

    # set seed for reproducibility
    set.seed(1234)

    # start tracking the duration of the function
    timeStart <- Sys.time()

    # shorten contrast list handler
    cl <- hm$contList
    # shorten regions handler
    regions <- rhm$regions

    message("Aggregating ", nrow(regions), " peaks into ", rhm$kl, " clusters...")

    # compute similarities
    cor(t(regions[, ..cl]), method = "spearman", use = "complete.obs") -> i_cor
    # compute dissimilarity
    as.dist(1 - i_cor) -> i_dist
    # cluster distances
    hclust(i_dist, method = "ward.D") -> i_hclust
    # cut the tree at the specified kl number
    cutree(i_hclust, k = rhm$kl) -> i_clusters
    # add the clusters to the data.table
    regions$cluster <- i_clusters
    # z-scale the normalised counts
    t(scale(t(regions[, ..cl]))) -> regions.mat

    # Prepare the genomic ranges of target regions
    makeGRangesFromDataFrame(as.data.frame(regions[, c("chrom", "start", "end", "peaks", "cluster")]),
                                           keep.extra.columns = TRUE) -> targets
    resize(targets, fix = "center", width = parList$bwExt * 2) -> targets.ext

    # Pickup colours for the clusters
    structure(colPal1(n = (length(unique(regions$cluster)))),
                   names = seq_along(unique(regions$cluster))) -> colClst

    enrichedPlots <- list()
    normMatrices  <- list()


    rowOrder <- NULL

    if (enriched == TRUE) {

    # Calculate the quantil of the reference bigwig to define colour ranges
    message("Calculating bigwig quantil based on ", cl[hm$orderBy], " ...")
    rtracklayer::import(con = paste0(parList$bwPath, cl[hm$orderBy], ".bw"),
                            format = "BigWig", selection = BigWigSelection(targets.ext)) -> ref_bw
    normalizeToMatrix(signal = ref_bw, target = resize(targets, fix = "center", width = 1),
                          extend = parList$bwExt, background = 0, keep = c(0, 0.99), value_column = "score") -> ref_nm
    enrMapCol <- circlize::colorRamp2(quantile(ref_nm, c(0, .90, .99)), hm$eCol)


    # read each bigwig file, construct into normalised matrix, and plot enriched heatmap

    # Load the content of the bigwig limited to the regions we are interested in.
    for (contrast in cl) {
        message(paste0("Processing ", contrast, "..."))
        rtracklayer::import(con = paste0(parList$bwPath, contrast, ".bw"),
                            format = "BigWig", selection = BigWigSelection(targets.ext)) -> i_bw
        normalizeToMatrix(signal = i_bw, target = resize(targets, fix = "center", width = 1),
                          extend = parList$bwExt, background = 0, keep = c(0, 0.99), value_column = "score") -> i_nm

        normMatrices[[contrast]]  <- list(i_nm)

        EnrichedHeatmap(i_nm,
                        col = enrMapCol,
                        name = contrast,
                        column_title = contrast,
                        use_raster = TRUE,
                        raster_device = "CairoPNG",
                        raster_quality = 5,
                        axis_name = c("-5k", "0", "5k"),
                        pos_line = FALSE,
                        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = colClst, lty = 1),
                                                           axis_param = list(side = "right", facing = "inside"),
                                                           ylim = c(0, 6.5))
                                                           ),
                        width = unit(2, "cm")
                        ) -> i_plot

       enrichedPlots[[contrast]] <- list(i_plot)
       }
       rowOrder = order(enriched_score(normMatrices[[hm$orderBy]][[1]]), decreasing = TRUE)
    }

    # build annodation for the heatmap
    hm.meta <- data.frame(Time = hm$Time, Cell = hm$Cell, Genotype = hm$Genotype)
    hm.col <- list("Cell"     = c("Trg" = colList$trg, "Th0" = colList$th0),
                   "Time"     = c("00h" = colList$T00, "06h" = colList$T06,
                                  "18h" = colList$T18, "72h" = colList$T72),
                   "Genotype" = c("WT"  = colList$wt, "KO"  = colList$ko))
    hm.anno <- HeatmapAnnotation(df = hm.meta, which = "col", col = hm.col,
                                annotation_width = unit(c(1, 4), "cm"),
                                annotation_name_side = "left"
                                )

    # build cluster annotation (as heatmap)
    message(paste0("Plotting cluster heatmap..."))
    ComplexHeatmap::Heatmap(targets$cluster,
           col = colClst,
           row_order =  rowOrder,
           use_raster = TRUE,
           cluster_columns = FALSE,
           cluster_rows = FALSE,
           name = "Cluster",
           show_row_names = FALSE,
           width = unit(3, "mm")) -> cluster.heatmap

    message("Plotting main heatmap...")

    # build main heatmap
    ComplexHeatmap::Heatmap(regions.mat,
                            cluster_columns = FALSE,
                            cluster_rows = FALSE,
                            row_split = (factor(targets$cluster, levels = 1:rhm$kl)),
                            column_split = hm$colSplit,
                            row_dend_reorder = FALSE,
                            show_row_names = FALSE,
                            top_annotation = hm.anno,
                            name = "Log2FC",
                            col = hm$hCol,
                            width = unit(5, "cm"),
                            border = TRUE,
                            use_raster = TRUE,
                            ) -> main.heatmap


    timeEnd <- Sys.time()
    duration <- difftime(timeEnd, timeStart, units = "mins")

    message(paste0("Done in ", round(duration, 2), " mins."))

    if (dend == TRUE) {
        message("Plotting the dendrogram of hierarchial clustering...")
        plot(as.dendrogram(i_hclust))
    }

    invisible(gc())

    maps <- list(cluster.heatmap, main.heatmap, unlist(enrichedPlots), regions$cluster)
    return(maps)

}

#' draw a \code{ComplexHeatmap} object constructed by \code{heatmapPlot}
#'
#'
#' @title heatmapPlot
#' @param hm a list of paramters of the heatmap, e.g. annotation, unit, width, order
#' @param rhm a list of paramters of the specific region, e.g. number of cluster
#' @param enriched a Boolean; whether to include enriched heatmaps
#' @return a \{ComplexHeatmap} object
#' @import data.table
#' @import ComplexHeatmap
#' @import EnrichedHeatmap
#' @export
#' @author Firas Sadiyah

drawHeatmap <- function(hm, rhm, enriched = TRUE) {

    # set seed for reproducibility
    set.seed(1234)
    ht_gap <- unit(2, "mm")

    # construct the list of heatmaps to plot
    drawList <- NULL
    # list the cluster and main heatmap
    drawList <- rhm$maps[[1]] + rhm$maps[[2]]

    if (enriched == TRUE) {
    # list each of the enirched heatmaps
    for (map in seq_along(hm$meta$contList)) {
        drawList = drawList + rhm$maps[[3]][[map]]
        }
    ht_gap <- hm$meta$gaps
    }



   message("Drawing combined heatmap...")

   draw(drawList,
        split = factor(rhm$maps[[4]], levels = 1:rhm$kl),
        column_title_gp = gpar(fontsize = 12, col = "black"),
        row_title_gp    = gpar(fontsize = 12, col = "black"),
        annotation_legend_side = "bottom",
        heatmap_legend_side    = "bottom",
        padding = unit(c(4, 4, 4, 4), "mm"),
        ht_gap = ht_gap
        ) -> drawObj

   print(drawObj)
   dev.copy(pdf, file = here::here("Plots", paste0("heatMaps_", rhm$name, ".pdf")), height = 12, width = 10)
   dev.off()

   return(drawObj)
}
