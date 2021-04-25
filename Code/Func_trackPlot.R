#' Plot customised linear genome tracks plot
#'
#'
#' @title groupTracksPlot
#' @param x a string; the region to plot e.g. "chr2:11,563,187-11,565,830"
#' @param maxValue an integer; the height of each genomic track
#' @param tickLength an integer; the length of the tick
#' @param tracksHeight an integer; the height of all genomic tracks
#' @param tracksTop an integer; the value of the ceiling of each track
#' @param tracksBottom an integer; the value of the floor of each track
#' @param rect.Col a string; the value of the colour of annotation box
#' @param back.Col a string; the value of the background of each track
#' @param track.Col1 a list; the list colours of the tracks in group 1
#' @param track.Col2 a list; the list colours of the tracks in group 2
#' @param chromHMM a string; the chrom_HMM track for the specifci cell type
#' @param tracks1 a list; the list of bigwig tracks to include in group 1
#' @param tracks2 a list; the list of bigwig tracks to include in group 2
#' @param groupLabel1 a string; the label of group 1
#' @param groupLabel2 a string; the label of group 2
#' @param plotTitle a string; the plot title
#' @return ggplot2 object
#' @import karyoploteR
#' @import regioneR
#' @import ggplot2
#' @export
#' @author Firas Sadiyah




groupTracksPlot <- function(
                      x,
                      maxValue = 14,
                      tickLength = 80,
                      tracksHeight = 600,
                      tracksTop = .99,
                      tracksBottom = .17,
                      rect.Col = "#AAFFCB80",
                      back.Col = "#F5F5F580",
                      track.Col1 = colPal1(12)[c(7, 7, 7, 7, 8, 8, 8, 8)],
                      track.Col2 = colPal1(12)[c(4, 4, 4)],
                      chromHMM = chrom_HMM,
                      tracks1,
                      tracks2,
                      groupLabel1,
                      groupLabel2,
                      plotTitle) {

## region coordinates
zoom.region <- regioneR::toGRanges(x)

# adjust the plotting parameters - aesthetic
pp <- karyoploteR::getDefaultPlotParams(plot.type = 1)
pp$leftmargin <- .2
pp$rightmargin <- .2
pp$ideogramheight <- 0.8
pp$data1height <- tracksHeight

# start by plotting gene track
kp <- karyoploteR::plotKaryotype(zoom = zoom.region, genome = "mm9", cex = 1, plot.params = pp)

# create a gene.data structure for kpPlotGenes
genes.data <- suppressMessages(karyoploteR::makeGenesDataFromTxDb(TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene,
                                                                  karyoplot = kp,
                                                                  plot.transcripts = TRUE,
                                                                  plot.transcripts.structure = TRUE))

# add gene symbol and merge mRNA transcripts
genes.data <- karyoploteR::addGeneNames(genes.data)
genes.data <- karyoploteR::mergeTranscripts(genes.data)

# add the base numbers along the chromosome ideograms
karyoploteR::kpAddBaseNumbers(kp, tick.dist = 1e4, minor.tick.dist = 5e3, add.units = TRUE, cex = 1, tick.len = 8)

# Plot the genes in this region
karyoploteR::kpPlotGenes(kp, data = genes.data, r0 = .01, r1 = .06, gene.name.cex = 1)

# add rectangle to draw background
karyoploteR::kpRect(kp,
                    chr = seqnames(zoom.region), x0  = start(zoom.region), x1  = end(zoom.region),
                    y0  = tracksBottom + .03, y1  = tracksTop - .01, col = back.Col, data.panel = "all", border = NA)

# add rectangle to highlight regions
if (!is.null(y)) {

## rectangle coordinates
rect.region <- regioneR::toGRanges(y)
## draw rectangle
karyoploteR::kpRect(kp,
                    chr = seqnames(rect.region), x0  = start(rect.region), x1  = end(rect.region),
                    y0  = tracksBottom + .03, y1  = tracksTop - .01, col = rect.Col, data.panel = "all", border = NA)
}

# Plot main title

# Plot chromHMM bed with label
karyoploteR::kpPlotRegions(kp, chromHMM, col = chromHMM$itemRgb,  r0 = .09, r1 = .12)
karyoploteR::kpAddLabels(kp, labels  =  "Chromatin\nState (HMM)", r0 = .08, r1 = .12, cex = 1)

# total tracks
total.tracks <- length(tracks1) + length(tracks2)

# tracks #1
out1.at <- karyoploteR::autotrack(seq_along(tracks1), total.tracks, margin = .05, r0 = tracksBottom, r1 = tracksTop)
# group1 label
karyoploteR::kpAddLabels(kp, labels = groupLabel1, r0 = out1.at$r0, r1 = out1.at$r1, side = "left", cex  = 1,
                             srt = 90, pos = 1, label.margin = .18)
# plot individual bigwig files
for (i in seq_along(tracks1)) {
  at <- karyoploteR::autotrack(i, length(tracks1), r0 = out1.at$r0, r1 = out1.at$r1, margin = .1)
  kp <- karyoploteR::kpPlotBigWig(kp, data = tracks1[i], ymax = maxValue, #"visible.region",
                                  r0 = at$r0, r1 = at$r1, col = track.Col1[i], border = track.Col1[i])
  computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
  karyoploteR::kpAxis(kp, ymin = 0, ymax = computed.ymax, tick.pos = c(0, computed.ymax),
                      tick.len = tickLength, r0 = at$r0, r1 = at$r1, cex = 1)
  # individual bigwig label
  karyoploteR::kpAddLabels(kp, labels = names(tracks1)[i], r0 = at$r0, r1 = at$r1, cex = 1, label.margin = .03)
}

if (!is.null(tracks2)) {
# tracks #2
out2.at <- karyoploteR::autotrack((length(tracks2) + 1):total.tracks, total.tracks, margin = .3, r0 = tracksBottom, r1 = tracksTop)
# group2 label
karyoploteR::kpAddLabels(kp, labels = groupLabel2, r0 = out2.at$r0, r1 = out2.at$r1, side = "left", cex = 1,
                             srt = 90, pos = 1, label.margin = .18)

for (i in seq_along(tracks2)) {
  at <- karyoploteR::autotrack(i, length(tracks2), r0 = out2.at$r0, r1 = out2.at$r1, margin = .1)
  kp <- karyoploteR::kpPlotBigWig(kp, data = tracks2[i], ymax = maxValue, #"visible.region"
                                  r0 = at$r0, r1 = at$r1, col = track.Col2[i], border = track.Col2[i])
  computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
        karyoploteR::kpAxis(kp,
                            ymin = 0, ymax = computed.ymax, tick.pos  =  c(0, computed.ymax), tick.len = tickLength,
                            r0 = at$r0, r1 = at$r1, cex = 1)
  # individual bigwig label
  karyoploteR::kpAddLabels(kp, labels  =  names(tracks2)[i], r0 = at$r0, r1 = at$r1, cex = 1, label.margin = .03)
  }
 }
}
