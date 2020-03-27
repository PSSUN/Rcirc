# ---
# *title*
# circular RNA summary and data preparation for circ_plot.R

# *param*
# bamfile: ribo-seq data mapping to virtual genome
# gff: circular annotation file corresponding virtual genome with junction site position in 10th columns

# *return*
# a circular RNA summary list containing: summary dataframe; trimmed circular RNA annotation; trimmed ribo-seq reads;
# a View object recording reads count at each site

# *details*
# we extract ribosome profiling reads which span circular RNA junction site and count reads on virtual genome aimming
# to facilitate the joining with circular RNA visualization with circ_plot.R

# *example*
#
# *author*
# Wang Haoming https://github.com/AlexWanghaoming
# ---

##' @export circ_summary
##' @import GenomicAlignments
##' @import GenomicFeatures
# bam <- "/Users/sunpeisen/Downloads/total.bam"
# gff <- "/Users/sunpeisen/gene.gff"

circ_summary <- function(bamfile, gff) {
  library(GenomicAlignments)
  library(GenomicFeatures)
  # make circ RNA Grange object
  gff_df <- read.table(gff, sep = "\t")
  colnames(gff_df) <-
    c(
      "chr",
      "source",
      "type",
      "start",
      "end",
      "r1",
      "strand",
      "r2",
      "description",
      "junction_site"
    )
  gr <- makeGRangesFromDataFrame(gff_df, keep.extra.columns = T)

  # make ribo-seq reads Grange object
  ga <- readGAlignments(bamfile)
  read.gr <- granges(ga)
  seqlevels(read.gr) <- "1" 
  seqnames(read.gr) <- sub("_CircularRNA", "", seqnames(read.gr))

  ## filter circRNA whose junction sites are spanned by reads
  # get +-1 bp around junction sites
  junc_site.gr <- GRanges(
    seqnames = 1,
    ranges = IRanges(
      start = gr$junction_site - 1,
      end = gr$junction_site + 1
    ),
    strand = "+"
  )
  keep_junc <-
    subsetByOverlaps(junc_site.gr,
      read.gr,
      ignore.strand = TRUE,
      minoverlap = 3
    )
  sites <- start(keep_junc) + 1
  trimmed.gr <- gr[which(gr$junction_site %in% sites)]


  ## filter ribo-seq reads
  mergeReads.gr <-
    GenomicRanges::reduce(read.gr) # merge ribo-seq reads together
  merged_reads_ovlwith_juncsites <-
    subsetByOverlaps(mergeReads.gr, junc_site.gr, ignore.strand = T)
  trimmed.reads <-
    subsetByOverlaps(read.gr, merged_reads_ovlwith_juncsites)

  ## calculate ribo-seq reads coverage length
  cov.rle <- coverage(trimmed.reads)
  circ_view <- Views(cov.rle, trimmed.gr)[[1]]
  coverage.length <-
    viewApply(circ_view, function(x)
      length(x[which(x != 0)]))

  ## calculate ribo-seq reads counts on junction sites
  readsCount_onCirc <-
    countOverlaps(trimmed.gr, trimmed.reads, ignore.strand = TRUE)

  ## calculate distances of reads coverage region border (3', 5') to junction sites

  dis.mat <- matrix(ncol = 2, nrow = length(sites))
  # head(dis.mat)
  # head(coverage.length)
  for (i in 1:length(sites)) {
    rec <- as.vector(circ_view[[i]])
    view_start <- start(circ_view[i])
    js <- sites[i] - view_start
    cov_region <- which(rec != 0)
    up_dis <- js - cov_region[1]
    down_dis <- cov_region[length(cov_region)] - js
    dis.mat[i, ] <- c(up_dis, down_dis)
  }

  # build circ RNA info dataframe
  name <- paste0("Circ", 1:NROW(trimmed.gr))
  length <- width(trimmed.gr)
  start <- start(trimmed.gr)
  end <- end(trimmed.gr)
  Reads_coverage_length <- coverage.length
  Reads_number <- readsCount_onCirc
  up_distance <- paste0("-", dis.mat[, 1])
  down_distance <- paste0("-", dis.mat[, 2])
  res.df <- data.frame(
    name = name,
    length = length,
    start = start,
    end = end,
    Reads_coverage_length = Reads_coverage_length,
    Reads_number = Reads_number,
    up_distance = up_distance,
    down_distance = down_distance
  )
  info <- list(summary = res.df, trimmed_circRNA = trimmed.gr, trimmed_reads = trimmed.reads, circ_view = circ_view)
  return(info)
}

# res <- circ_summary(bamfile = bam, gff = gff)
