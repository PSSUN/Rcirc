
##' @export
##' @import ggplot2
##' @import dplyr
##' @import RColorBrewer
##' @import circlize
##' @author Haoming Wang & Peisen Sun
##' @param circbed circRNA bed file
##' @param gff genome annotation file in gff format
##' @param genomefasta genome fasta file

showOverview <- function(circbed, gff, genomefasta, ribo, rna) {
  message("Make sure that the chromosomal information of your two annotation files is corresponding, and the non-corresponding chromosomalinformation will cause an error.")
  ribo <- read.delim(ribo, header = FALSE)
  rna <- read.delim(rna, header = FALSE)
  bed <- read.delim(circbed,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(bed) < 3) {
    stop("circbed must contain at least three columns (Chr, Start, End).")
  }
  colnames(bed)[1:3] <- c("Chr", "Start", "End")
  if (ncol(bed) >= 4) {
    colnames(bed)[4] <- "Name"
  } else {
    bed$Name <- paste0("circ_", seq_len(nrow(bed)))
  }
  total_genome <- Biostrings::readDNAStringSet(genomefasta, format = "fasta")
  windowSize <- 50000


  # divide circRNA into four class:
  # intergenic exon-intron exon intron
  genomegff <- read.delim(gff,
    header = FALSE,
    comment.char = "#",
    stringsAsFactors = FALSE
  )
  
  # Extract the first transcript to culculate
  genomegff <- dplyr::filter(genomegff, stringr::str_extract(genomegff$V9, "[.]1;") == ".1;")
  
  gene <- genomegff[genomegff$V3 == "gene", ]
  exon <- genomegff[genomegff$V3 == "exon", ]

  gene <- data.frame(gene$V1, gene$V4, gene$V5, gene$V9)
  id <- sapply(gene$gene.V9, function(x) stringr::str_match(as.character(x), "gene:(.*?);")[, 2])
  gene$gene.V9 <- id

  exon <- data.frame(exon$V1, exon$V4, exon$V5, exon$V9)
  id2 <- sapply(exon$exon.V9, function(x) stringr::str_match(as.character(x), "exon_id=(.*?);")[, 2])
  exon$exon.V9 <- id2

  tmp_chr <- bed$Chr

  gr_gene <- GenomicRanges::GRanges(
    seqnames = gene$gene.V1,
    meta = gene$gene.V9,
    IRanges::IRanges(start = gene$gene.V4, end = gene$gene.V5)
  )
  gr_exon <- GenomicRanges::GRanges(
    seqnames = exon$exon.V1,
    meta = exon$exon.V9,
    IRanges::IRanges(start = exon$exon.V4, end = exon$exon.V5)
  )

  # Find out the complement of gene and exon,
  # and set the complement as intron
  gr_intron <- GenomicRanges::setdiff(x = gr_gene, y = gr_exon)
  # set the meta information as 'intron'
  gr_intron$meta <- "intron"

  gr_circ <- GenomicRanges::GRanges(
    seqnames = tmp_chr,
    meta = bed$Name,
    IRanges::IRanges(
      start = bed$Start,
      end = bed$End
    )
  )

  chromosome <- subset(genomegff, genomegff$V3 == "chromosome")
  chromosome_gr <- GenomicRanges::GRanges(
    seqnames = chromosome$V1,
    meta = chromosome$V3,
    ranges = IRanges::IRanges(
      start = chromosome$V4,
      end = chromosome$V5
    )
  )
  df <- data.frame(
    name = chromosome$V1,
    start = chromosome$V4,
    end = chromosome$V5,
    stringsAsFactors = FALSE
  )
  if (nrow(df) == 0) {
    stop("No chromosome records found in the annotation file.")
  }
  intergenic_gr <- GenomicRanges::setdiff(chromosome_gr, gr_gene)

  intergenic_gr$meta <- "intergenic"
  ######
  total <- c(intergenic_gr, gr_exon, gr_intron)
  total <- sort(total)
  total1 <- as.data.frame(total)

  total <- dplyr::distinct(total1)
  total <- GenomicRanges::GRanges(
    seqnames = total$seqnames,
    meta = total$meta,
    IRanges::IRanges(start = total$start, end = total$end)
  )
  #####
  circ_start <- GenomicRanges::GRanges(
    seqnames = tmp_chr,
    meta = bed$Name,
    IRanges::IRanges(start = bed$Start, end = bed$Start)
  )

  circ_end <- GenomicRanges::GRanges(
    seqnames = tmp_chr,
    meta = bed$Name,
    IRanges::IRanges(start = bed$End, end = bed$End)
  )
  #############
  start <- findOverlaps(circ_start, total, select = "arbitrary")
  end <- findOverlaps(circ_end, total, select = "arbitrary")
  #############
  gr_circ$left <- total[start]$meta
  gr_circ$right <- total[end]$meta
  circ_data <- as.data.frame(gr_circ)
  intron_circ <- subset(
    circ_data,
    circ_data$right == "intron" &
      circ_data$left == "intron"
  )
  intron_circ$type <- "intron_circRNA"
  interg_circ <- subset(
    circ_data,
    circ_data$right == "intergenic" &
      circ_data$left == "intergenic"
  )
  interg_circ$type <- "intergenic_circRNA"
  interg_intron_circ <- subset(
    circ_data,
    (circ_data$right == "intergenic" & circ_data$left == "intron") |
      (circ_data$left == "intergenic" & circ_data$right == "intron")
  )
  interg_intron_circ$type <- "intergenic_intron_circRNA"

  same_exon_circ <- subset(
    circ_data,
    grepl("exon", circ_data$left) &
      grepl("exon", circ_data$right) &
      circ_data$left == circ_data$right
  )

  same_exon_circ$type <- "same_exon_circRNA"
  diff_exon_circ <- subset(
    circ_data,
    grepl("exon", circ_data$left) &
      grepl("exon", circ_data$right) &
      circ_data$left != circ_data$right
  )
  diff_exon_circ$type <- "different_exon_circRNA"

  exon_intron_circ <- subset(
    circ_data,
    (grepl("exon", circ_data$left) & circ_data$right == "intron") |
      (grepl("exon", circ_data$right) & circ_data$left == "intron")
  )
  exon_intron_circ$type <- "exon_intron_circRNA"

  interg_exon_circ <- subset(
    circ_data,
    (grepl("exon", circ_data$left) & circ_data$right == "intergenic") |
      (grepl("exon", circ_data$right) & circ_data$left == "intergenic")
  )
  interg_exon_circ$type <- "intergenic_exon_circRNA"

  final_data <- rbind(
    same_exon_circ,
    diff_exon_circ,
    exon_intron_circ,
    intron_circ,
    interg_circ,
    interg_exon_circ,
    interg_intron_circ
  )

  cla_circ <- final_data
  if (nrow(cla_circ) == 0) {
    stop("No circRNAs available after classification; nothing to plot.")
  }

  bed1 <- data.frame(
    chr = cla_circ$seqnames,
    start = cla_circ$start,
    end = cla_circ$end,
    value = cla_circ$width,
    type = cla_circ$type,
    stringsAsFactors = FALSE
  )

  bed1$type <- factor(
    bed1$type,
    levels = c(
      "same_exon_circRNA",
      "different_exon_circRNA",
      "exon_intron_circRNA",
      "intron_circRNA",
      "intergenic_circRNA",
      "intergenic_exon_circRNA",
      "intergenic_intron_circRNA"
    )
  )
  bed1$type_id <- as.integer(bed1$type)
  bed1 <- bed1[!is.na(bed1$type_id), ]
  if (nrow(bed1) == 0) {
    stop("No circRNAs passed classification; nothing to plot.")
  }
  bed1$value <- log(bed1$value, 2)

  palette_length <- max(nrow(df), 1)
  c1 <- colorRampPalette(brewer.pal(9, "Set3"))(palette_length)

  circos.par(gap.degree = 2, track.margin = c(0, 0))

  ######################################################################

  bedDF <- data.frame()
  genome_fa <- Rsamtools::FaFile(genomefasta)
  Rsamtools::open(genome_fa)
  on.exit(Rsamtools::close(genome_fa))
  for (i in seq_along(total_genome)) {
    seq_name <- names(total_genome)[i]
    if (is.na(seq_name)) {
      seq_name <- paste0("seq", i)
    }
    chr_label <- strsplit(seq_name, " ")[[1]][1]
    seq_width <- Biostrings::width(total_genome[i])
    ran <- IRanges::IRanges(
      start = seq(1, seq_width, by = windowSize),
      width = windowSize
    )
    IRanges::end(ran[length(ran)]) <- seq_width
    totalgenome_bins <- GenomicRanges::GRanges(
      seqnames = rep(chr_label, length(ran)),
      ranges = ran,
      strand = "*"
    )
    bins_seq <- Biostrings::getSeq(
      genome_fa,
      totalgenome_bins
    )
    GC_content <- as.numeric(Biostrings::letterFrequency(bins_seq, "GC")) /
      as.numeric(Biostrings::letterFrequency(bins_seq, "ATCG"))
    tmpDF <- data.frame(
      chr = chr_label,
      start = IRanges::start(ran),
      end = IRanges::end(ran),
      gc_content = GC_content,
      stringsAsFactors = FALSE
    )
    bedDF <- rbind(bedDF, tmpDF)
  }

  f2 <- circlize::colorRamp2(breaks = c(0.2, 0.8), colors = c("white", "red"))
  #####################################################################


  circos.genomicInitialize(df)
  circos.track(
    ylim = c(0, 1),
    bg.col = c1,
    bg.border = NA, track.height = 0.05
  )
  circos.genomicTrackPlotRegion(bed1,
    ylim = c(min(bed1$value), max(bed1$value)),
    numeric.column = c("value", "type_id"),
    panel.fun = function(region, value, ...) {
      type_palette <- c(
        "#FFFF00",
        "#6633FF",
        "#3366CC",
        "#009999",
        "#33CCFF",
        "#66FF99",
        "#FF6633"
      )
      idx <- value[[2]]
      idx[is.na(idx)] <- 1
      col <- type_palette[pmax(pmin(idx, length(type_palette)), 1)]
      circos.genomicPoints(region, value, col = col, cex = 0.2, pch = 16, numeric.column = c("value"))
    }, track.height = 0.2,
#    circos.yaxis(side = 'left',at = seq(min(bed1$value),
#                                        max(bed1$value),
#                                        by=(max(bed1$value)-min(bed1$value))/3
#                                        )
#                 )

  )

  circos.genomicDensity(bed1, track.height = 0.1, bg.border = NA)

  ### GC###
  # Genome annotation dataframe: genomegff
  circos.genomicTrackPlotRegion(bedDF[, c("chr", "start", "end", "gc_content")],
    track.height = 0.1,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value,
        col = f2(value[[1]]),
        border = NA, ...
      )
    }
  )
  ### END###
  circos.genomicDensity(rna, col = "#99DBFF", track.height = 0.15)
  circos.genomicDensity(ribo, col = "#F97C7C", track.height = 0.15)
  title("Overview of total circRNAs")
  circos.clear()
}

