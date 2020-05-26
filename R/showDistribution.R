##' @export
##' @import GenomicRanges
##' @import IRanges
##' @import grid
##' @import gtable
##' @import ggplot2
##' @import stringr
##' @param gbed The genome annotation file: gtf or gff3
##' @param cbed The circRNA bed file


showDistribution <- function(gbed, cbed) {
  bed <- read.csv(cbed, header = FALSE)
  # divide circRNA into four class:
  # intergenic exon-intron exon intron
  genomegff <- read.table(gbed, sep = '\t')
  # genomegff <- filter(genomegff,str_extract(gff$V9,'[.]1;')=='.1;')
  gene <- genomegff[genomegff$V3 == "gene", ]
  exon <- genomegff[genomegff$V3 == "exon", ]
  mRNA <- genomegff[genomegff$V3 == "mRNA", ]

  gene <- data.frame(gene$V1, gene$V4, gene$V5, gene$V9)
  id <- sapply(gene$gene.V9, function(x) str_match(as.character(x), "gene:(.*?);")[, 2])
  gene$gene.V9 <- id

  exon <- data.frame(exon$V1, exon$V4, exon$V5, exon$V9)
  id2 <- sapply(exon$exon.V9, function(x) str_match(as.character(x), "exon_id=(.*?);")[, 2])
  exon$exon.V9 <- id2

  tmp_chr <- sapply(bed$V1, function(x) gsub("Chr", "", x))

  gr_gene <- GRanges(
    seqnames = gene$gene.V1,
    meta = gene$gene.V9,
    IRanges(start = gene$gene.V4, end = gene$gene.V5)
  )
  gr_exon <- GRanges(
    seqnames = exon$exon.V1,
    meta = exon$exon.V9,
    IRanges(start = exon$exon.V4, end = exon$exon.V5)
  )

  # Find out the complement of gene and exon,
  # and set the complement as intron
  gr_intron <- setdiff(x = gr_gene, y = gr_exon)
  # set the meta information as 'intron'
  gr_intron$meta <- "intron"

  gr_circ <- GRanges(
    seqnames = tmp_chr,
    #                     meta = bed$Name,
    IRanges(
      start = bed$V2,
      end = bed$V3
    )
  )

  chromosome <- subset(genomegff, genomegff$V3 == "chromosome")
  chromosome
  chromosome_gr <- GRanges(
    seqnames = chromosome$V1,
    meta = chromosome$V3,
    ranges = IRanges(
      start = chromosome$V4,
      end = chromosome$V5
    )
  )

  intergenic_gr <- setdiff(chromosome_gr, gr_gene)
  intergenic_gr$meta <- "intergenic"
  ######
  total <- c(intergenic_gr, gr_exon, gr_intron)
  total <- sort(total)
  total1 <- as.data.frame(total)
  total <- distinct(total1)
  total <- GRanges(
    seqnames = total$seqnames,
    meta = total$meta,
    IRanges(start = total$start, end = total$end)
  )
  circ_start <- GRanges(
    seqnames = tmp_chr,
    #                        meta = bed$Name,
    IRanges(start = bed$V2, end = bed$V2)
  )

  circ_end <- GRanges(
    seqnames = tmp_chr,
    #                      meta = bed$Name,
    IRanges(start = bed$V3, end = bed$V3)
  )

  start <- findOverlaps(circ_start, total, select = "arbitrary")
  end <- findOverlaps(circ_end, total, select = "arbitrary")


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
    circ_data$right == "intergenic" &
      circ_data$left == "intron" |
      circ_data$left == "intergenic" &
        circ_data$right == "intron"
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
    grepl("exon", circ_data$left) &
      circ_data$right == "intron" |
      grepl("exon", circ_data$right) &
        circ_data$left == "intron"
  )
  exon_intron_circ$type <- "exon_intron_circRNA"

  interg_exon_circ <- subset(
    circ_data,
    grepl("exon", circ_data$left) &
      circ_data$right == "intergenic" |
      grepl("exon", circ_data$right) &
        circ_data$left == "intergenic"
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

  ggplot(final_data, aes(x = seqnames, fill = type)) + geom_bar(stat = "count") +
    labs(x = "chromosome", title = "circRNA types distribution on chromosome") +
    scale_fill_brewer(palette = "Set2") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
}


