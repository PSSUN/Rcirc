
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
  print("Make sure that the chromosomal information of your two annotation files is corresponding, and the non-corresponding chromosomalinformation will cause an error.")
  genomegff <- gff
  ribo <- read.csv(ribo, sep = "\t", header = FALSE)
  rna <- read.csv(rna, sep = "\t", header = FALSE)
  ribo <- ribo
  rna <- rna
  bed <- read.csv(bed, sep = "	", col.names = c("Chr", "Start", "End"))
  total_genome <- readDNAStringSet(genomefasta, format = "fasta")
  windowSize <- 50000


  # divide circRNA into four class:
  # intergenic exon-intron exon intron
  genomegff <- read.csv(genomegff, sep = "	", header = FALSE)
  
  # Extract the first transcript to culculate
  genomegff <- filter(genomegff,str_extract(genomegff$V9,'[.]1;')=='.1;')
  
  gene <- genomegff[genomegff$V3 == "gene", ]
  exon <- genomegff[genomegff$V3 == "exon", ]
  mRNA <- genomegff[genomegff$V3 == "mRNA", ]

  gene <- data.frame(gene$V1, gene$V4, gene$V5, gene$V9)
  id <- sapply(gene$gene.V9, function(x) str_match(as.character(x), "gene:(.*?);")[, 2])
  gene$gene.V9 <- id

  exon <- data.frame(exon$V1, exon$V4, exon$V5, exon$V9)
  id2 <- sapply(exon$exon.V9, function(x) str_match(as.character(x), "exon_id=(.*?);")[, 2])
  exon$exon.V9 <- id2

  tmp_chr <- sapply(bed$Chr, function(x) gsub("Chr", "", x))

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
  gr_intron <- GenomicRanges::setdiff(x = gr_gene, y = gr_exon)
  # set the meta information as 'intron'
  gr_intron$meta <- "intron"

  gr_circ <- GRanges(
    seqnames = tmp_chr,
    meta = bed$Name,
    IRanges(
      start = bed$Start,
      end = bed$End
    )
  )

  chromosome <- subset(genomegff, genomegff$V3 == "chromosome")
  chromosome_gr <- GRanges(
    seqnames = chromosome$V1,
    meta = chromosome$V3,
    ranges = IRanges(
      start = chromosome$V4,
      end = chromosome$V5
    )
  )
  intergenic_gr <- GenomicRanges::setdiff(chromosome_gr, gr_gene)

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
  total
  #####
  circ_start <- GRanges(
    seqnames = tmp_chr,
    meta = bed$Name,
    IRanges(start = bed$Start, end = bed$Start)
  )

  circ_end <- GRanges(
    seqnames = tmp_chr,
    meta = bed$Name,
    IRanges(start = bed$End, end = bed$End)
  )
  #############
  start <- findOverlaps(circ_start, total, select = "arbitrary")
  end <- findOverlaps(circ_end, total, select = "arbitrary")
  #############
  circ_start
  total

  gr_circ$left <- total[start]$meta
  gr_circ$right <- total[end]$meta
  gr_circ
  #####
  start_info <- subsetByOverlaps(total, circ_start)
  intersect(total, circ_start)
  start_info
  end_info <- subsetByOverlaps(total, circ_end)
  end_info
  circ_data <- as.data.frame(gr_circ)
  circ_data
  intron_circ <- subset(
    circ_data,
    circ_data$right == "intron" &
      circ_data$left == "intron"
  )
  intron_circ$type <- "intron_circRNA"
  intron_circ
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



  cla_circ <- final_data

  gene <- genomegff[genomegff$V3 == "gene", ]
  exon <- genomegff[genomegff$V3 == "exon", ]
  chr <- subset(genomegff, genomegff$V3 == "chromosome")
  name <- chr$V1
  start <- chr$V4
  end <- chr$V5
  c1 <- sample(colors(), length(name))
  df <- data.frame(
    name = name,
    start = start,
    end = end
  )
  tmp_chr <- sapply(bed$Chr, function(x) gsub("Chr", "", x))

  bed1 <- data.frame(
    chr = cla_circ$seqnames,
    start = cla_circ$start,
    end = cla_circ$end,
    value = cla_circ$width,
    type = cla_circ$type
  )

  bed1$type <- as.numeric(bed1$type)
  newpalette <- colorRampPalette(brewer.pal(9, "Set3"))(length(levels(bed1$chr)))
  c1 <- newpalette

  bed1$value <- log(bed1$value, 2)

  circos.par(gap.degree = 2, track.margin = c(0, 0))

  ######################################################################

  bedDF <- data.frame()
  for (i in c(1:length(total_genome))) {
    print(i)
    chrLabels <- unlist(strsplit(names(total_genome[i]), " "))[1]
    ran <- IRanges(
      start = seq(1, width(total_genome[i]), by = windowSize),
      width = windowSize
    )
    end(ran[length(ran)]) <- width(total_genome[i])
    totalgenome_bins <- GRanges(seqnames = rep(chrLabels, length(ran)), ranges = ran, strand = "*")
    bins_seq <- getSeq(
      FaFile(genomefasta),
      totalgenome_bins
    )
    GC_content <- (letterFrequency(bins_seq, "GC") / letterFrequency(bins_seq, "ATCG"))
    tmpDF <- data.frame("chr" = i, "start" = ran@start, "end" = ran@start + ran@width - 1, "gc_content" = GC_content)
    bedDF <- rbind(bedDF, tmpDF)
  }

  f2 <- colorRamp2(breaks = c(0.2, 0.8), colors = c("white", "red"))
  #####################################################################


  circos.genomicInitialize(df)
  circos.track(
    ylim = c(0, 1),
    bg.col = c1,
    bg.border = NA, track.height = 0.05
  )
  circos.genomicTrackPlotRegion(bed1,
    ylim = c(min(bed1$value), max(bed1$value)),
    numeric.column = c("value", "type"),
    panel.fun = function(region, value, ...) {


      col <- ifelse(value[[2]] == 7, "#FF6633",
        ifelse(value[[2]] == 6, "#66FF99",
          ifelse(value[[2]] == 5, "#33CCFF",
            ifelse(value[[2]] == 4, "#009999",
              ifelse(value[[2]] == 3, "#3366CC",
                ifelse(value[[2]] == 2, "#6633FF", "#FFFF00")
              )
            )
          )
        )
      )
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
  circos.genomicTrackPlotRegion(bedDF[, c("chr", "start", "end", "G.C")],
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

