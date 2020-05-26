

###
##' @export
##' @author sun peisen
##' @import stringr
##' @import dplyr
##' @import GenomicRanges
##' @import IRanges
##' @import circlize
##' @import RColorBrewer
classByType <- function(circbed, gff, file = "./result.csv") {

  #    print("Make sure that the chromosomal information of your two annotation files is corresponding, and the non-corresponding chromosomalinformation will cause an error.")

  bed <- read.csv(circbed, sep = "	", col.names = c("Chr", "Start", "End"))
  # divide circRNA into four class:
  # intergenic exon-intron exon intron

  genomegff <- read.csv(gff, sep = "	", header = FALSE)
  
  genomegff <- filter(gff,str_extract(genomegff$V9,'[.]1;')=='.1;')
  
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
  #    tmp_chr <- bed$Chr
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
    meta = bed$Start + bed$End,
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
  # many package has function 'setdiff'
  # here we use GenomicRanges::setdiff
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

  start <- findOverlaps(circ_start, total, select = "arbitrary")
  end <- findOverlaps(circ_end, total, select = "arbitrary")

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
    circ_data$right == "intron" & circ_data$left == "intron"
  )
  intron_circ$type <- "intron_circRNA"
  intron_circ
  interg_circ <- subset(
    circ_data,
    circ_data$right == "intergenic" & circ_data$left == "intergenic"
  )
  interg_circ$type <- "intergenic_circRNA"
  #
  interg_intron_circ <- subset(
    circ_data,
    circ_data$right == "intergenic" & circ_data$left == "intron" |
      circ_data$left == "intergenic" & circ_data$right == "intron"
  )
  interg_intron_circ$type <- "intergenic_intron_circRNA"

  same_exon_circ <- subset(
    circ_data,
    grepl("exon", circ_data$left) & grepl("exon", circ_data$right) &
      circ_data$left == circ_data$right
  )
  same_exon_circ$type <- "same_exon_circRNA"

  diff_exon_circ <- subset(
    circ_data,
    grepl("exon", circ_data$left) & grepl("exon", circ_data$right) &
      circ_data$left != circ_data$right
  )
  diff_exon_circ$type <- "different_exon_circRNA"

  #-------May be used one day-------
  exon_circ <- subset(
    circ_data,
    grepl("exon", circ_data$left) & grepl("exon", circ_data$right)
  )
  diff_exon_circ$type <- "different_exon_circRNA"
  #---------------------------------

  exon_intron_circ <- subset(
    circ_data,
    grepl("exon", circ_data$left) & circ_data$right == "intron" |
      grepl("exon", circ_data$right) & circ_data$left == "intron"
  )
  exon_intron_circ$type <- "exon_intron_circRNA"

  interg_exon_circ <- subset(
    circ_data,
    grepl("exon", circ_data$left) & circ_data$right == "intergenic" |
      grepl("exon", circ_data$right) & circ_data$left == "intergenic"
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


  ##### CLASSIFY######
  cla_circ <- final_data
  ######## PLOT#######

  genomegff <- read.csv(gff, sep = "	", header = FALSE)
  gene <- genomegff[genomegff$V3 == "gene", ]
  exon <- genomegff[genomegff$V3 == "exon", ]
  chr <- subset(genomegff, genomegff$V3 == "chromosome")
  name <- chr$V1
  start <- chr$V4
  end <- chr$V5

  newpalette <- colorRampPalette(brewer.pal(9, "Set3"))(length(name))
  c1 <- newpalette

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
  circos.clear()
  par(mai = c(0.1, 0.1, 0.1, 0.9))
  circos.par(gap.degree = 0.5, track.margin = c(0, 0))
  circos.genomicInitialize(df)
  print("Circlize size:")
  print("-----------------")
  print(df)
  print("-----------------")
  print("Start drawing...")
  circos.track(
    ylim = c(0, 1),
    bg.col = c1,
    bg.border = NA, track.height = 0.05
  )
  #
  circos.genomicDensity(bed1[bed1$type == "same_exon_circRNA", ],
    track.height = 0.08, bg.border = NA, col = c("#5BC0EB")
  )
  circos.genomicDensity(bed1[bed1$type == "different_exon_circRNA", ],
    track.height = 0.08, bg.border = NA, col = c("#FDE74C")
  )
  circos.genomicDensity(bed1[bed1$type == "exon_intron_circRNA", ],
    track.height = 0.08, bg.border = NA, col = c("#9BC53D")
  )
  circos.genomicDensity(bed1[bed1$type == "intergenic_circRNA", ],
    track.height = 0.08, bg.border = NA, col = c("#E55934")
  )
  circos.genomicDensity(bed1[bed1$type == "intergenic_exon_circRNA", ],
    track.height = 0.08, bg.border = NA, col = c("#FA7921")
  )
  circos.genomicDensity(bed1[bed1$type == "intergenic_intron_circRNA", ],
    track.height = 0.08, bg.border = NA, col = c("#3CC453")
  )
  circos.genomicDensity(bed1[bed1$type == "intron_circRNA", ],
    track.height = 0.08, bg.border = NA, col = c("#197EEA")
  )
  legend(0.8, -0.4,
    legend = c(
      "same_exon_circRNA",
      "different_exon_circRNA",
      "exon_intron_circRNA",
      "intergenic_circRNA",
      "intergenic_exon_circRNA",
      "intergenic_intron_circRNA",
      "intron_circRNA"
    ),
    fill = c(
      "#5BC0EB",
      "#FDE74C",
      "#9BC53D",
      "#E55934",
      "#FA7921",
      "#3CC453",
      "#197EEA"
    ),
    border = FALSE, cex = 0.7, bty = "n", ncol = 1, adj = 0, xpd = T
  )
  circos.clear()
  print("Save result data to:")
  print(file)
  write.csv(bed1, file = file, row.names = FALSE)
  print("Finished! Thank you for using Rcirc.")
}

