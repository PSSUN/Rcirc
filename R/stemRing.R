
##' @export
##' @author Sun Peisen
##' @import Biostrings
##' @import IRanges
##' @import dplyr
##' @import IRanges


stemRing <- function(circbed, genomefasta, out, len =1000) {
  ttt <- read.csv(circbed, sep = "\t", header = F)
  c <- data.frame(V1 = ttt$V1, V2 = ttt$V2, V3 = ttt$V3)
  circbed <- c
  circbed$V1 <- gsub(circbed$V1, pattern = "Chr", replacement = "")
  circbed$V1 <- gsub(circbed$V1, pattern = "chr", replacement = "")
  fasta <- readDNAStringSet(genomefasta)
  # len is the length of upstream and down stream
  len <- len
  # n is the index of fasta names
  n <- 1
  name_list <- unique(circbed$V1)
  for (i in names(fasta)) {
    i <- substring(i, 1, 1)
    names(fasta)[n] <- i
    n <- n + 1
    print(i)
  }
  circ <- GRanges(
    seqnames = circbed$V1,
    id = paste(circbed$V2, "-", circbed$V3, sep = ""),
    IRanges(start = circbed$V2, end = circbed$V3)
  )

  circ_up <- GRanges(
    seqnames = circbed$V1,
    id = paste(circbed$V2, "-", circbed$V3, sep = ""),
    IRanges(start = circbed$V2 - len, end = circbed$V2 - 1)
  )
  circ_down <- GRanges(
    seqnames = circbed$V1,
    id = paste(circbed$V2, "-", circbed$V3, sep = ""),
    IRanges(start = circbed$V3 + 1, end = circbed$V3 + len)
  )

  total_list <- c()
  summary_table <- data.frame()
  for (i in names(fasta)) {
    if (i %in% name_list) {
      tmp_up <- circ_up[circ_up@seqnames == i]
      tmp_up_str <- substring(
        fasta[i],
        tmp_up@ranges@start,
        tmp_up@ranges@start + tmp_up@ranges@width
      )
      tmp_down <- circ_down[circ_down@seqnames == i]
      tmp_down_str <- substring(
        fasta[i],
        tmp_down@ranges@start,
        tmp_down@ranges@start + tmp_down@ranges@width
      )
      tmp_al <- pairwiseAlignment(
        type = "local",
        pattern = reverseComplement(DNAStringSet(tmp_up_str)),
        subject = DNAStringSet(tmp_down_str)
      )
      tmp_table <- data.frame(
        chr = tmp_down@seqnames,
        range = tmp_up@elementMetadata,
        up_stream = paste("-", (len - tmp_al@pattern@range@start), sep = ""),
        down_stream = paste("+", tmp_al@subject@range@start, sep = ""),
        width = tmp_al@subject@range@width
      )
      summary_table <- rbind(summary_table, tmp_table)
      # Filter the result by alignment length
    }

    else {
      print(paste(i, "is not in circRNA bed file"))
    }
  }
  filtered_table <- filter(summary_table, width >= 50)

  tmp_start_fun <- function(x) {
    unlist(strsplit(x, "-")[[1]][1])
  }
  tmp_end_fun <- function(x) {
    unlist(strsplit(x, "-")[[1]][2])
  }
  tmp_forward_fun <- function(x) {
    unlist(strsplit(x, "-")[[1]][2])
  }
  tmp_back_fun <- function(x) {
    unlist(strsplit(x, "+ ")[[1]][1])
  }

  tmp_s <- as.numeric(unlist(lapply(filtered_table$id, tmp_start_fun)))
  tmp_e <- as.numeric(unlist(lapply(filtered_table$id, tmp_end_fun)))
  tmp_f <- as.numeric(unlist(lapply(as.character(filtered_table$up_stream), tmp_forward_fun)))
  tmp_b <- as.numeric(unlist(lapply(as.character(filtered_table$down_stream), tmp_back_fun)))

  filtered_table$up_stream_start <- as.numeric(tmp_s) - tmp_f
  filtered_table$up_stream_end <- as.numeric(tmp_s) - tmp_f + as.numeric(filtered_table$width)
  filtered_table$down_stream_start <- as.numeric(tmp_e) + tmp_b
  filtered_table$down_stream_end <- as.numeric(tmp_e) + tmp_b + as.numeric(filtered_table$width)

  write.table(filtered_table, file = out, sep = "\t", row.names = F, quote = F)
}


