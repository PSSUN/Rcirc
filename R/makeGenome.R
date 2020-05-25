
##' @import Biostrings
##' @author Peisen Sun

##' @export makeGenome
makeGenome <- function(fastaIn, fastaOut, gffOut, Nnum = 200) {
  seqs <- readDNAStringSet(filepath = fastaIn, format = "fasta")
  Nnum <- Nnum
  N <- ""
  for (i in 1:Nnum) {
    N <- paste(N, "N", sep = "")
  }

  increase_length <- 0
  start_list <- c()
  end_list <- c()
  total_seq <- ""
  n <- 1

  for (i in 1:length(seqs)) {
    seq_str <- toString(seqs[i])
    len <- nchar(seq_str)
    if (len <= 500) {
      changed_seq <- paste(seq_str, seq_str, N, sep = "")
      total_seq <- paste(total_seq, changed_seq, sep = "")
      increase_length <- increase_length + len * 2 + Nnum
    }

    else {
      changed_seq <- paste(reverse(reverse(seq_str)[1:500]),
        seq_str[1:500],
        N,
        sep = ""
      )
      increase_length <- increase_length + 500 * 2 + Nnum
    }


    if (n == 1) {
      start_position <- 1
      end_position <- increase_length - Nnum
      next_start <- increase_length + 1
    }
    else {
      start_position <- next_start
      end_position <- increase_length - Nnum
      next_start <- increase_length + 1
    }
    n <- n + 1
    start_list <- c(start_list, start_position)
    end_list <- c(end_list, end_position)
    print(n)
  }

  gff <- data.frame(
    chr = 1,
    second = "araport11",
    gene = "gene",
    start = start_list,
    end = end_list,
    third = ".",
    forth = "+",
    fifth = ".",
    sixth = "."
  )

  write.table(
    x = gff,
    file = gffOut,
    sep = "\t",
    col.names = F,
    row.names = F,
    quote = F
  )

  final_seq <- DNAStringSet(x = total_seq)
  names(final_seq)[1] <- "newgenome"
  out <- fastaOut
  writeXStringSet(x = final_seq, filepath = out, format = "fasta")
}

# fastaIn <- "/home/sun/arabidopsis_Fan_seq.fa"
# gffOut <- "/home/sun/Rcirc/test.gff"
# fastaOut <- "/home/sun/Rcirc/test.fa"

# makeGenome(fastaIn, fastaOut, gffOut)
