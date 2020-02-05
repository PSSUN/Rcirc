##' @export
##' @import ggplot2
##' @import Biostrings
##' @param fasta fasta file of circRNAs
##' @param max max types to show in barplot
##' @param title title of barplot


showJunction <- function(data, max = 20, title = 0) {
  circ <- readDNAStringSet(data)
  tmp <- c()
  for (i in c(1:length(circ))) {
    head <- subseq(circ[i], 1, 2)
    tail <- subseq(
      circ[i],
      as.numeric(circ[i]@ranges@width) - 1,
      as.numeric(circ[i]@ranges@width)
    )
    s <- paste(head, "/", tail)
    tmp <- c(tmp, s)
  }
  tmp <- table(tmp)
  tmp <- as.data.frame(tmp)
  tmp <- tmp[with(tmp, order(-Freq)), ]
  tmp2 <- tmp[1:max, ]
  if (title == 0) {
    tit <- paste0("Top ", max, " splice-signal type")
  }
  else {
    tit <- title
  }
  ggplot(tmp2, aes(x = tmp, y = Freq)) +
    geom_bar(stat = "identity", fill = "lightblue", colour = "black") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.background = element_blank()
    ) +
    xlab("Type") + ylab("Number") +
    ggtitle(tit) +
    scale_x_discrete(limits = tmp2$tmp)
}


