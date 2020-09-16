##' @export
##' @import Biostrings
##' @import ggplot2
##' @param x

showCodon <- function(x) {
circ <- readDNAStringSet(x)

codon <- trinucleotideFrequency(circ)
codon <- data.frame(codon)
codon <- apply(codon, MARGIN = 2, sum)

codon <- data.frame(codon, type = names(codon))
codon <- codon[sort(codon$codon,index.return=TRUE)$ix,]
x <- codon$type
y <- codon$codon
ggplot(codon) +
  geom_bar(stat = "identity", aes(x = x, y = y), fill = "lightblue") +
  coord_flip() + scale_x_discrete(limits = rev(x)) +
  labs(title = "Distribution of 3-nt pattern", x = "Type", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank())+
  scale_y_continuous(expand = c(0,0))

}
