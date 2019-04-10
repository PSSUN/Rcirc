#!/bin/bash

# Usage: ./sh find_intron_without_circRNA.sh genome.gff/bed circRNA.gff/bed 
# extra exon and gene
grep 'exon' $1 > exon.gff
grep 'gene' $1 > gene.gff

bedtools subtract -a gene.gff -b exon.gff > intron.gff
rm exon.gff
bedtools interset -a intron.gff -wa > result
diff intron.gff result | grep '<' > intron_without_circRNA.gff
rm intron.gff
rm result

