#!/bin/bash

# Usage: ./sh find_intron_without_circRNA.sh genome.gff/bed circRNA.gff/bed 
# extra exon and gene
grep 'exon' $1 > exon.gff
grep 'gene' $1 > gene.gff

# extra intron
bedtools subtract -a gene.gff -b exon.gff > intron.gff

# find intron covered circRNA
bedtools interset -a intron.gff -wa > result

# find intron without circRNA
diff intron.gff result | grep '<' > intron_without_circRNA.gff

# remove tmp file
rm exon.gff
rm intron.gff
rm result
