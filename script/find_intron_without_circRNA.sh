#!/bin/bash

# Usage: ./sh find_intron_without_circRNA.sh genome.gff/bed circRNA.gff/bed 
if [ ! $1 ];then
	echo 'please input genome bed file'
	exit
else
	pass
fi
if [ ! $2 ];then
        echo 'please input circRNA bed file'
	exit
else
        pass
fi

# extra exon and gene
grep 'exon' $1 > exon.gff
grep 'gene' $1 > gene.gff

bedtools subtract -a gene.gff -b exon.gff > intron.gff
rm exon.gff
bedtools intersect -a intron.gff -wa > result
diff intron.gff result | grep '<' > intron_without_circRNA.gff
rm intron.gff
rm result
