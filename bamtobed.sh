#!/bin/bash
for var in $1/ctcf $1/dnase $1/h3k9me3 $1/rna
do
	samtools sort ${var}.bam -o ${var}.sort.bam
	samtools index ${var}.sort.bam
	samtools view -b ${var}.sort.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 > ${var}.sort.filter.bam
	bedtools coverage -sorted -a 1000.genome.sort.filter.bed -b ${var}.sort.filter.bam -counts > ${var}.bed
	rm ${var}.sort.bam
	rm ${var}.sort.filter.bam
done
