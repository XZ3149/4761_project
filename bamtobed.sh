#!/bin/bash
cd /home/project/samples/try_2
for var in ctcf dnase h3k9me3 rna
do
	samtools sort ${var}.bam -o ${var}.sort.bam
	samtools index ${var}.sort.bam
	samtools view -b ${var}.sort.bam chr{1..22} > ${var}.sort.filter.bam
	bedtools coverage -sorted -a ../1000.genome.sort.filter.bed -b ${var}.sort.filter.bam -counts > ${var}.bed
done
