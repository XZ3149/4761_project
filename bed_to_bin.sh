#!/bin/bash
awk '{print $1"\t"$2"\t"$3"\t"1}' $1 > temp.bed

bedtools sort -i temp.bed > temp2.bed
bedtools coverage -a $2 -b temp.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' > ${1}_computed.bin



