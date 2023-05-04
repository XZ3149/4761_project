#!/bin/bash
for file in `ls .`
do
	if [ -d ./$file ]
	then
		for bed in ./$file/*.bed
		do
			sh bed_to_bin.sh  $bed ./1000.genome.bed
		done
		echo '$file is done'
	fi
done
