#!/bin/bash
for file in `ls .`
do
	if [ -d ./$file ]
	then
		for bin in ./$file/dna*.bed
		do
			cp $bin ../combined_files/DNase/${file}_DNase.bin
		done
	
	fi
done
