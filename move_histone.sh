#!/bin/bash
for file in `ls .`
do
	if [ -d ./$file ]
	then
		for bin in ./$file/h3k9*.bed
		do
			cp $bin ../combined_files/h3k9me3/${file}_h3k9me3.bin
		done
	
	fi
done
