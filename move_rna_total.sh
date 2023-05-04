#!/bin/bash
for file in `ls .`
do
	if [ -d ./$file ]
	then
		for bin in ./$file/rna*bed
		do
			cp $bin ../combined_files/rna/${file}_rna.bin
		done
	
	fi
done
