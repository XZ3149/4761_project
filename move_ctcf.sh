#!/bin/bash
for file in `ls .`
do
	if [ -d ./$file ]
	then
		for bin in ./$file/ct*.bed
		do
			cp $bin ../combined_files/ctcf/${file}_ctcf.bin
		done
	
	fi
done
