#!/bin/bash
for file in `ls .`
do
	if [ -d ./$file ]
	then
		bigWigToWig ./$file/*minus*.bigWig ./$file/rnaminus.wig
		bigWigToWig ./$file/*plus*.bigWig ./$file/rnaplus.wig
		wig2bed --zero-indexed < ./$file/rnaminus.wig > ./$file/rnaminus.bed
		wig2bed --zero-indexed < ./$file/rnaplus.wig > ./$file/rnaplus.bed
	fi
done
