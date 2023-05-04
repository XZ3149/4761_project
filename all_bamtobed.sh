#!/bin/bash
for file in `ls .`
do
	if [ -d ./$file ]
	then
		sh bamtobed_ver2.sh $file && touch $file.finished
		
	fi
done
