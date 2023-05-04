#!/bin/bash
for file in `ls .`
do
	if [ -d ./$file ]
	then
		sh bamtobed.sh $file && touch $file.finished
		
	fi
done
