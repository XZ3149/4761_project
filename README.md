# 4761_project


This is the GitHub repo for 4761 final project by


- Step 1. data downloading from ENCODE. We create dictionary for each cell line. Then, we download data through bash script using wget in Linux server. The download.sh is an example script that we use to download cell line A549 bam files. We run the script as 

```
sh download.sh
```

- Step 2. Reference hg38 bed genertaion
We downloaded hg38.fa as our reference genome and create a reference bin file using following commands:

```
# bin the chromosome based on reference 
# use -w to indicate the bin width

samtools dict hg38.fa > hg38.fa.dict
cat hg38.fa.dict|grep -v '^@HD'|sed 's/:/\t/g'|cut -f 3,5 >genome.chr.ln
bedtools makewindows -g genome.chr.ln -w 1000 > 1000.genome.bed
```


- Step 2. data binirization

Bed and BigWig files were downloaded from ENCODE. Firstly, we use bash script to transfrom all BigWig format into Bed format using bigwig_to_bed.sh. 
 
```
cd woring_directory/ 
sh bigwig_to_bed.sh
```

Then, we process all bed files to binary bin with another bash script all_bash_to_bin.sh. The all_bash_to_bin.sh have a dependency on bed_to.bin.sh.

```
cd woring_directory/ 
sh all_bed_to_bin.sh
```













```
test
```
