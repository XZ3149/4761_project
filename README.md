# 4761_project


This is the GitHub repo for 4761 final project by


- Step 1. data downloading from ENCODE. We create dictionary for each cell line. Then, we download data through bash script using wget in Linux server. The download.sh is an example script that we use to download cell line A549 bam files. We run the script as 

```
sh download.sh
```


- Step 2. data binirization

Bed and BigWig files were downloaded from ENCODE. Firstly, we use bash script to transfrom all BigWig format into Bed format using bigwig_to_bed.sh. 
 
```
cd woring_directory/ 
sh bigwig_to_bed.sh
```















```
test
```
