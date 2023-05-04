# 4761_project


This is the GitHub repo for 4761 final project by

Dependencies:
python 3.9
anaconda
jupyternotebook
samtools 
bedtools 
wget 
matplotlib 
numpy
seaborn 
scipy
pandas
numpy
sklearn




- Step 1. data downloading from ENCODE for cell lines described in cell_lines.txt. We download data from https://www.encodeproject.org/matrix/?type=Experiment&control_type!=*&status=released&perturbed=false. We create dictionary for each cell line. Then, we download data through bash script using wget in Linux server. The download.sh is an example script that we use to download cell line A549 bam files. We run the script as 

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


- Step 3. Data binirization

Bed and BigWig files were downloaded from ENCODE as described earlier. Firstly, we use bash script to transfrom all BigWig format into bed format using bigwig_to_bed.sh. 
 
```
cd woring_directory/ 
sh bigwig_to_bed.sh
```

An example working directory will look like this:

![My Image](Images/Sample_directory.png)



Then, we process all bed files to binary bin with another bash script all_bash_to_bin.sh. The all_bash_to_bin.sh have dependencis on bed_to.bin.sh and 1000.genome.bed.

```
cd woring_directory/ 
sh all_bed_to_bin.sh
```



- Step 4. Data aggregation

We use bash script move_ctcf.sh, move_DNase.sh, move_histone.sh, move_datarna_total.sh to moved respected seqencing data together and use combined.sh to merge them together as a single file. The script require the input of cell_lines_header file


```
cd woring_directory/data_type/
sh combined.sh
```

We will obtain data in format like this: 


![My Image](Images/Combined_file.png)



- Step 5. Bam files process 

After trying binirized data, we decide to process bam to bin. 




```
test
```
