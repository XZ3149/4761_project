#!/bin/bash
cp ../result.bin ./

awk '{print $4}' A549* | paste result.bin - > result_1.bin 
awk '{print $4}' A673* | paste result_1.bin - > result.bin
awk '{print $4}' Caco-2* | paste result.bin - > result_1.bin
awk '{print $4}' GM12878* | paste result_1.bin - > result.bin
awk '{print $4}' H1* | paste result.bin - > result_1.bin
awk '{print $4}' H9* | paste result_1.bin - > result.bin
awk '{print $4}' HCT116* | paste result.bin - > result_1.bin
awk '{print $4}' HepG2* | paste result_1.bin - > result.bin
awk '{print $4}' IMR-90* | paste result.bin - > result_1.bin
awk '{print $4}' K562* | paste result_1.bin - > result.bin
awk '{print $4}' MCF-7* | paste result.bin - > result_1.bin
awk '{print $4}' OCI-LY7* | paste result_1.bin - > result.bin
awk '{print $4}' PC-3* | paste result.bin - > result_1.bin
awk '{print $4}' PC-9* | paste result_1.bin - > result.bin
awk '{print $4}' Panc1* | paste result.bin - > result_1.bin
cat ../cell_lines_header result_1.bin > result_final.bin



