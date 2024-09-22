#!/bin/bash
#bash getUMI.sh outputDir/PreprocessDir/clean_CCS.out.txt outputDir/PreprocessDir/clean_CCS.out.cut.txt outputDir/PreprocessDir/clean_CCS.umi.txt

infile=$1
outfile=$2
umifile=$3

awk -vOFS='\t' -vUMI=$umifile '{split($2,box,":");$2 = box[1]":"box[2]":"box[3];print $0; print box[1]":"box[2]":"box[3]"\t"box[4]"\t"box[5]"\t"box[6] > UMI }' $infile > $outfile


