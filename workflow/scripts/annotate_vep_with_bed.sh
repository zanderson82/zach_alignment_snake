#!/bin/bash

# needs bedtools

set +o pipefail

SERVER="McClintock"
annotation=hg38.mapping.filter.sorted.mergephenos.bed

while getopts "i:o:a:n:" option; do
  case $option in
    i) FILE="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    a) annotation="$OPTARG" ;;
    n) annotationname="$OPTARG" ;;
  esac
done

if [ -z {$annotationname+x} ]
then
  annotationname=$annotation
fi

tail -n+2 $FILE | cut -d ',' -f5,6 | awk -F ',' '{print $1, $2, $2+1, NR}' | tr ' ' '\t' > temp.bed

bedtools intersect -a temp.bed -b $annotation -loj | cut -f1,2,3,4,8 > temp.intersect.bed

echo  "$annotationname" > temp.intersect.2.bed

awk 'NR==1{i=$4;pheno=$5;chr=$1;start=$2;stop=$3}\
NR>1{if($4==i){pheno=pheno","$5}\
else{\
print chr,start,stop,pheno;\
chr=$1;start=$2;stop=$3;i=$4;pheno=$5}\
}END{\
print chr,start,stop,pheno}' temp.intersect.bed >> temp.intersect.2.bed

paste -d ',' <(cat $FILE) <(cat temp.intersect.2.bed | cut -d ' ' -f4) | tr -s ',' > $OUTPUT


rm -rf temp.intersect.bed
rm -rf temp.intersect.2.bed
rm -rf temp.bed