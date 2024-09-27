#!/bin/bash

# needs bedtools

set +o pipefail

ORIGINALDIR=$(pwd)

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

cd "$parent_path"

SERVER="McClintock"
TEMPLOC=/tmp

CLAIRQUALITY=20.0
MINIMUMDEPTH=4

while getopts "i:o:q:d:T:" option; do
  case $option in
    i) FILE="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    q) CLAIRQUALITY="$OPTARG" ;;
    d) MINIMUMDEPTH="$OPTARG" ;;
    T) TEMPLOC="$OPTARG" ;;
  esac
done

bash annotate_omim.sh -i $FILE -o $TEMPLOC/temp.csv

declare -A annotationFiles

annotationFiles["Mappability"]=/n/dat/annotationData/catted.ok.Mappability.merged.bed
annotationFiles["GRCExclusions"]=/n/dat/annotationData/catted.ok.Exclusions.merged.bed
annotationFiles["Repeats"]=/n/dat/annotationData/catted.ok.AllRepeats.merged.bed
annotationFiles["UCSC_Unusual"]=/n/dat/annotationData/catted.ok.UCSC_unusual.sorted.merged.bed
annotationFiles["SegmentalDuplications"]=/n/dat/annotationData/catted.ok.UCSC_segdup.sorted.merged.bed
annotationFiles["Homopolymers"]=/n/dat/annotationData/catted.ok.homopolymer.sorted.merged.bed

keys=( "Mappability" "GRCExclusions" "Repeats" "UCSC_Unusual" "SegmentalDuplications" "Homopolymers" )
for key in ${keys[@]}
do
    bash annotate_vep_with_bed.sh -i $TEMPLOC/temp.csv -o $TEMPLOC/temp2.csv -a ${annotationFiles[$key]} -n $key
    mv $TEMPLOC/temp2.csv $TEMPLOC/temp.csv
done

bash prioritizeVep.sh -i $TEMPLOC/temp.csv -o ${OUTPUT%*.csv} -c -q $CLAIRQUALITY -d $MINIMUMDEPTH
mv $TEMPLOC/temp.csv $OUTPUT

echo "running makeRMD.sh -i ${OUTPUT%*.csv}.prioritized.csv -o ${OUTPUT%*.csv}.prioritized.report.html -l ${OUTPUT%*.csv}.log"
bash makeRMD.sh -i ${OUTPUT%*.csv}.prioritized.csv -o ${OUTPUT%*.csv}.prioritized.report.html -l ${OUTPUT%*.csv}.log

cd $ORIGINALDIR

exit 0