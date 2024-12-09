#!/bin/bash

RAWDF=tmp.rawdf.tsv
BASECALLPATH="/data/prealign_qc/dorado_summary/qual_only"
EXPLICIT=0

while getopts "c:l:o:r:b:e" option; do
  case $option in
    c) CRAMINO="$OPTARG" ;;
    l) LIBRARY="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    r) RAWDF="$OPTARG" ;;
    b) BASECALLPATH="$OPTARG" ;;
    b) EXPLICIT=1 ;;
  esac
done

echo -e "Stage\tStat\tValue" > $OUTPUT
echo "$OUTPUT"
if [[ $EXPLICIT -eq 0 ]]
then
    libraries=( $(ls $LIBRARY) )
else
    bams=( $(echo "$LIBRARY" |  tr '?' ' ') )
    libraries=()
    for bam in ${bams[@]}
    do
        libraries+=( $(echo "${bam%%.*}" ) )
    done
fi

# Get sequencing run reports
seqPath1=/waldo/original_data
seqPath2=/data/incoming-seqdata/mv.to.waldo

reportpaths=()
for library in ${libraries[@]}
do
    libraryname=${library##*/}
    librarytrunk=${libraryname%*.dorado*}
    if compgen -G "$seqPath1/$librarytrunk*/*/sequencing_summary*txt" > /dev/null
    then
        reportpaths+=( $(ls $seqPath1/$librarytrunk*/*/sequencing_summary*txt) )
    elif compgen -G "$seqPath2/$librarytrunk*/*/sequencing_summary*txt" > /dev/null
    then
        reportpaths+=( $(ls $seqPath2/$librarytrunk*/*/sequencing_summary*txt) )
    else
        echo "No sequencing stats were found for $librarytrunk"
    fi
done

rawoutput=$RAWDF
for report in ${reportpaths[@]}
do
    cut -f11 $report | tail -n+2 >> $rawoutput
done

seqRead=$( wc -l $rawoutput | tr ' ' '\t' | cut -f1)
totalDur=$( awk 'NR==1{tot=$1}NR>1{if($1>0.5){tot+=$1}}END{printf "%.2f", tot}' $rawoutput)
seqYield=$( echo "$totalDur * 362" | bc)
halfDur=$(echo "$totalDur / 2" | bc)
seqN50b=$( sort -nrk1,1 $rawoutput | awk -v tl="$halfDur" 'NR==1{sum=$1}NR>1{if($1>0.5){sum=sum+$1};if(sum>tl){print $1; exit 0;}}')
seqN50=$(echo "scale=2; $seqN50b * 362" | bc)

echo -e "Sequencing\tYield\t$seqYield" >> $OUTPUT
echo -e "Sequencing\tReads\t$seqRead" >> $OUTPUT
echo -e "Sequencing\tN50\t$seqN50" >> $OUTPUT

# Get basecalling reports
bcPath=$BASECALLPATH
bcYieldGB=0
bcRead=0
bcN50=0
for library in ${libraries[@]}
do
    libraryname=${library##*/}
    librarytrunk=${libraryname%*report*}
    if [ -z "$librarytrunk" ]; then
    librarytrunk=${libraryname%*dorado*}
    fi
    if compgen -G $bcPath/$librarytrunk*html > /dev/null
    then
        report=$(ls $bcPath/$librarytrunk*html)
        yieldparts=( $(grep -A1 "Yield (GB)" $report) )
        yieldstr=${yieldparts[2]}
        newYield=$(echo "$bcYieldGB + ${yieldstr//[^0-9.]/}" | bc)
        bcYieldGB=$newYield
        readparts=( $(grep -A1 "Reads" $report) )
        readstr=${readparts[1]}
        newReads=$(echo "$bcRead + ${readstr//[^0-9]/}" | bc)
        bcRead=$newReads
        n50parts=( $(grep -A1 "full N50" $report) )
        n50str=${n50parts[3]}
        newN50=$(echo "$bcN50 + ${n50str//[^0-9.]/}" | bc)
        bcN50=$newN50
    else
        echo "No basecalling stats were found for $librarytrunk"
    fi
done

bcYield=$( echo "$bcYieldGB * 1000000000" | bc)

averageBCN50n=$(echo "scale=2;$bcN50 / ${#libraries[@]}" | bc)
averageBCN50=$(echo "$averageBCN50n * 1000" | bc)

echo -e "Basecalling\tYield\t$bcYield" >> $OUTPUT
echo -e "Basecalling\tReads\t$bcRead" >> $OUTPUT
echo -e "Basecalling\tN50\t$averageBCN50" >> $OUTPUT

# Get alignment reports

cramstats=$CRAMINO
alignYieldgb=$(sed '4q;d' $cramstats | cut -f 2)
alignYield=$( echo "$alignYieldgb * 1000000000" | bc)
alignReads=$(sed '2q;d' $cramstats | cut -f 2)
alignN50=$(sed '7q;d' $cramstats | cut -f 2)

echo -e "Alignment\tYield\t$alignYield" >> $OUTPUT
echo -e "Alignment\tReads\t$alignReads" >> $OUTPUT
echo -e "Alignment\tN50\t$alignN50" >> $OUTPUT