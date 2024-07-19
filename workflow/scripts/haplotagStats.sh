#!/bin/bash

RUSTSCRIPT="workflow/scripts/./haplotagParse"
HEADER="File,Region,Coords,Depth,HP1,HP2,PercentAssigned,Phaseblocks"
FMR="chrX:147900000-148000000"
COL1A1="chr17:50150000-50250000"
LONGFORM=0;


while getopts "i:r:b:l" option; do
    case $option in 
        i) INPUT=$OPTARG ;;
        r) TARGET=$OPTARG ;;
        b) BEDFILE=$OPTARG ;;
        l) LONGFORM=1 ;;
    esac
done

if [ $LONGFORM -eq 1 ]
then
    RUSTSCRIPT="workflow/scripts/./haplotagParseLong"
    HEADER="File,Region,Chromosome,Position,Depth,HP1,HP2,PercentAssigned"
fi

if [ -z ${BEDFILE+x} ]
then
# make temp for FMR
    coord="${FMR##*:}"
    chr="${FMR%:*}"
    seq ${coord%-*} 1000 ${coord##*-} | awk -v chrom=$chr '{print chrom,$1,$1+1}' | tr ' ' '\t' > temp.positions

    region=$FMR

    samtools mpileup -r $region --positions temp.positions --ff SUPPLEMENTARY,UNMAP,SECONDARY,QCFAIL,DUP -q 1 -Q 1 --output-extra "PS,HP" $INPUT > temp.pileup.tsv

    FMR_RESULTS=$( cut -f8 temp.pileup.tsv | $RUSTSCRIPT )

    rm temp.positions
    rm temp.pileup.tsv

    ###make temp for COL1A1

    coord="${COL1A1##*:}"
    chr="${COL1A1%:*}"
    seq ${coord%-*} 1000 ${coord##*-} | awk -v chr=$chr '{print chr,$1,$1+1}' | tr ' ' '\t' > temp.positions

    region=$COL1A1

    samtools mpileup -r $region --positions temp.positions --ff SUPPLEMENTARY,UNMAP,SECONDARY,QCFAIL,DUP -q 1 -Q 1 --output-extra "PS,HP" $INPUT > temp.pileup.tsv

    COL1A1_RESULTS=$( cut -f8 temp.pileup.tsv | $RUSTSCRIPT )

    rm temp.positions
    rm temp.pileup.tsv

    ###make temp for $TARGET

    coord="${TARGET##*:}"
    chr="${TARGET%:*}"
    seq ${coord%-*} 1000 ${coord##*-} | awk -v chr=$chr '{print chr,$1,$1+1}' | tr ' ' '\t' > temp.positions

    region=$TARGET

    samtools mpileup -r $region --positions temp.positions --ff SUPPLEMENTARY,UNMAP,SECONDARY,QCFAIL,DUP -q 1 -Q 1 --output-extra "PS,HP" $INPUT > temp.pileup.tsv

    TARGET_RESULTS=$( cut -f8 temp.pileup.tsv | $RUSTSCRIPT )

    rm temp.positions
    rm temp.pileup.tsv


    ###output results

    echo -e "$HEADER"
    echo -e "$INPUT\tFMR\t$FMR\t"${FMR_RESULTS[@]} | tr ' ' '\t'
    echo -e "$INPUT\tCOL1A1\t$COL1A1\t"${COL1A1_RESULTS[@]} | tr ' ' '\t'
    echo -e "$INPUT\tTARGET\t$TARGET\t"${TARGET_RESULTS[@]} | tr ' ' '\t'

else
    printresults=( "$HEADER" )
    while read line
    do
        coordparts=( $line )
        seq ${coordparts[2]} 1000 ${coordparts[3]} | awk -v chr=${coordparts[1]} '{print chr,$1,$1+1}' | tr ' ' '\t' > temp.positions
        region=${coordparts[1]}":"${coordparts[2]}"-"${coordparts[3]}
        samtools mpileup -r $region --positions temp.positions --ff SUPPLEMENTARY,UNMAP,SECONDARY,QCFAIL,DUP -q 1 -Q 1 --output-extra "PS,HP" $INPUT > temp.pileup.tsv
        if [ "$LONGFORM" -eq 1 ]
        then
            TARGET_RESULTS=( $( cut -f8 temp.pileup.tsv | $RUSTSCRIPT) )
            numTargets="${#TARGET_RESULTS[@]}"
            n=( $(seq 1 "$numTargets") )
            printresults+=( $(paste -d ',' <(paste -d ',' <( printf "$INPUT,${coordparts[0]}\n"%.0s "${n[@]}" ) <( cut -f1,2 temp.pileup.tsv | tr '\t' ',' ) ) <( echo ${TARGET_RESULTS[@]} | tr ' ' '\n') ))
        else
            TARGET_RESULTS=$( cut -f8 temp.pileup.tsv | $RUSTSCRIPT)
            PRINT_TARGETS=$(echo ${TARGET_RESULTS[@]} | tr ' ' ',')
            set -o noglob
            PS=( $( cut -f7 temp.pileup.tsv | tr ',' ' ' | tr ' ' '\n' | sort | uniq) )
            numPS="${#PS[@]}"
            if [ "${PS[0]}" == '*' ]
            then
                numPS=$(echo "$numPS - 1" | bc)
            fi
            set +o noglob
            result=$(echo "$INPUT,${coordparts[0]},$region,$PRINT_TARGETS,$numPS")
            printresults+=( $result ) 
        fi
        #rm temp.positions
        #rm temp.pileup.tsv
           
    done <<< $(cat $BEDFILE)

    echo -e "${printresults[@]}" | tr ' ' '\n'
fi

