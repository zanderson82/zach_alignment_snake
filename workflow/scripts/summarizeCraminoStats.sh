#!/bin/bash

while getopts "i:" option; do
  case $option in
    i) FILE="$OPTARG";;
  esac
done


GB=$(sed '4q;d' $FILE | cut -f 2)
READS=$(sed '2q;d' $FILE | cut -f 2)
N50=$(sed '7q;d' $FILE | cut -f 2)
MEANRL=$(sed '9q;d' $FILE | cut -f 2)
GB25kb=$(sed '6q;d' $FILE | cut -f 2)
IDENTITY=$(sed '12q;d' $FILE | cut -f 2)
COVALL=$(echo "scale=3; $GB / 3.1" | bc)
COV25=$(echo "scale=3; $GB25kb / 3.1" | bc)

echo "$FILE,$GB,$READS,$N50,$MEANRL,$GB25kb,$IDENTITY,$COVALL,$COV25"
