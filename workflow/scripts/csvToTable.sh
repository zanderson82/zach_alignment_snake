#!/bin/bash

while getopts "i:o:f:" option; do
  case $option in
    i) FILE="$OPTARG";;
    o) OUTPUT="$OPTARG" ;;
    f) FIELDS="$OPTARG" ;;
  esac
done

tempout=${FILE%*.csv}.temp.csv

cut -d ',' -f"$FIELDS" $FILE > $tempout

echo '<table class="autotable">' > $OUTPUT
line=0

while read INPUT
do
    if [ $line -lt 1 ]
    then
        echo "<thead><tr><th>${INPUT//,/</th><th>}</th></tr></thead><tbody>" >> $OUTPUT
        line=$(echo "$line+1" | bc )
    else
        echo "<tr><td>${INPUT//,/</td><td>}</td></tr>" >> $OUTPUT
    fi
done < $tempout
echo "</tbody></table>" >> $OUTPUT

rm -rf $tempout