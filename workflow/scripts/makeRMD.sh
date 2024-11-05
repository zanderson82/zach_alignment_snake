#!/bin/bash

set +o pipefail

while getopts "i:o:l:" option; do
  case $option in
    i) FILE="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    l) LOG="$OPTARG" ;;
  esac
done

expression=$( echo 'rmarkdown::render("workflow/resources/one_table_vep.Rmd",params=list(vepfile="'$FILE'", logfile="'$LOG'"), output_file="'$OUTPUT'")' )

Rscript -e "$expression"; echo $?

exit 0