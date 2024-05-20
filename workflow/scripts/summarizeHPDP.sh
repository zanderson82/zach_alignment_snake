#!/bin/bash

while getopts "i:" option; do
  case $option in
    i) INPUTFILE="$OPTARG";;
  esac
done

paste <(tail -n+2 $INPUTFILE | sed -r '/^\s*$/d' | cut -f1 | cut -d "-" -f1) \
<( paste <(tail -n+2 $INPUTFILE | sed -r '/^\s*$/d') <(tail -n+2 $INPUTFILE | sed -r '/^\s*$/d' \
| cut -f3 | cut -d ":" -f1) ) | tr '\t' ','

