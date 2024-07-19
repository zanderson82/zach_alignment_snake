#!/bin/bash

OLD=1

while getopts "i:n" option; do
  case $option in
    i) INPUTFILE="$OPTARG";;
    n) OLD=0;;
  esac
done

if [ "$OLD" -eq 1]
then
  paste <(tail -n+2 $INPUTFILE | sed -r '/^\s*$/d' | cut -f1 | cut -d "-" -f1) \
  <( paste <(tail -n+2 $INPUTFILE | sed -r '/^\s*$/d') <(tail -n+2 $INPUTFILE | sed -r '/^\s*$/d' \
  | cut -f3 | cut -d ":" -f1) ) | tr '\t' ','
else
  echo tail -n+2 $INPUTFILE
done
