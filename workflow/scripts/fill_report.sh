#!/bin/bash

LOGDIR="logs"
READUNTIL=0

while getopts "o:l:b:c:d:rh:" option; do
  case $option in
    o) OUTPUT="$OPTARG" ;;
    l) LIBRARYNAME="$OPTARG" ;;
    b) BAMFILES="$OPTARG";;
    c) CRAMINO="$OPTARG" ;;
    d) LOGDIR="$OPTARG" ;;
    r) READUNTIL=1 ;;
    h) HPDP="$OPTARG" ;;
  esac
done

echo "settings:"
echo "output: $OUTPUT"
echo "library name: $LIBRARYNAME"
echo "bamfiles: $BAMFILES"
echo "cramino: $CRAMINO"
echo "log dir: $LOGDIR"
echo "readuntil: $READUNTIL"
echo "HPDP: $HPDP"
echo ""

firstfile=$(ls -t $LOGDIR/$LIBRARYNAME* | tail -n1)
lastfile=$(ls -t $LOGDIR/$LIBRARYNAME* | head -n1)

startdate=$(date -r $firstfile)
enddate=$(date -r $lastfile)

startSeconds=$(date -d "$startdate" +%s)
endSeconds=$(date -d "$enddate" +%s)
diffSeconds=$(($endSeconds-$startSeconds))

printElapsed=$( date -d @${diffSeconds} +"%H:%M:%S" -u)

bamPrint=$( echo $BAMFILES| tr ',' '\n' | awk '{print "<li>"$0"</li>"}')

# Fill header


sed -i -e "s/LIBRARYNAME/$LIBRARYNAME/g" $OUTPUT
sed -i -e "s/STARTDATE/$startdate/g" $OUTPUT
sed -i -e "s/STOPDATE/$enddate/g" $OUTPUT
sed -i -e "s/TIMEELAPSED/$printElapsed/g" $OUTPUT
sed -i -e "s@<!--BAMFILES-->@$bamPrint@g" $OUTPUT

# Fill table

declare -A SUMMARYSTATS

SUMMARYSTATS["STATGB"]=$(sed '4q;d' $CRAMINO | cut -f 2)
N50=$(sed '7q;d' $CRAMINO | cut -f 2)
SUMMARYSTATS["STATN50"]=$(echo "scale=2; $N50 / 1000" | bc)
SUMMARYSTATS["STATREADS"]=$(sed '2q;d' $CRAMINO | cut -f 2)
SUMMARYSTATS["STATCOV"]=$(sed '5q;d' $CRAMINO | cut -f 2)
SUMMARYSTATS["STATMODAL"]=$(sed '13q;d' $CRAMINO | cut -f 2)

for stat in ${!SUMMARYSTATS[@]}
do
  sed -i -z -e "s/$stat/${SUMMARYSTATS[$stat]}/" $OUTPUT
done

## Other tables
declare -A TABLES

TABLES["HPDPTABLE"]="$LIBRARYNAME/$LIBRARYNAME.hp_dp.snippet.html"
TABLES["VEPTARGETTABLE"]="$LIBRARYNAME/$LIBRARYNAME.vep.target.snippet.html"
TABLES["VEPPATHTABLE"]="$LIBRARYNAME/$LIBRARYNAME.vep.pathogenic.snippet.html"

for table in ${!TABLES[@]}
do
  printf '%s\n' "/$table/r ${TABLES[$table]}" w | ed $OUTPUT
done

#temporary SV table fill

declare -A SVs
declare -A SVCOUNTS
SVs["CUTESVCOUNT"]="$LIBRARYNAME/$LIBRARYNAME.sv_cutesv.phased.vcf"
SVs["SNIFFLESCOUNT"]="$LIBRARYNAME/$LIBRARYNAME.sv_sniffles.phased.vcf"
SVs["SVIMCOUNT"]="$LIBRARYNAME/$LIBRARYNAME.sv_svim.phased.vcf"

for sv in ${!SVs[@]}
do
  if compgen -G ${SVs[$sv]} > /dev/null
  then
    count=$(wc -l "${SVs[$sv]}" | tr -s ' ' | cut -d ' ' -f1)
  else
    count=0
  fi
  sed -i -e "s/$sv/$count/g" $OUTPUT
done


# Find the pictures
declare -A PICS
PICS["LENGTHPLOT"]="$LIBRARYNAME/$LIBRARYNAME.plot_readlengths.png"
PICS["INDELPLOT"]="$LIBRARYNAME/$LIBRARYNAME.plot_indel_quality.png"
PICS["SNVPLOT"]="$LIBRARYNAME/$LIBRARYNAME.plot_snv_quality.png"
PICS["COVERAGEPLOT"]="$LIBRARYNAME/$LIBRARYNAME.plot_depth_coverage.png"
PICS["EFFREADS"]="$LIBRARYNAME/$LIBRARYNAME.read_efficiency.png"
PICS["EFFYIELD"]="$LIBRARYNAME/$LIBRARYNAME.yield_efficiency.png"
PICS["EFFLENGTH"]="$LIBRARYNAME/$LIBRARYNAME.n50_efficiency.png"
PICS["ZYGOPLOT"]="$LIBRARYNAME/$LIBRARYNAME.clair3.notPhased.snpPlot.png"
PICS["PHASEPLOT"]="$LIBRARYNAME/$LIBRARYNAME.clair3.phased.whatshap_plot.png"
PICS["CNVPLOT"]="$LIBRARYNAME/$LIBRARYNAME.called_cnv.detail_plot.png"

# Fill pictures
for pic in ${!PICS[@]}
do
  echo '<img src="data:img/png;base64, '"$( base64 -w0 ${PICS[$pic]})"'">' > tmp.base64
  printf '%s\n' "/$pic/r tmp.base64" w | ed $OUTPUT
done


#haplotagging plots

HAPLOS=( $(ls $LIBRARYNAME/$LIBRARYNAME.hp_dp.detail_plot.png.*.png) )
rm -rf tmp.base64

for pic in ${HAPLOS[@]}
do
  echo '<div class="tile-wrapper tile"><img src="data:img/png;base64, '"$( base64 -w0 $pic)"'"></div>' >> tmp.base64
done

printf '%s\n' "/HAPLOPICS/r tmp.base64" w | ed $OUTPUT

rm -rf tmp.base64

if [ $READUNTIL -eq 1 ]
then
  rucramino=${CRAMINO%*.cramino.stats}.target.cramino.stats
  declare -A RUSUMMARYSTATS

  RUSUMMARYSTATS["RUGB"]=$(sed '4q;d' $rucramino | cut -f 2)
  ruN50=$(sed '7q;d' $rucramino | cut -f 2)
  RUSUMMARYSTATS["RUN50"]=$(echo "scale=2; $ruN50 / 1000" | bc)
  RUSUMMARYSTATS["RUREADS"]=$(sed '2q;d' $rucramino | cut -f 2)

  # get coverage from hpdp file
  rucov=$( awk -F "," 'NR==1{count=0;depth=0}NR>1{if($2 != "FMR" && $2 != "COL1A1"){count+=1;depth+=$4}}END{if(count > 0){print depth/count}else{print 0}}' $HPDP)

  #RUSUMMARYSTATS["RUCOV"]=$(sed '5q;d' $rucramino | cut -f 2)
  RUSUMMARYSTATS["RUCOV"]=$rucov
  coords=( $( tail -n+2 $HPDP | cut -d ',' -f3 ))
  totlength=$(echo ${coords[@]} | tr ' ' '\n' | tr ':' '\t' | tr '-' '\t' | awk 'BEGIN{lengths=0}{lengths+=($3-$2)}END{print lengths}')
  RUSUMMARYSTATS["RUGENCOV"]="$( echo "scale=4;$totlength/3100000000" | bc ) %"
  RUSUMMARYSTATS["RUPERCREADS"]=$( echo "scale=4; ${RUSUMMARYSTATS['RUREADS']} / ${SUMMARYSTATS['STATREADS']}" | bc)
  RUSUMMARYSTATS["RUPERCGB"]=$( echo "scale=4; ${RUSUMMARYSTATS['RUGB']} / ${SUMMARYSTATS['STATGB']}" | bc)

  for stat in ${!RUSUMMARYSTATS[@]}
  do
    sed -i -z -e "s/$stat/${RUSUMMARYSTATS[$stat]}/" $OUTPUT
  done

  declare -A RUPICS
  RUPICS["RUPLOTLENGTH"]="$LIBRARYNAME/$LIBRARYNAME.target.plot_readlengths.png"
  RUPICS["RUPLOTINDEL"]="$LIBRARYNAME/$LIBRARYNAME.target.plot_indel_quality.png"
  RUPICS["RUPLOTSNV"]="$LIBRARYNAME/$LIBRARYNAME.target.plot_snv_quality.png"

  #Fill pictures
  for pic in ${!RUPICS[@]}
  do
    echo '<img src="data:img/png;base64, '"$( base64 -w0 ${RUPICS[$pic]})"'">' > tmp.base64
    printf '%s\n' "/$pic/r tmp.base64" w | ed $OUTPUT
  done
fi
