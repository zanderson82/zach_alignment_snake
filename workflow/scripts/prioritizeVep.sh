#!/bin/bash

# the g6* model of clair should use a quality threshold around 10. dorado trained models should use a cutoff of 20.

set +o pipefail

CLEAN=0
QUALTHRESH=20
MINIMUMDEPTH=4

while getopts "i:o:cq:d:" option; do
  case $option in
    i) FILE="$OPTARG" ;;
    o) PREFIX="$OPTARG" ;;
    c) CLEAN=1 ;;
    q) QUALTHRESH="$OPTARG" ;;
    d) MINIMUMDEPTH="$OPTARG" ;;
    l) LOGFILE="$OPTARG" ;;
  esac
done

if [ -z ${LOGFILE+x} ]
then
  LOGFILE=$PREFIX.log
fi

header=$(head -n1 $FILE)

countVariants() {
  local infile=$1
  lines=$(wc -l $infile | tr -s ' ' | cut -d ' ' -f1)
  variants=$(echo "$lines-1" | bc)
  echo $variants
}

easy() {
  local infile=$1
  outfile=${infile%*.csv}.easy.csv
  awk -F ',' 'NR==1{print $0}{if($95=="." && $96=="." && $97=="." && $98=="." && $99=="." && $100=="."){print $0}}' $infile > $outfile
  echo "-$(countVariants $outfile) with no mapping issues, homopolymers, repeats, or segdups" | tee -a $LOGFILE
}

keepOMIM() {
  local infile=$1
  outfile=${infile%*.csv}.phenotyped.csv
  awk -F ',' 'NR==1{print $0}NR>1{if($94!="."){print $0}}' $infile > $outfile
  echo "--$(countVariants $outfile) with OMIM record" | tee -a $LOGFILE
  #echo $(wc -l $outfile) " with OMIM record"
}

noOMIM() {
  local infile=$1
  outfile=${infile%*.csv}.nopheno.csv
  awk -F ',' 'NR==1{print $0}NR>1{if($94=="."){print $0}}' $infile > $outfile
  echo "--$(countVariants $outfile) without OMIM record" | tee -a $LOGFILE
}

noCLINVAR() {
  local infile=$1
  outfile=${infile%*.csv}.noclinvar.csv
  awk -F ',' 'NR==1{print $0}{if($91=="."){print $0}}' $infile > $outfile
  echo "---$(countVariants $outfile) with known pathogenic variant" | tee -a $LOGFILE
}

knownVariant() {
  local infile=$1
  outfile=${infile%*.csv}.knownpathogenic.csv
  awk -F ',' 'NR==1{print $0}{if($91=="Pathogenic" || $91=="Pathogenic/Likely_pathogenic"){print $0}}' $infile > $outfile
  echo "---$(countVariants $outfile) with known pathogenic variant" | tee -a $LOGFILE
}

knownVUS() {
  local infile=$1
  outfile=${infile%*.csv}.vus.csv
  awk -F ',' 'NR==1{print $0}{if($91=="Uncertain significance" || $91=="Conflicting_classifications_of_pathogenicity" || $91=="Likely_pathogenic"){print $0}}' $infile > $outfile
  echo "---$(countVariants $outfile) with VUS" | tee -a $LOGFILE
}

dominant() {
  local infile=$1
  outfile=${infile%*.csv}.dominant.csv
  echo $header > $outfile
  awk -F ',' 'NR>1{if($94 ~ "dominant"){print $0}}' $infile | sort -t ',' -k77,77nr >> $outfile
  echo "----$(countVariants $outfile) with dominant inheritance" | tee -a $LOGFILE
}

recessive() {
  local infile=$1
  outfile=${infile%*.csv}.recessive.csv
  echo $header > $outfile
  awk -F ',' 'NR>1{if($94 !~ "dominant"){print $0}}' $infile | sort -t ',' -k17 >> $outfile
  echo "----$(countVariants $outfile) with recessive inheritance" | tee -a $LOGFILE
}

homozygous() {
  local infile=$1
  outfile=${infile%*.csv}.homozygous.csv
  echo $header > $outfile
  awk -F ',' 'NR>1{if($10 == "1/1"){print $0}}' $infile | sort -t ',' -k77,77nr >> $outfile
  echo "----$(countVariants $outfile) homozygous" | tee -a $LOGFILE
}

heterozygous() {
  local infile=$1
  outfile=${infile%*.csv}.heterozygous.csv
  echo $header > $outfile
  awk -F ',' 'NR>1{if($10 == "0/1" || $10=="0|1" || $10=="1|0"){print $0}}' $infile | sort -t ',' -k77,77nr >> $outfile
  echo "----$(countVariants $outfile) heterozygous" | tee -a $LOGFILE
}


hetcompound() {
  local infile=$1
  outfile=${infile%*.csv}.heterozygous_compound.csv
  echo $header > $outfile
  genes=( $(awk -F ',' 'NR==2{gene=$17}NR>2{if(gene==$17){print $17};gene=$17}' $infile | sort | uniq ) )
  if [[ "${#genes[@]}" -gt 0 ]]
  then
    for gene in ${genes[@]}
    do
      grep $gene $infile >> $outfile
    done
  fi
  echo "----$(countVariants $outfile) with heterozygous compound inheritance" | tee -a $LOGFILE
}

# pHaplo >= 0.86
haplosensitive() {
  local infile=$1
  outfile=${infile%*.csv}.haplosensitive.csv
  echo $header > $outfile
  awk -F ',' 'NR>1{if($73 >= 0.86){print $0}}' $infile | sort -t ',' -k77,77nr >> $outfile
  echo "-----$(countVariants $outfile) with haplosensitivity score >= 0.86" | tee -a $LOGFILE
}

inheritance() {
  local infile=$1
  dominant $infile
  recessive $infile
  homozygous ${infile%*.csv}.recessive.csv
  hetcompound ${infile%*.csv}.recessive.csv
}

strongPredPathogenic() {
  local infile=$1
  outfile=${infile%*.csv}.predicted_pathogenic.strong.csv
  awk -F ',' 'NR==1{print $0}NR>1{if($91!="Pathogenic" && $42 ~ "deleterious" && $42 !~ "low" && $43 ~ "damaging" && $43 ~ "probably" && $75=="likely_pathogenic"){print $0}}' $infile > $outfile
  echo "--$(countVariants $outfile) with strongly predicted pathogenicity" | tee -a $LOGFILE
}

midPredPathogenic() {
  local infile=$1
  outfile=${infile%*.csv}.predicted_pathogenic.mid.csv
  awk -F ',' 'NR==1{print $0}NR>1{if(($91=="." || $91=="not_provided") && (($42 ~ "deleterious" && $42 !~ "low")|| $43 ~ "probably_damaging" || $75=="likely_pathogenic")){print $0}}' $infile > $outfile
  echo "--$(countVariants $outfile) with some predicted pathogenicity" | tee -a $LOGFILE
}

highQual() {
  local infile=$1
  outfile=${infile%*.csv}.high_qual.csv
  echo $header > $outfile
  awk -F ',' -v quality=$QUALTHRESH 'NR>1{if($9>quality && $11 > 10 && ($10 == "1/1" || $10=="1|0" || $10=="0|1" || $10=="0/1")){print$0}}' $infile | sort -t ',' -k77,77nr >> $outfile
  echo "--$(countVariants $outfile) with QUAL > $QUALTHRESH, not multiallelic" | tee -a $LOGFILE
}

highCADD() {
  local infile=$1
  outfile=${infile%*.csv}.high_cadd.csv
  echo $header > $outfile
  awk -F ',' 'NR>1{if($77>29){print$0}}' $infile | sort -t ',' -k77,77nr >> $outfile
  echo "--$(countVariants $outfile) with CADD (phred) >= 30, not multiallelic" | tee -a $LOGFILE
}

snv() {
  local infile=$1
  outfile=${infile%*.csv}.snv.csv
  echo $header > $outfile
  awk -F ',' 'NR>1{if($2=="SNV"){print$0}}' $infile | sort -t ',' -k77,77nr >> $outfile
  echo "$(countVariants $outfile) SNVs" | tee -a $LOGFILE
}

indel(){
  local infile=$1
  outfile=${infile%*.csv}.indel.csv
  echo $header > $outfile
  awk -F ',' 'NR>1{if($2=="INDEL"){print$0}}' $infile | sort -t ',' -k77,77nr >> $outfile
  echo "$(countVariants $outfile) INDELs" | tee -a $LOGFILE
}

anyFlag(){
  local infile=$1
  outfile=${infile%*.csv}.flagged.csv
  echo $header > $outfile
  awk -F ',' 'NR>1{if($95 != "." || $96!="." || $97!="." || $98!="." || $99!="." || $100!="."){print$0}}' $infile | sort -t ',' -k77,77nr >> $outfile
  echo "-$(countVariants $outfile) with any flagged issue" | tee -a $LOGFILE
}

allCoding(){
  local SUBFIX=$1
  keepOMIM $SUBFIX.csv
  knownVariant $SUBFIX.phenotyped.csv
  knownVUS $SUBFIX.phenotyped.csv
  inheritance $SUBFIX.phenotyped.vus.csv
  strongPredPathogenic $SUBFIX.phenotyped.csv
  midPredPathogenic $SUBFIX.phenotyped.csv
  inheritance $SUBFIX.phenotyped.predicted_pathogenic.mid.csv
  noCLINVAR $SUBFIX.phenotyped.csv
  highQual $SUBFIX.phenotyped.noclinvar.csv
  highCADD $SUBFIX.phenotyped.noclinvar.high_qual.csv
  highImpact $SUBFIX.phenotyped.noclinvar.high_qual.high_cadd.csv
  noOMIM $SUBFIX.csv
  knownVariant $SUBFIX.nopheno.csv
  knownVUS $SUBFIX.nopheno.csv
  strongPredPathogenic $SUBFIX.nopheno.csv
  midPredPathogenic $SUBFIX.nopheno.csv
  highQual $SUBFIX.nopheno.predicted_pathogenic.mid.csv
  homozygous $SUBFIX.nopheno.predicted_pathogenic.mid.high_qual.csv
  heterozygous $SUBFIX.nopheno.predicted_pathogenic.mid.high_qual.csv
  hetcompound $SUBFIX.nopheno.predicted_pathogenic.mid.high_qual.heterozygous.csv
  haplosensitive $SUBFIX.nopheno.predicted_pathogenic.mid.high_qual.heterozygous.csv
}



## splice AI score exists and at least one DS is >= 0.5
probableSplice(){
  local infile=$1
  outfile=${infile%*.csv}.spliceai_gt05.csv
  echo $header > $outfile
  awk -F ',' 'NR>1{if($83 >= 0.5 || $84 >=0.5 || $85 >=0.5 || $86 >=0.5){print $0}}' $infile | sort -t ',' -k77,77nr >> $outfile
  echo "--$(countVariants $outfile) with probable splice variant" | tee -a $LOGFILE
}

## no spliceAI score exists
noSpliceAI(){
  local infile=$1
  outfile=${infile%*.csv}.nospliceai.csv
  echo $header > $outfile
  awk -F ',' 'NR>1{if($83=="."){print $0}}' $infile | sort -t ',' -k77,77nr >> $outfile
  echo "--$(countVariants $outfile) with no spliceAI scores" | tee -a $LOGFILE
}

## impact is HIGH
highImpact(){
  local infile=$1
  outfile=${infile%*.csv}.high_impact.csv
  echo $header > $outfile
  awk -F ',' 'NR>1{if($16=="HIGH"){print $0}}' $infile | sort -t ',' -k77,77nr >> $outfile
  echo "---$(countVariants $outfile) with high impact" | tee -a $LOGFILE
}

## CADD >= 10
cadd10(){
  local infile=$1
  outfile=${infile%*.csv}.cadd_gt10.csv
  echo $header > $outfile
  awk -F ',' 'NR>1{if($77 >= 10){print $0}}' $infile | sort -t ',' -k77,77nr >> $outfile
  echo "---$(countVariants $outfile) with CADD greater than 10" | tee -a $LOGFILE
}

## CADD >= 20
cadd20(){
  local infile=$1
  outfile=${infile%*.csv}.cadd_gt20.csv
  echo $header > $outfile
  awk -F ',' 'NR>1{if($77 >= 20){print $0}}' $infile | sort -t ',' -k77,77nr >> $outfile
  echo "---$(countVariants $outfile) with CADD greater than 20" | tee -a $LOGFILE
}

echo "$(countVariants $FILE) total starting variants" | tee -a $LOGFILE

# exclude intronic variants, upstream gene variants, downstream gene variants, DP > 1

awk -F ',' -v depth=$MINIMUMDEPTH '{if((($15 !~ "splice" && $15 !~ "intron" && $15 !~ "upstream" && $15 !~ "synonymous" \
&& $15 !~ "UTR" && $15 !~ "downstream" && $15 !~ "intergenic" && $15 !~ "non_coding") || $15 ~ "stop" || $15 ~ "frameshift") && $3 >= depth){print $0}}' $FILE > $PREFIX.coding.csv

echo "$(countVariants $PREFIX.coding.csv) coding variants" | tee -a $LOGFILE

## no mappability issues, no grcexclusions, no repeats, no unusual UCSC, no seg dups, no homopolymers
easy $PREFIX.coding.csv
SUBFIX=$PREFIX.coding.easy

allCoding $SUBFIX

# any mappability flag allowed
anyFlag $PREFIX.coding.csv

# high quality any flag
highQual $PREFIX.coding.flagged.csv

SUBFIX=$PREFIX.coding.flagged.high_qual
allCoding $SUBFIX


##### cleanup

rm -rf $PREFIX.coding.csv
rm -rf $PREFIX.coding.easy.csv
rm -rf $PREFIX.coding.easy.nopheno.csv
rm -rf $PREFIX.coding.easy.nopheno.predicted_pathogenic.mid.csv
rm -rf $PREFIX.coding.easy.nopheno.predicted_pathogenic.mid.high_qual.csv
rm -rf $PREFIX.coding.easy.phenotyped.csv
rm -rf $PREFIX.coding.easy.phenotyped.predicted_pathogenic.mid.csv
rm -rf $PREFIX.coding.easy.phenotyped.vus.csv
rm -rf $PREFIX.coding.flagged.csv
rm -rf $PREFIX.coding.flagged.high_qual.csv
rm -rf $PREFIX.coding.flagged.high_qual.nopheno.csv
rm -rf $PREFIX.coding.flagged.high_qual.nopheno.predicted_pathogenic.mid.csv
rm -rf $PREFIX.coding.flagged.high_qual.nopheno.predicted_pathogenic.mid.high_qual.csv
rm -rf $PREFIX.coding.flagged.high_qual.phenotyped.csv
rm -rf $PREFIX.coding.flagged.high_qual.phenotyped.predicted_pathogenic.mid.csv
rm -rf $PREFIX.coding.flagged.high_qual.phenotyped.vus.csv

########################################## noncoding: exclude stops, frame shifts, missense, lost starts (includes splice variants), DP  > 1
old=$FILE
new=$PREFIX.noncoding.csv
awk -F ',' -v depth=$MINIMUMDEPTH '{if((($15 !~ "frameshift" && $15 !~ "inframe" && $15 !~ "missense" && $15 !~ "stop" && $15 !~ "start") || $15 ~ "non-coding") && $3 >= depth){print $0}}' $old > $new
echo "$(countVariants $new) non coding variants (excluding stops, frameshifts, missense, lost starts)" | tee -a $LOGFILE

####### no mappability issues, no grcexclusions, no repeats, no unusual UCSC, no seg dups, no homopolymers

easy $PREFIX.noncoding.csv

#spliceAI likely

probableSplice $PREFIX.noncoding.easy.csv
## spliceAI + OMIM
keepOMIM $PREFIX.noncoding.easy.spliceai_gt05.csv
## no spliceAI, high impact
noSpliceAI $PREFIX.noncoding.easy.csv

highImpact $PREFIX.noncoding.easy.nospliceai.csv
keepOMIM $PREFIX.noncoding.easy.nospliceai.high_impact.csv

# high CADD
cadd20 $PREFIX.noncoding.easy.nospliceai.csv
keepOMIM $PREFIX.noncoding.easy.nospliceai.cadd_gt20.csv

######## high quality any flag
anyFlag $PREFIX.noncoding.csv
highQual $PREFIX.noncoding.flagged.csv

# probable splice variants
probableSplice $PREFIX.noncoding.flagged.high_qual.csv
# with phenotype
keepOMIM $PREFIX.noncoding.flagged.high_qual.spliceai_gt05.csv

## no spliceAI
noSpliceAI $PREFIX.noncoding.flagged.high_qual.csv
#
highImpact $PREFIX.noncoding.flagged.high_qual.nospliceai.csv

keepOMIM $PREFIX.noncoding.flagged.high_qual.nospliceai.high_impact.csv
# high CADD

cadd20 $PREFIX.noncoding.flagged.high_qual.nospliceai.csv
keepOMIM $PREFIX.noncoding.flagged.high_qual.nospliceai.cadd_gt20.csv
inheritance $PREFIX.noncoding.flagged.high_qual.nospliceai.cadd_gt20.phenotyped.csv

##### cleanup

rm -rf $PREFIX.noncoding.csv
rm -rf $PREFIX.noncoding.easy.csv
rm -rf $PREFIX.noncoding.easy.nospliceai.csv
rm -rf $PREFIX.noncoding.flagged.csv
rm -rf $PREFIX.noncoding.flagged.high_qual.csv
rm -rf $PREFIX.noncoding.flagged.high_qual.nospliceai.cadd_gt20.phenotyped.csv
rm -rf $PREFIX.noncoding.flagged.high_qual.nospliceai.csv

#### concatenate

addOutput(){
  local inFile=$1
  local text=$2  
  local priority=$3

  awk -v priority=$priority -v text="REPLACETEXT" 'NR>1{print $0","text","priority}' $inFile >> $finalfile
  sed -i -e "s/REPLACETEXT/$text/g" $finalfile
  echo -n "."
}

finalfile=$PREFIX.prioritized.csv
echo "$header,priority_reason,priority" > $finalfile

echo -n "writing outputs to final file"

addOutput $PREFIX.coding.easy.phenotyped.knownpathogenic.csv "known pathogenic variant with associated phenotype in a coding non-error prone region" 1
addOutput $PREFIX.noncoding.easy.spliceai_gt05.phenotyped.csv "high probability splice variant with an associated phenotype in a non-error prone region" 2
addOutput $PREFIX.coding.easy.phenotyped.noclinvar.high_qual.high_cadd.high_impact.csv "High impact high cadd high quality variant in a non error prone coding region with an OMIM phenotype and no ClinVar entry" 2
addOutput $PREFIX.coding.easy.nopheno.knownpathogenic.csv "known pathogenic variant in a coding non-error prone region" 3
addOutput $PREFIX.coding.easy.phenotyped.vus.dominant.csv "documented VUS in a coding non error prone region with an associated phenotype and dominant inheritance pattern" 3
addOutput $PREFIX.coding.easy.phenotyped.vus.recessive.heterozygous_compound.csv "documented VUS in a coding non error prone region with an associated phenotype and compound heterozygous inheritance pattern" 3
addOutput $PREFIX.coding.easy.phenotyped.vus.recessive.homozygous.csv "documented VUS in a coding non error prone region with an associated phenotype and homozygous recessive inheritance pattern" 3
addOutput $PREFIX.coding.flagged.high_qual.phenotyped.knownpathogenic.csv "known pathogenic variant of high quality with associated phenotype in a coding error prone region" 3
addOutput $PREFIX.noncoding.easy.nospliceai.cadd_gt20.phenotyped.csv "high CADD score non-coding variant with associated phenotype in a non-error prone region" 3
addOutput $PREFIX.noncoding.easy.nospliceai.high_impact.phenotyped.csv "high impact non-coding variant with associated phenotype in a non-error prone region" 3
addOutput $PREFIX.noncoding.flagged.high_qual.spliceai_gt05.phenotyped.csv "high probability splice variant with an associated phenotype in an error prone region" 2
addOutput $PREFIX.coding.easy.nopheno.vus.csv "documented VUS in a coding non error prone region with no associated phenotype" 4
addOutput $PREFIX.coding.easy.phenotyped.predicted_pathogenic.strong.csv "variant in a coding non error prone region with an associated phenotype and consistently predicted pathogenicity" 4
addOutput $PREFIX.coding.flagged.high_qual.nopheno.knownpathogenic.csv "known pathogenic variant of high quality in a coding error prone region with no associated phenotype" 4
addOutput $PREFIX.coding.flagged.high_qual.phenotyped.vus.dominant.csv "documented VUS of high quality in a coding error prone region with an associated phenotype and dominant inheritance" 4
addOutput $PREFIX.coding.flagged.high_qual.phenotyped.vus.recessive.heterozygous_compound.csv "documented VUS of high quality in a coding error prone region with an associated phenotype and heterozygous compound inheritance" 4
addOutput $PREFIX.coding.flagged.high_qual.phenotyped.vus.recessive.homozygous.csv "documented VUS of high quality in a coding error prone region with an associated phenotype and homozygous recessive inheritance" 4
addOutput $PREFIX.noncoding.flagged.high_qual.nospliceai.cadd_gt20.phenotyped.dominant.csv "high quality variant in a non-coding error prone region with a high CADD score associated phenotype and dominant inheritance" 4
addOutput $PREFIX.noncoding.flagged.high_qual.nospliceai.cadd_gt20.phenotyped.recessive.heterozygous_compound.csv "high quality variant in a non-coding error prone region with a high CADD score associated phenotype and heterozygous compound inheritance" 4
addOutput $PREFIX.noncoding.flagged.high_qual.nospliceai.cadd_gt20.phenotyped.recessive.homozygous.csv "high quality variant in a non-coding error prone region with a high CADD score associated phenotype and homozygous recessive inheritance" 4
addOutput $PREFIX.noncoding.flagged.high_qual.nospliceai.high_impact.phenotyped.csv "high quality variant in a non-coding error prone region with a high impact score and associated phenotype" 4
addOutput $PREFIX.coding.easy.nopheno.predicted_pathogenic.strong.csv "variant in a coding non error prone region with consistently predicted pathogenicity and no phenotype association" 5
addOutput $PREFIX.coding.easy.phenotyped.predicted_pathogenic.mid.dominant.csv "variant in a coding non error prone region with an associated phenotype and some prediction of pathogenicity and dominant inheritance" 5
addOutput $PREFIX.coding.easy.phenotyped.predicted_pathogenic.mid.recessive.heterozygous_compound.csv "variant in a coding non error prone region with an associated phenotype and some prediction of pathogenicity and heterozygous compound inheritance" 5
addOutput $PREFIX.coding.easy.phenotyped.predicted_pathogenic.mid.recessive.homozygous.csv "variant in a coding non error prone region with an associated phenotype and some prediction of pathogenicity and homozygous recessive inheritance" 5
addOutput $PREFIX.coding.flagged.high_qual.nopheno.vus.csv "document VUS of high quality in a coding error prone region with no associated phenotype" 5
addOutput $PREFIX.coding.flagged.high_qual.phenotyped.predicted_pathogenic.strong.csv "variant of high quality in a coding error prone region with an associated phenotype and consistently predicted pathogenicity" 5
addOutput $PREFIX.noncoding.easy.nospliceai.high_impact.csv "high impact non-coding variant in a non-error prone region with no associated phenotype" 5
addOutput $PREFIX.noncoding.easy.spliceai_gt05.csv "high probability splice variant in a non-error prone region with no associated phenotype" 6
addOutput $PREFIX.noncoding.flagged.high_qual.nospliceai.high_impact.csv "high impact and high quality non-coding variant in an error prone region with no associated phenotype" 5
addOutput $PREFIX.noncoding.flagged.high_qual.spliceai_gt05.csv "high quality high probability splice variant in an error prone region with no associated phenotype" 6
addOutput $PREFIX.coding.easy.nopheno.predicted_pathogenic.mid.high_qual.heterozygous.haplosensitive.csv "variant in a coding non error prone region with no associated phenotype and some evidence of pathogenicity and a heterozygous inheritance pattern with high probable haploinsufficiency" 7
addOutput $PREFIX.coding.easy.nopheno.predicted_pathogenic.mid.high_qual.heterozygous.heterozygous_compound.csv "variant in a coding non error prone region with no associated phenotype and some evidence of pathogenicity and a heterozygous compound inheritance pattern" 7
addOutput $PREFIX.coding.easy.nopheno.predicted_pathogenic.mid.high_qual.homozygous.csv "variant in a coding non error prone region with no associated phenotype and some evidence of pathogenicity and a homozygous inheritance pattern" 7
addOutput $PREFIX.coding.easy.phenotyped.vus.recessive.csv "documented VUS in a coding non error prone region with an associated phenotype a heterozygous recessive inheritance pattern " 7
addOutput $PREFIX.coding.flagged.high_qual.nopheno.predicted_pathogenic.strong.csv "variant of high quality in a coding error prone region with no associated phenotype and consistently predicted pathogenicity" 7
addOutput $PREFIX.coding.flagged.high_qual.phenotyped.predicted_pathogenic.mid.dominant.csv "variant of high quality in a codding error prone region with an associated phenotype and some evidence of pathogenicity and a dominant inheritance pattern" 7
addOutput $PREFIX.coding.flagged.high_qual.phenotyped.predicted_pathogenic.mid.recessive.heterozygous_compound.csv "variant of high quality in a codding error prone region with an associated phenotype and some evidence of pathogenicity and a compound heterozygous inheritance pattern" 7
addOutput $PREFIX.coding.flagged.high_qual.phenotyped.predicted_pathogenic.mid.recessive.homozygous.csv "variant of high quality in a codding error prone region with an associated phenotype and some evidence of pathogenicity and a compound homozygous recessive inheritance pattern" 7
addOutput $PREFIX.noncoding.easy.nospliceai.cadd_gt20.csv "non coding variant in a non error prone region with a high CADD score " 7
addOutput $PREFIX.noncoding.flagged.high_qual.nospliceai.cadd_gt20.phenotyped.recessive.csv "high quality variant in a non-coding error prone region with a high CADD score associated phenotype and heterozygous recessive inheritance" 7
addOutput $PREFIX.coding.easy.nopheno.predicted_pathogenic.mid.high_qual.heterozygous.csv "variant in a coding non error prone region with no associated phenotype and some evidence of pathogenicity and a heterozygous recessive inheritance pattern " 8
addOutput $PREFIX.coding.easy.phenotyped.predicted_pathogenic.mid.recessive.csv "variant in a coding non error prone region with an associated phenotype and some prediction of pathogenicity and heterozygous recessive inheritance" 8
addOutput $PREFIX.coding.flagged.high_qual.nopheno.predicted_pathogenic.mid.high_qual.heterozygous.haplosensitive.csv "high quality variant in a coding error prone region with no associated phenotype and some evidence of pathogenicity and a heterozygous inheritance pattern with high probable haploinsufficiency" 8
addOutput $PREFIX.coding.flagged.high_qual.nopheno.predicted_pathogenic.mid.high_qual.heterozygous.heterozygous_compound.csv "high quality variant in a coding error prone region with no associated phenotype and some evidence of pathogenicity and a compound heterozygous inheritance pattern" 8
addOutput $PREFIX.coding.flagged.high_qual.nopheno.predicted_pathogenic.mid.high_qual.homozygous.csv "high quality variant in a coding error prone region with no associated phenotype and some evidence of pathogenicity and a homozygous inheritance pattern" 8
addOutput $PREFIX.coding.flagged.high_qual.phenotyped.vus.recessive.csv "documented VUS of high quality in a coding error prone region with an associated phenotype and heterozygous recessive inheritance pattern" 8
addOutput $PREFIX.noncoding.flagged.high_qual.nospliceai.cadd_gt20.csv "non-coding variant of high quality in an error prone region with a high CADD score" 8
addOutput $PREFIX.coding.flagged.high_qual.nopheno.predicted_pathogenic.mid.high_qual.heterozygous.csv "high quality variant in a coding error prone region with no associated phenotype and some evidence of pathogenicity and a heterozygous inheritance pattern" 9
addOutput $PREFIX.coding.flagged.high_qual.phenotyped.predicted_pathogenic.mid.recessive.csv "high quality variant in a coding error prone region with an associated phenotype and some prediction of pathogenicity and a heterozygous recessive inheritance pattern" 9

echo "$(countVariants $finalfile) variants in final output" | tee -a $LOGFILE

tiers=( $(awk -F ',' 'NR>1{tiers[$102]+=1}END{print tiers[1],tiers[2],tiers[3],tiers[4],tiers[5],tiers[6],tiers[7]}' $finalfile) )

for i in $(seq 1 ${#tiers[@]})
do
  let loc=$i-1
  echo "Tier $i: ${tiers[$loc]} variants" | tee -a $LOGFILE
done

if [[ $CLEAN -eq 1 ]]
then
  rm -rf $PREFIX.coding*csv
  rm -rf $PREFIX.noncoding*csv
fi

exit 0