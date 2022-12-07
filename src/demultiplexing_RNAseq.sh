#!/bin/bash

echo "Merging Illumina lanes into one file per sample and read (for paired reads)
Using simple cat function"

for fn in OV40{1..9} OV41{0..3}
do
  lanes=`echo ${fn}*`
  i=1
    for lane in $lanes
  do
    echo "Processing lane $lane"
    if [[ "$lane" == *"L001"* ]];
    then
      temp=`ls ${lane}/*R1*`
      sample=`echo $(basename "$temp")`
      sample=`cut -d "_" -f 1 <<<  "$sample"`
      if [[ -d $sample ]];
      then
       rm -r $sample
      fi
      mkdir -v $sample
    fi

    echo "concatinating sample ${sample}"
    gunzip $lane/*
    cat -v $lane/*R1*.fastq >> $sample/${sample}_R1.fastq
    cat -v $lane/*R2*.fastq >> $sample/${sample}_R2.fastq
  done
done