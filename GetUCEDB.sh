#!/bin/bash

## Get 353 UCE db from Johnson et al. Syst. Biol. 0(0):1–13, 2018
## https://datadryad.org/resource/doi:10.5061/dryad.s3h9r6j

## Tools
MACSE=/Users/pouchonc/PhyloAlps/OrthoSkim/TOOLS/macse_v2.03.jar

## Localisation of UCE sequences
SEQ_LOC=/Users/pouchonc/Downloads/sequences/exon
RES=/Users/pouchonc/Downloads/sequences/
NAME='ref_UCE'

mkdir ${RES}Alignments

for f in ${SEQ_LOC}/*.FNA
do
  java -jar ${MACSE} -prog alignSequences -seq $f
  mv ${SEQ_LOC}/$(basename ${f%%.*})_AA.FNA ${RES}/Alignments/
  mv ${SEQ_LOC}/$(basename ${f%%.*})_NT.FNA ${RES}/Alignments/
  `dirname $0`/src/UCEfilter.py -i ${RES}/Alignments/$(basename ${f%%.*})_AA.FNA -p ${RES} -o ${NAME} -r 0.5
done
