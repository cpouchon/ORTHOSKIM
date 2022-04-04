#!/bin/bash


while getopts ":m:c:s:" opt; do
  case $opt in
    m) mode="$OPTARG"
    ;;
    c) config="$OPTARG"
    ;;
    s) seeds+=("$OPTARG")
    ;;
  esac
done


## Help display of function
if [ "$1" == "-h" ]; then
    ##cat `dirname $0`/README.txt
  echo "Usage: `basename $0` [-h] [-m mode] [-c config_file.txt]

Extraction of targeted sequences from genomic annotations following the ORTHOSKIM architecture
           where:
            -m (mode) perform extraction according to the chosen mode:
                  - chloroplast, mitochondrion, nucrdna
            -c (config_file) set the ORTHOSKIM parameter file
            -s (rdna seeds) rdna seeds used including the ITS1 and ITS2 sequences
            "
  exit 0
fi


## OrthoSkim process
source $config
source `dirname $0`/tools.sh

set -e

if [ ${VERBOSE} -eq 1 ]; then
      set -x
fi

mkdir -p ${RES}/Assembly/Samples
mkdir -p ${RES}/Extraction

if [ $mode == 'chloroplast' ]; then
  echo "INFO: Extraction of chloroplast sequences from annotations"
  start=$(date +%s)
  mkdir -p ${RES}/tmp/unclean
  start=$(date +%s)
  echo '*** extraction of all genes from annotations ***'
  echo "CMD: `dirname $0`/src/Ano_extraction.py --single -in ${CHLORO_ANNOTATIONS} -o ${RES}/tmp/unclean -m chloroplast -fmt ${CHLORO_DB_FMT} --codon `dirname $0`/resources/tRNA_codons.tab"
  `dirname $0`/src/Ano_extraction.py --single -in ${CHLORO_ANNOTATIONS} -o ${RES}/tmp/unclean -m chloroplast -fmt ${CHLORO_DB_FMT} --codon `dirname $0`/resources/tRNA_codons.tab

  echo '*** alignment of CDS into seeds and extraction of final sequences ***'
  for file in `find ${RES}/tmp/unclean/chloroplast_CDS/ -type f -name \*.fa`;
  do
  	echo "CMD: ${EXONERATE} --model protein2genome -s ${EXO_SCORE} -q ${SEEDS_CHLORO_CDS} -t ${file} --showquerygff yes --showtargetgff yes --showvulgar no --showcigar no --showalignment no | awk '!/^Hostname:|^Command line:|^-- completed exonerate analysis|#/ {print \$0}' | awk '/Target/ { if (ok) print entry; ok=1;entry=\$0;next;}{entry = entry \"\n\" \$0;}/intron\t-/{ok=0;}END{if (ok) print entry;}' > ${RES}/tmp/tmp.exonerate"
  	${EXONERATE} --model protein2genome -s ${EXO_SCORE} -q ${SEEDS_CHLORO_CDS} -t ${file} --showquerygff yes --showtargetgff yes --showvulgar no --showcigar no --showalignment no | awk '!/^Hostname:|^Command line:|^-- completed exonerate analysis|#/ {print $0}' | awk '/Target/ { if (ok) print entry; ok=1;entry=$0;next;}{entry = entry "\n" $0;}/intron\t-/{ok=0;}END{if (ok) print entry;}' > ${RES}/tmp/tmp.exonerate
  	echo "CMD: `dirname $0`/src/ExoComp.py -i ${file} -s ${SEEDS_CHLORO_CDS} -g ${RES}/tmp/tmp.exonerate -m chloroplast_CDS -o ${RES}/Extraction --threads ${THREADS}"
  	`dirname $0`/src/ExoComp.py -i ${file} -s ${SEEDS_CHLORO_CDS} -g ${RES}/tmp/tmp.exonerate -m chloroplast_CDS -o ${RES}/Extraction/ --threads ${THREADS}
  done
  rm -rf ${RES}/tmp/unclean/chloroplast_CDS/

  echo '*** alignment of rRNA into seeds and extraction of final sequences ***'
  for file in `find ${RES}/tmp/unclean/chloroplast_rRNA/ -type f -name \*.fa`;
  do
  	echo "CMD: ${EXONERATE} --model genome2genome -s ${EXO_SCORE} -q ${SEEDS_CHLORO_rRNA} -t ${file} --showquerygff yes --showtargetgff yes --showvulgar no --showcigar no --showalignment no | awk '!/^Hostname:|^Command line:|^-- completed exonerate analysis|#/ {print \$0}' | awk '/Target/ { if (ok) print entry; ok=1;entry=\$0;next;}{entry = entry \"\n\" \$0;}/intron\t-/{ok=0;}END{if (ok) print entry;}' > ${RES}/tmp/tmp.exonerate"
  	${EXONERATE} --model genome2genome -s ${EXO_SCORE} -q ${SEEDS_CHLORO_rRNA} -t ${file} --showquerygff yes --showtargetgff yes --showvulgar no --showcigar no --showalignment no | awk '!/^Hostname:|^Command line:|^-- completed exonerate analysis|#/ {print $0}' | awk '/Target/ { if (ok) print entry; ok=1;entry=$0;next;}{entry = entry "\n" $0;}/intron\t-/{ok=0;}END{if (ok) print entry;}' > ${RES}/tmp/tmp.exonerate
  	echo "CMD: `dirname $0`/src/ExoComp.py -i ${file} -s ${SEEDS_CHLORO_rRNA} -g ${RES}/tmp/tmp.exonerate -m chloroplast_rRNA -o ${RES}/Extraction --threads ${THREADS}"
  	`dirname $0`/src/ExoComp.py -i ${file} -s ${SEEDS_CHLORO_rRNA} -g ${RES}/tmp/tmp.exonerate -m chloroplast_rRNA -o ${RES}/Extraction/ --threads ${THREADS}
  done
  rm -rf ${RES}/tmp/unclean/chloroplast_rRNA/

  echo '*** alignment of tRNA into seeds and extraction of final sequences ***'
  for file in `find ${RES}/tmp/unclean/chloroplast_tRNA/ -type f -name \*.fa`;
  do
  	echo "CMD: ${EXONERATE} --model genome2genome -s ${EXO_SCORE} -q ${SEEDS_CHLORO_tRNA} -t ${file} --showquerygff yes --showtargetgff yes --showvulgar no --showcigar no --showalignment no | awk '!/^Hostname:|^Command line:|^-- completed exonerate analysis|#/ {print \$0}' | awk '/Target/ { if (ok) print entry; ok=1;entry=\$0;next;}{entry = entry \"\n\" \$0;}/intron\t-/{ok=0;}END{if (ok) print entry;}' > ${RES}/tmp/tmp.exonerate"
  	${EXONERATE} --model genome2genome -s ${EXO_SCORE} -q ${SEEDS_CHLORO_tRNA} -t ${file} --showquerygff yes --showtargetgff yes --showvulgar no --showcigar no --showalignment no | awk '!/^Hostname:|^Command line:|^-- completed exonerate analysis|#/ {print $0}' | awk '/Target/ { if (ok) print entry; ok=1;entry=$0;next;}{entry = entry "\n" $0;}/intron\t-/{ok=0;}END{if (ok) print entry;}' > ${RES}/tmp/tmp.exonerate
  	echo "CMD: `dirname $0`/src/ExoComp.py -i ${file} -s ${SEEDS_CHLORO_tRNA} -g ${RES}/tmp/tmp.exonerate -m chloroplast_tRNA -o ${RES}/Extraction --threads ${THREADS}"
  	`dirname $0`/src/ExoComp.py -i ${file} -s ${SEEDS_CHLORO_tRNA} -g ${RES}/tmp/tmp.exonerate -m chloroplast_tRNA -o ${RES}/Extraction/ --threads ${THREADS}
  done
  rm -rf ${RES}/tmp/unclean/chloroplast_tRNA/
  end=`date +%s`
  runtime=$( echo "$end - $start" | bc -l )
  hours=$((runtime / 3600))
  minutes=$(( (runtime % 3600) / 60 ))
  seconds=$(( (runtime % 3600) % 60 ))
  echo "INFO: time ellapsed for $mode step: $hours:$minutes:$seconds (hh:mm:ss)"
  echo "STATUS: done"
  echo ""
  exit 0

elif [ $mode == 'mitochondrion' ]; then
  echo "INFO: Extraction of mitochondrion sequences from annotations"
  start=$(date +%s)
  mkdir -p ${RES}/tmp/unclean
  start=$(date +%s)
  echo '*** extraction of all genes from annotations ***'
  echo "CMD: `dirname $0`/src/Ano_extraction.py --single -in ${MITO_ANNOTATIONS} -o ${RES}/tmp/unclean -m mitochondrion -fmt ${MITO_DB_FMT} --codon `dirname $0`/resources/tRNA_codons.tab"
  `dirname $0`/src/Ano_extraction.py --single -in ${MITO_ANNOTATIONS} -o ${RES}/tmp/unclean -m mitochondrion -fmt ${MITO_DB_FMT} --codon `dirname $0`/resources/tRNA_codons.tab

  echo '*** alignment of CDS into seeds and extraction of final sequences ***'
  for file in `find ${RES}/tmp/unclean/mitochondrion_CDS/ -type f -name \*.fa`;
  do
  	echo "CMD: ${EXONERATE} --model protein2genome -s ${EXO_SCORE} -q ${SEEDS_MITO_CDS} -t ${file} --showquerygff yes --showtargetgff yes --showvulgar no --showcigar no --showalignment no | awk '!/^Hostname:|^Command line:|^-- completed exonerate analysis|#/ {print \$0}' | awk '/Target/ { if (ok) print entry; ok=1;entry=\$0;next;}{entry = entry \"\n\" \$0;}/intron\t-/{ok=0;}END{if (ok) print entry;}' > ${RES}/tmp/tmp.exonerate"
  	${EXONERATE} --model protein2genome -s ${EXO_SCORE} -q ${SEEDS_MITO_CDS} -t ${file} --showquerygff yes --showtargetgff yes --showvulgar no --showcigar no --showalignment no | awk '!/^Hostname:|^Command line:|^-- completed exonerate analysis|#/ {print $0}' | awk '/Target/ { if (ok) print entry; ok=1;entry=$0;next;}{entry = entry "\n" $0;}/intron\t-/{ok=0;}END{if (ok) print entry;}' > ${RES}/tmp/tmp.exonerate
  	echo "CMD: `dirname $0`/src/ExoComp.py -i ${file} -s ${SEEDS_MITO_CDS} -g ${RES}/tmp/tmp.exonerate -m mitochondrion_CDS -o ${RES}/Extraction --threads ${THREADS}"
  	`dirname $0`/src/ExoComp.py -i ${file} -s ${SEEDS_MITO_CDS} -g ${RES}/tmp/tmp.exonerate -m mitochondrion_CDS -o ${RES}/Extraction/ --threads ${THREADS}
  done
  rm -rf ${RES}/tmp/unclean/mitochondrion_CDS/

  echo '*** alignment of rRNA into seeds and extraction of final sequences ***'
  for file in `find ${RES}/tmp/unclean/mitochondrion_rRNA/ -type f -name \*.fa`;
  do
  	echo "CMD: ${EXONERATE} --model genome2genome -s ${EXO_SCORE} -q ${SEEDS_MITO_rRNA} -t ${file} --showquerygff yes --showtargetgff yes --showvulgar no --showcigar no --showalignment no | awk '!/^Hostname:|^Command line:|^-- completed exonerate analysis|#/ {print \$0}' | awk '/Target/ { if (ok) print entry; ok=1;entry=\$0;next;}{entry = entry \"\n\" \$0;}/intron\t-/{ok=0;}END{if (ok) print entry;}' > ${RES}/tmp/tmp.exonerate"
  	${EXONERATE} --model genome2genome -s ${EXO_SCORE} -q ${SEEDS_MITO_rRNA} -t ${file} --showquerygff yes --showtargetgff yes --showvulgar no --showcigar no --showalignment no | awk '!/^Hostname:|^Command line:|^-- completed exonerate analysis|#/ {print $0}' | awk '/Target/ { if (ok) print entry; ok=1;entry=$0;next;}{entry = entry "\n" $0;}/intron\t-/{ok=0;}END{if (ok) print entry;}' > ${RES}/tmp/tmp.exonerate
  	echo "CMD: `dirname $0`/src/ExoComp.py -i ${file} -s ${SEEDS_MITO_rRNA} -g ${RES}/tmp/tmp.exonerate -m mitochondrion_rRNA -o ${RES}/Extraction --threads ${THREADS}"
  	`dirname $0`/src/ExoComp.py -i ${file} -s ${SEEDS_MITO_rRNA} -g ${RES}/tmp/tmp.exonerate -m mitochondrion_rRNA -o ${RES}/Extraction/ --threads ${THREADS}
  done
  rm -rf ${RES}/tmp/unclean/mitochondrion_rRNA/
  end=`date +%s`
  runtime=$( echo "$end - $start" | bc -l )
  hours=$((runtime / 3600))
  minutes=$(( (runtime % 3600) / 60 ))
  seconds=$(( (runtime % 3600) % 60 ))
  echo "INFO: time ellapsed for $mode step: $hours:$minutes:$seconds (hh:mm:ss)"
  echo "STATUS: done"
  echo ""
  exit 0

elif [ $mode == 'nucrdna' ]; then
  echo "INFO: Extraction of rdna sequences from annotations"
  start=$(date +%s)
  mkdir -p ${RES}/tmp/unclean
  start=$(date +%s)
  rdna_seeds=$seeds
  echo '*** extraction of all genes from annotations ***'
  echo "CMD: `dirname $0`/src/Ano_extraction_nucrdna.py --single -in ${NRDNA_ANNOTATIONS} -o ${RES}/tmp/unclean -m nucrdna -fmt ${NRDNA_DB_FMT}"
  `dirname $0`/src/Ano_extraction_nucrdna.py --single -in ${NRDNA_ANNOTATIONS} -o ${RES}/tmp/unclean -m nucrdna -fmt ${NRDNA_DB_FMT}

  echo '*** alignment of sequences into seeds and computing final sequences ***'
  for file in `find ${RES}/tmp/unclean/nucrdna/ -type f -name \*.fa`;
  do
  	echo "CMD: ${EXONERATE} --model genome2genome -s ${EXO_SCORE} -q ${rdna_seeds} -t ${file} --showquerygff yes --showtargetgff yes --showvulgar no --showcigar no --showalignment no | awk '!/^Hostname:|^Command line:|^-- completed exonerate analysis|#/ {print \$0}' | awk '/Target/ { if (ok) print entry; ok=1;entry=\$0;next;}{entry = entry \"\n\" \$0;}/intron\t-/{ok=0;}END{if (ok) print entry;}' > ${RES}/tmp/tmp.exonerate"
  	${EXONERATE} --model genome2genome -s ${EXO_SCORE} -q ${rdna_seeds} -t ${file} --showquerygff yes --showtargetgff yes --showvulgar no --showcigar no --showalignment no | awk '!/^Hostname:|^Command line:|^-- completed exonerate analysis|#/ {print $0}' | awk '/Target/ { if (ok) print entry; ok=1;entry=$0;next;}{entry = entry "\n" $0;}/intron\t-/{ok=0;}END{if (ok) print entry;}' > ${RES}/tmp/tmp.exonerate
  	echo "CMD: `dirname $0`/src/ExoComp.py -i ${file} -s ${rdna_seeds} -g ${RES}/tmp/tmp.exonerate -m nucrdna -o ${RES}/Extraction --threads ${THREADS}"
  	`dirname $0`/src/ExoComp.py -i ${file} -s ${rdna_seeds} -g ${RES}/tmp/tmp.exonerate -m nucrdna -o ${RES}/Extraction/ --threads ${THREADS}
  done
  rm -rf ${RES}/tmp/unclean/nucrdna/
  end=`date +%s`
  runtime=$( echo "$end - $start" | bc -l )
  hours=$((runtime / 3600))
  minutes=$(( (runtime % 3600) / 60 ))
  seconds=$(( (runtime % 3600) % 60 ))
  echo "INFO: time ellapsed for $mode step: $hours:$minutes:$seconds (hh:mm:ss)"
  echo "STATUS: done"
  echo ""
  exit 0

fi

exit 0
