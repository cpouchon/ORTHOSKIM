#!/bin/bash

while getopts ":m:c:p:" opt; do
  case $opt in
    m) mode="$OPTARG"
    ;;
    c) config="$OPTARG"
    ;;
    p) path="$OPTARG"
    ;;
  esac
done


## Help display of function
if [ "$1" == "-h" ]; then
    ##cat `dirname $0`/README.txt
  echo "Usage: `basename $0` [-h] [-m mode] [-c configfile.txt] [-p path]

____________________________________________________
   ___       _   _          ____  _  ___
' / _ \ _ __| |_| |__   ___/ ___|| |/ (_)_ __ ___  '
'| | | |  __| __|  _ \ / _ \___ \|   /| |  _   _  \'
'| |_| | |  | |_| | | | (_) |__) | . \| | | | | | |'
' \___/|_|   \__|_| |_|\___/____/|_|\_\_|_| |_| |_|'
____________________________________________________


OrthoSkim: Skimming of orthologous regions from shotgun sequencing
           where:
            -m (mode)  perform pipeline according to the choosen mode:
                  - [indexing] (make paramfile from chloroplastic annotated files)
                  - [SPAdes_assembly]
                  - [SPAdes_reformate] (reformates SPAdes output according to OrthoSkim pipeline)
                  - [nucleus] (perform alignment and extraction of nuclear regions from a reference list
                              and assemblies)
                  - [chloroplast] (extract regions of interest from annotation of chloroplast assemblies)
                  - [nucrdna] (extract regions of interest from annotation of rdna assemblies)
                  - [get_mitoRef] (make a bank of mitochondrial CDS as a reference)
                  - [mitochondrion] (perform alignment and extraction of mitochondrial regions from a
                                    reference list and assemblies)
                  - [stat_chloro] (output summary statistic for chloroplast assemblies)
                  - [stat_rdna] (output summary statistic for rdna assemblies)
                  - [stat_SPAdes] (check assemblies by generating summary statistics)
                  - [get_stat] (get statistics of gene extraction for files in -p {PATH} directory)
            -c (config)  set the config file
            -p (path)  set the path to get out statistic of gene extraction for [get_stat] mode
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


echo "Processing OrthoSkim program"

mkdir -p ${RES}/${PATHNAME_ASSEMBLY}/Samples

  if [ $mode == 'indexing' ]; then
		echo "INFO: mode=$mode - indexing of annotation files"
    echo "CMD: `dirname $0`/src/MakeParam.py -p ${PATH_FIND_CHLORO} -o ${RES}/${PATHNAME_ASSEMBLY}/ -cfind > ${LIST_FILES}"
		`dirname $0`/src/MakeParam.py -p ${PATH_FIND_CHLORO} -o ${RES}/${PATHNAME_ASSEMBLY}/ -cfind > ${LIST_FILES}
		echo "STATUS: done"
		exit 0

	elif [ $mode == 'SPAdes_reformate' ]; then
		echo "INFO: mode=$mode - preprocessing of assemblies files from SPAdes for nucleus/mitochondrion modes"
		while read f;
		do
			samplename=`echo ${f} | awk '{print $3}'`
			cp ${RES}/${PATHNAME_ASSEMBLY}/${samplename}/scaffolds.fasta ${RES}/${PATHNAME_ASSEMBLY}/Samples/${samplename}.fa
		done <${LIST_FILES}
		echo "STATUS: done"
		exit 0

	elif [ $mode == 'SPAdes_assembly' ]; then
		echo "INFO: mode=$mode - SPAdes assembly run"
		while read f;
		do
			samplename=`echo ${f} | awk '{print $3}'`
			r1=`echo ${f} | awk '{print $4}'`
			r2=`echo ${f} | awk '{print $5}'`
			outp=`echo ${f} | awk '{print $6}'`
      echo "${SPADES} -1 ${r1} -2 ${r2} --cov-cutoff auto -o ${outp}/${samplename} -t ${THREADS} -m ${MEMORY} --only-assembler --careful"
			${SPADES} -1 ${r1} -2 ${r2} --cov-cutoff auto -o ${outp}/${samplename} -t ${THREADS} -m ${MEMORY} --only-assembler --careful
			echo ${samplename} >> ${RES}/assembly_mode.log
		done <${LIST_FILES}
		echo "STATUS: done"
		exit 0

	elif [ $mode == 'nucleus' ]; then
		echo "INFO: mode=$mode - Extraction of nuclear genes from mapping assemblies into reference"
		mkdir ${RES}/nuc_mapping
		echo "*** make DIAMOND formatted reference database ***"
		refdb=${NUC_REF/.fasta/}
    echo "CMD: ${DIAMOND} makedb --in ${NUC_REF} -d ${refdb} --threads ${THREADS}"
		${DIAMOND} makedb --in ${NUC_REF} -d ${refdb} --threads ${THREADS}
		echo "*** mapping and extraction of assemblies into reference ***"
		for f in `find ${RES}/${PATHNAME_ASSEMBLY}/Samples/ -type f -name \*.fa`;
		do
			lib=`basename $f | perl -pe 's/.fa//'`
      echo "CMD: ${DIAMOND} blastx --outfmt 6 qseqid sseqid pident length mismatch gapopen qframe qstart qend sstart send evalue bitscore slen -d ${refdb} -q ${f} -o ${RES}/nuc_mapping/matches_${lib} --evalue ${EVALUE} --threads ${THREADS}"
			${DIAMOND} blastx --outfmt 6 qseqid sseqid pident length mismatch gapopen qframe qstart qend sstart send evalue bitscore slen -d ${refdb} -q ${f} -o ${RES}/nuc_mapping/matches_${lib} --evalue ${EVALUE} --threads ${THREADS}
			awk '{print $1}' ${RES}/nuc_mapping/matches_${lib} | sort | uniq > ${RES}/nuc_mapping/hits_${lib}
			awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $f | perl -pe 's@>@@' | awk ' NR==FNR {a[$1]=$1;next} {if($1 in a) {print ">"$1"\n"$NF}}' ${RES}/nuc_mapping/hits_${lib} - > ${RES}/nuc_mapping/contigs_hits_${lib}.fasta
      echo "CMD: ${EXONERATE} --showcigar no --showtargetgff no --showvulgar no -q ${NUC_REF} -t ${RES}/nuc_mapping/contigs_hits_${lib}.fasta --ryo ">%qi gene=%qi ref=%ti length=%tl qalnlen=%qal qbal=%qab qeal=%qae qlen=%ql alnlen=%tal baln=%tab ealn=%tae score=%s\n%tas" --showalignment no | gawk '!/^Hostname:|^Command line:|^-- completed exonerate analysis/ {print $0}' > ${RES}/nuc_mapping/out_${lib}.fasta"
      ${EXONERATE} --showcigar no --showtargetgff no --showvulgar no -q ${NUC_REF} -t ${RES}/nuc_mapping/contigs_hits_${lib}.fasta --ryo ">%qi gene=%qi ref=%ti length=%tl qalnlen=%qal qbal=%qab qeal=%qae qlen=%ql alnlen=%tal baln=%tab ealn=%tae score=%s\n%tas" --showalignment no | gawk '!/^Hostname:|^Command line:|^-- completed exonerate analysis/ {print $0}' > ${RES}/nuc_mapping/out_${lib}.fasta
      echo "CMD: `dirname $0`/src/ExoSamp.py -i ${RES}/nuc_mapping/out_${lib}.fasta -m ${mode} -o ${RES}/$mode -c ${MAXCONT} -n ${lib} -l ${MINLENGTH}"
      `dirname $0`/src/ExoSamp.py -i ${RES}/nuc_mapping/out_${lib}.fasta -m ${mode} -o ${RES}/$mode -c ${MAXCONT} -n ${lib} -l ${MINLENGTH}
			echo ${f} >> ${RES}/nucleus_mode.log
		done
		echo "STATUS: done"
		exit 0

	elif [ $mode == 'chloroplast' ]; then
		echo "INFO: mode=$mode - Extraction of chloroplastic genes from annotated files"
		while read f;
		do
			samplename=`echo ${f} | awk '{print $3}'`
			file=`echo ${f} | awk '{print $1}'`
      echo "CMD: `dirname $0`/src/AnoSamp.py -t ${file} -n ${samplename} -m "chloroplast" -g ${CHLORO_GENES} -fmt ${ANNOFMT} -o ${RES}"
			`dirname $0`/src/AnoSamp.py -t ${file} -n ${samplename} -m "chloroplast" -g ${CHLORO_GENES} -fmt ${ANNOFMT} -o ${RES}
			echo $file >> ${RES}/chloroplast_mode.log
		done <${LIST_FILES}
		echo "STATUS: done"
		exit 0

	elif [ $mode == 'nucrdna' ]; then
		echo "INFO: mode=$mode - Extraction of nuclear rDNA genes from annotated files"
		while read f;
		do
			samplename=`echo ${f} | awk '{print $3}'`
			file=`echo ${f} | awk '{print $1}' | perl -pe 's/.chloro.embl/.rdnanuc.embl/'`
      echo "CMD: `dirname $0`/src/AnoSamp.py -t ${file} -n ${samplename} -m "nucrdna" -g ${NRDNA_GENES} -fmt ${ANNOFMT} -o ${RES}"
			`dirname $0`/src/AnoSamp.py -t ${file} -n ${samplename} -m "nucrdna" -g ${NRDNA_GENES} -fmt ${ANNOFMT} -o ${RES}
			echo ${file} >> ${RES}/nucrdna_mode.log
		done <${LIST_FILES}
		echo "STATUS: done"
		exit 0

	elif [ $mode == 'get_mitoRef' ]; then
		echo "INFO: mode=$mode - Creation of mitochondria ref database"
		mkdir ${RES}/mito_ref
		find ${MITO_GENBANK_LOC} -type f -name \*.gb > ${RES}/mito_ref/listTaxa_mito
		echo '*** extraction of all CDS from annotated genbank files ***'
    echo "CMD: `dirname $0`/src/AnoRef.py -in ${RES}/mito_ref/listTaxa_mito -m mitochondrion -t CDS -fmt genbank -o ${RES}/mito_ref/refCDS -protonly"
		`dirname $0`/src/AnoRef.py -in ${RES}/mito_ref/listTaxa_mito -m mitochondrion -t CDS -fmt genbank -o ${RES}/mito_ref/refCDS -protonly
		echo '*** alignment of CDS into mito seeds, and extraction of final CDS of reference ***'
		for file in `find ${RES}/mito_ref/refCDS/ -type f -name \*.fa`;
		do
      echo "CMD: ${EXONERATE} --showcigar no --showtargetgff no --showvulgar no -q ${SEEDS_MITO} -t ${file} --ryo ">%qi gene=%qi ref=%ti length=%tl qalnlen=%qal qbal=%qab qeal=%qae qlen=%ql alnlen=%tal baln=%tab ealn=%tae score=%s\n%tas" --showalignment no | gawk '!/^Hostname:|^Command line:|^-- completed exonerate analysis/ {print $0}' > ${RES}/mito_ref/tmp.exonerate"
			${EXONERATE} --showcigar no --showtargetgff no --showvulgar no -q ${SEEDS_MITO} -t ${file} --ryo ">%qi gene=%qi ref=%ti length=%tl qalnlen=%qal qbal=%qab qeal=%qae qlen=%ql alnlen=%tal baln=%tab ealn=%tae score=%s\n%tas" --showalignment no | gawk '!/^Hostname:|^Command line:|^-- completed exonerate analysis/ {print $0}' > ${RES}/mito_ref/tmp.exonerate
      echo "CMD: `dirname $0`/src/ExoRef.py -i ${RES}/mito_ref/tmp.exonerate -m intra -p ${RES}/mito_ref/ -o refCDS_mito.fa -r ${LENGTH_RATIO}"
      `dirname $0`/src/ExoRef.py -i ${RES}/mito_ref/tmp.exonerate -m intra -p ${RES}/mito_ref/ -o refCDS_mito.fa -r ${LENGTH_RATIO}
		done
		cp ${RES}/mito_ref/refCDS_mito.fa ${MITO_REF}
		echo "STATUS: done"
		exit 0

	elif [ $mode == 'mitochondrion' ]; then
		echo 'INFO: mode=$mode - Extraction of mitchondrial genes from assemblies'
		mkdir ${RES}/mito_mapping
		refdb=${MITO_REF/.fasta/}
		echo "*** make DIAMOND formatted reference database ***"
    echo "CMD: ${DIAMOND} makedb --in ${MITO_REF} -d ${refdb} --threads ${THREADS}"
		${DIAMOND} makedb --in ${MITO_REF} -d ${refdb} --threads ${THREADS}
		echo "*** mapping and extraction of assemblies into reference ***"
		for f in `find ${RES}/${PATHNAME_ASSEMBLY}/Samples/ -type f -name \*.fa`;
		do
			lib=`basename $f | perl -pe 's/.fa//'`
      echo "CMD: ${DIAMOND} blastx --outfmt 6 qseqid sseqid pident length mismatch gapopen qframe qstart qend sstart send evalue bitscore slen -d ${refdb} -q ${f} -o ${RES}/mito_mapping/matches_${lib} --evalue ${EVALUE} --threads ${THREADS}"
			${DIAMOND} blastx --outfmt 6 qseqid sseqid pident length mismatch gapopen qframe qstart qend sstart send evalue bitscore slen -d ${refdb} -q ${f} -o ${RES}/mito_mapping/matches_${lib} --evalue ${EVALUE} --threads ${THREADS}
			awk '{print $1}' ${RES}/mito_mapping/matches_${lib} | sort | uniq > ${RES}/mito_mapping/hits_${lib}
			awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $f | perl -pe 's@>@@' | awk ' NR==FNR {a[$1]=$1;next} {if($1 in a) {print ">"$1"\n"$NF}}' ${RES}/mito_mapping/hits_${lib} - > ${RES}/mito_mapping/contigs_hits_${lib}.fasta
      echo "CMD: ${EXONERATE} --showcigar no --showtargetgff no --showvulgar no -q ${MITO_REF} -t ${RES}/mito_mapping/contigs_hits_${lib}.fasta --ryo ">%qi gene=%qi ref=%ti length=%tl qalnlen=%qal qbal=%qab qeal=%qae qlen=%ql alnlen=%tal baln=%tab ealn=%tae score=%s\n%tas" --showalignment no | gawk '!/^Hostname:|^Command line:|^-- completed exonerate analysis/ {print $0}' > ${RES}/mito_mapping/out_${lib}.fasta"
			${EXONERATE} --showcigar no --showtargetgff no --showvulgar no -q ${MITO_REF} -t ${RES}/mito_mapping/contigs_hits_${lib}.fasta --ryo ">%qi gene=%qi ref=%ti length=%tl qalnlen=%qal qbal=%qab qeal=%qae qlen=%ql alnlen=%tal baln=%tab ealn=%tae score=%s\n%tas" --showalignment no | gawk '!/^Hostname:|^Command line:|^-- completed exonerate analysis/ {print $0}' > ${RES}/mito_mapping/out_${lib}.fasta
      echo "CMD: `dirname $0`/src/ExoSamp.py -i ${RES}/mito_mapping/out_${lib}.fasta -m $mode -o ${RES}/$mode -c ${MAXCONT} -n ${lib} -l ${MINLENGTH}"
			`dirname $0`/src/ExoSamp.py -i ${RES}/mito_mapping/out_${lib}.fasta -m $mode -o ${RES}/$mode -c ${MAXCONT} -n ${lib} -l ${MINLENGTH}
			rm ${RES}/mito_mapping/contigs_hits_${lib}.fasta ${RES}/mito_mapping/matches_${lib} ${RES}/mito_mapping/hits_${lib} ${RES}/mito_mapping/out_${lib}.fasta
			echo ${f} >> ${RES}/mitochondrion_mode.log
		done
		echo "STATUS: done"
		exit 0


	elif [ $mode == 'stat_chloro' ]; then
		echo "INFO: mode=$mode - Summary statistics of chloroplast assemblies"
		${QUAST} `awk '{print $1}' ${LIST_FILES} | perl -pe 's/.embl/.fasta/'  | perl -pe 's/\n/ /'` -o ${RES}/report_chloro_assemblies
		rm -rf ${RES}/report_chloro_assemblies/basic_stats/
		#rm -rf ${RES}/report_chloro_assemblies/icarus_viewers/
		rm ${RES}/report_chloro_assemblies/report.tex ${RES}/report_chloro_assemblies/report.tsv ${RES}/report_chloro_assemblies/report.txt ${RES}/report_chloro_assemblies/transposed_report.tex ${RES}/report_chloro_assemblies/transposed_report.tsv
		echo "STATUS: done"
		exit 0

	elif [ $mode == 'stat_rdna' ]; then
		echo "INFO: mode=$mode - Summary statistics of rdna assemblies"
		${QUAST} `awk '{print $1}' ${LIST_FILES} | perl -pe 's/.chloro.embl/.rdnanuc.fasta/'  | perl -pe 's/\n/ /'` -o ${RES}/report_rdnanuc_assemblies
		rm -rf ${RES}/report_rdnanuc_assemblies/basic_stats/
		#rm -rf ${RES}/report_rdnanuc_assemblies/icarus_viewers/
		rm ${RES}/report_rdnanuc_assemblies/report.tex ${RES}/report_rdnanuc_assemblies/report.tsv ${RES}/report_rdnanuc_assemblies/report.txt ${RES}/report_rdnanuc_assemblies/transposed_report.tex ${RES}/report_rdnanuc_assemblies/transposed_report.tsv
		echo "STATUS: done"
		exit 0

	elif [ $mode == 'stat_SPAdes' ]; then
		echo "INFO: mode=$mode - Summary statistics of nuclear assemblies"
		${QUAST} ${RES}/${PATHNAME_ASSEMBLY}/Samples/*.fa -o ${RES}/report_SPAdes_assemblies
		rm -rf ${RES}/report_SPAdes_assemblies/basic_stats/
		#rm -rf ${RES}/report_SPAdes_assemblies/icarus_viewers/
		rm ${RES}/report_SPAdes_assemblies/report.tex ${RES}/report_SPAdes_assemblies/report.tsv ${RES}/report_SPAdes_assemblies/report.txt ${RES}/report_SPAdes_assemblies/transposed_report.tex ${RES}/report_SPAdes_assemblies/transposed_report.tsv
		echo "STATUS: done"
		exit 0

  elif [ $mode == 'get_stat' ]; then
    echo "INFO: mode=$mode - Statistics of genes extraction according to given path"
    echo "CMD: `dirname $0`/src/ExoStat.py -p $path -pfind > $path/report.log"
    `dirname $0`/src/ExoStat.py -p $path -pfind > $path/report.log
    echo "STATUS: done"
    exit 0

	fi

exit 0
