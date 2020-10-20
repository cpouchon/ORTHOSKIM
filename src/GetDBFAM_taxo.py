#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
#from Bio.Alphabet import *
from Bio.SeqRecord import *
#from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA;
from Bio.Seq import *
from Bio.SeqUtils import *
from joblib import Parallel, delayed
import multiprocessing
import glob
import argparse
import sys
import os, errno

from ete3 import NCBITaxa


# to write NCBI accessions from SILVA
# if [ -s `dirname $0`/ressources/Silva_NCBI_accessions_to_taxid.txt ]; then
#     echo ""
# else
#     echo "extraction of NCBI accessions IDs from SILVA database"
#     for f in ${RNA_CONTA_DB}/silva_ids_acc_tax/silva-*.txt;
#     do
#       awk '{print $2}' $f >> ${RNA_CONTA_DB}/silva_ids_acc_tax/access.tmp
#     done
#     sort ${RNA_CONTA_DB}/silva_ids_acc_tax/access.tmp | uniq > ${RNA_CONTA_DB}/silva_ids_acc_tax/accessions_NCBI.txt
#     awk 'NR==FNR{seen[$1]; next} $1 in seen' ${RNA_CONTA_DB}/silva_ids_acc_tax/accessions_NCBI.txt ${BLAST_NT_ACCESSION_TAXID} >> `dirname $0`/ressources/Silva_NCBI_accessions_to_taxid.txt
#     awk 'NR==FNR{seen[$1]; next} $1 in seen' ${RNA_CONTA_DB}/silva_ids_acc_tax/accessions_NCBI.txt ${BLAST_NT_WGS_ACCESSION_TAXID} >> `dirname $0`/ressources/Silva_NCBI_accessions_to_taxid.txt
#     rm ${RNA_CONTA_DB}/silva_ids_acc_tax/accessions_NCBI.txt ${RNA_CONTA_DB}/silva_ids_acc_tax/access.tmp
# fi
#

#grep ">" rfam-5.8s-database-id98.fasta | awk '{split($0,a,">"); print a[2]}' | awk '{split($0,a,"/"); print a[1]}' >> access_rfam.tmp
#grep ">" rfam-5s-database-id98.fasta | awk '{split($0,a,">"); print a[2]}' | awk '{split($0,a,"/"); print a[1]}' >> access_rfam.tmp

parser = argparse.ArgumentParser(description='Export the NCBI taxonomy of DBFAM annotations. Script was writen by C. Pouchon (2020).')
parser.add_argument('--outdir', help='sequence names of rfam', type=str)
parser.add_argument("--path", help="path to DBFAM accessions",
                    type=str)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

class ContinueR(Exception):
    pass

outp=args.outdir
path = args.path

in_files = []
for r, d, f in os.walk(path):
    for file in f:
        if 'DBFAM' in file:
            if os.path.join(r, file) not in in_files:
                in_files.append(os.path.join(r,file))

ncbi = NCBITaxa()
taxo_list=list()
continue_record = ContinueR()

for f in in_files:
    seq=dict()
    outf=os.path.basename(f).split(".")[0]+".fasta"
    cur_genome = SeqIO.parse(f, "genbank")
    for record in cur_genome:
        try:
            for feat in record.features:
                if feat.type == 'source':
                    org=feat.qualifiers['organism'][0].split(" ")
                    db_info=feat.qualifiers['db_xref']
                    for subp in db_info:
                        if 'taxon:' in subp:
                            taxid=subp.split('taxon:')[1]
                        else:
                            raise continue_record
            sequence=str(record.seq)
            seqid=str(record.id)
            accession=str(record.name)

            if len(ncbi.get_taxid_translator([taxid]))>0:
                sub_lineages=ncbi.get_lineage(int(taxid))
                sub_names=ncbi.get_taxid_translator(sub_lineages)
                sub_list=[val.lower() for val in list(ncbi.get_taxid_translator(sub_lineages).values())]
                seq[seqid]=sequence
                toprint=list()
                toprint.append(seqid)
                toprint.append(accession)
                toprint.append(",".join(sub_list))
                taxo_list.append(toprint)
            else:
                if len(ncbi.get_name_translator([org[0]]))>0:
                    newtaxid=list(ncbi.get_name_translator([org[0]]).values())[0][0]
                    sub_lineages=ncbi.get_lineage(newtaxid)
                    sub_names=ncbi.get_taxid_translator(sub_lineages)
                    sub_list=[val.lower() for val in list(ncbi.get_taxid_translator(sub_lineages).values())]
                    seq[seqid]=sequence
                    toprint=list()
                    toprint.append(seqid)
                    toprint.append(accession)
                    toprint.append(",".join(sub_list))
                    taxo_list.append(toprint)
                else:
                    pass

        except ContinueR:
            continue

    for s in seq.keys():
        header=">"+str(s)
        if os.path.isfile(os.path.join(outp, outf)):
            with open(os.path.join(outp, outf), 'a+') as file:
                old_headers = []
                end_file=file.tell()
                file.seek(0)
                for line in file:
                    if line.startswith(">"):
                        old_headers.append(line.rstrip())
                if not header in old_headers:
                    file.seek(end_file)
                    file.write(header+'\n')
                    file.write(str(seq[s])+'\n')
                else:
                    pass
        else :
            with open(os.path.join(outp, outf), 'w') as out:
                out.write(header+'\n')
                out.write(str(seq[s])+'\n')

with open(os.path.join(outp, "DBFAM_NCBI_taxonomy.txt"),'a+') as taxfile:
    for l in taxo_list:
        taxfile.write("\t".join(l)+"\n")
