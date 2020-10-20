#!/usr/bin/env python

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

parser = argparse.ArgumentParser(description='Identification of contaminant according to taxonomic assignment of contigs into rRNA database. Script was writen by C. Pouchon (2020).')
parser.add_argument("--database", help="path to rRNA database",
                    type=str)
parser.add_argument('--out', help='taxonomy outfile', type=str)
parser.add_argument("--accessions", help="NCBI accessions to taxid for rRNA databases",
                    type=str)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

path=args.database
NCBI_accessions=args.accessions
outfile=args.out

ncbi = NCBITaxa()

accessions={}

access_dict={}
with open(NCBI_accessions) as acessf:
    for l in acessf:
        tab=l.rstrip().split('\t')
        if 'accession' in tab:
            pass
        else:
            code=tab[0]
            code_txid=tab[2]
            access_dict[code]=int(code_txid)

in_files = []
for r, d, f in os.walk(os.path.join(path,"silva_ids_acc_tax")):
    for file in f:
        if str("_accession_taxonomy.txt") in file:
            in_files.append(os.path.join(r, file))


for file in in_files:
    with open(file) as t:
        for l in t:
            toprint=list()
            line=l.rstrip().split("\t")
            access=line[0]
            ncbi_access=line[1]
            if ncbi_access in access_dict.keys():
                if len(ncbi.get_taxid_translator([access_dict[ncbi_access]]))>0:
                    sub_lineages=ncbi.get_lineage(access_dict[ncbi_access])
                    sub_names=ncbi.get_taxid_translator(sub_lineages)
                    sub_list=[val.lower() for val in list(ncbi.get_taxid_translator(sub_lineages).values())]
                else:
                    pass
                toprint.append(access)
                toprint.append(ncbi_access)
                toprint.append(",".join(sub_list))
                with open(outfile,"a+") as file:
                    file.write("\t".join(toprint)+"\n")
            else:
                pass
