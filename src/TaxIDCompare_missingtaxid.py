#!/usr/bin/env python

import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.Alphabet import *
from Bio.SeqRecord import *
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA;
from Bio.Seq import *
from Bio.SeqUtils import *
from ete3 import NCBITaxa

import argparse
import sys
import os, errno



parser = argparse.ArgumentParser(description='Check errors according to taxonomic assignment of best hits and query taxid. Script was writen by C. Pouchon (2020).')
parser.add_argument('-b', '--barcodes', help='list of barcodes', type=str)
parser.add_argument('--families', help='list of families for new taxid', type=str)
parser.add_argument("--local", help="use local corresponding for query taxid and families",
                    action="store_true")
parser.add_argument("-o","--outdir", help="out directory path",
                    type=str)
parser.add_argument("-t","--taxa", help="list of taxa to include in the final concatenated file",
                    type=str)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()


#list_accessions = args.accession
outpath=args.outdir
giventaxa = args.taxa
list_newtaxid=args.families

ncbi = NCBITaxa()

# access_dict={}
# with open(list_accessions) as f:
#     for l in f:
#         tab=l.rstrip().split('\t')
#         if 'accession' in tab:
#             pass
#         else:
#             code=tab[0]
#             code_txid=tab[1]
#             access_dict[code]=int(code_txid)

if args.local:
    list_newtaxid=args.families
    missing_families={}
    with open(list_newtaxid) as f:
        for l in f:
            tab=l.rstrip().split(' ')
            missing_families[tab[0]]=tab[1]

list_barcodes = [str(item) for item in args.barcodes.split(',')]

taxa_dict={}
with open(giventaxa) as f:
    for l in f:
        taxa_dict.setdefault(l.rstrip(), {})
        for b in list_barcodes:
            taxa_dict[l.rstrip()][b]="NA"

for b in list_barcodes:
    with open(str(outpath+"/"+b+"_besthits"+".out")) as out:
        for l in out:
            sub_error=0
            tab=l.rstrip().split('\t')
            query_name=tab[0].replace(";","")
            sub_name=tab[1]
            taxid=int(tab[2].split(";")[0])
            if len(ncbi.get_taxid_translator([taxid]))==0:
                sub_error=sub_error+1
            else:
                sub_lineages=ncbi.get_lineage(taxid)
                sub_names=ncbi.get_taxid_translator(sub_lineages)
                sub_ranks=ncbi.get_rank(sub_names.keys())
                list_fam=[key  for (key, value) in sub_ranks.items() if value == u'family']
                if len(list_fam)==0:
                    res="NA"
                    with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                        file.write("No family rank for taxid on subject: "+sub_name+'\n')
                        file.close()
                else:
                    sub_fam=str(sub_names[[key  for (key, value) in sub_ranks.items() if value == u'family'][0]])

            if len([int(i) for i in query_name.split("_") if i.isdigit()])==0:
                res="NA"
                with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                    file.write("No taxid on query: "+query_name+'\n')
                    file.close()
            else:
                query_taxid=[int(i) for i in query_name.split("_") if i.isdigit()][0]
                if sub_error==0:
                    if args.local:
                        if str(query_taxid) in missing_families.keys():
                            query_fam=missing_families[str(query_taxid)]
                            if len(query_fam)==0:
                                res="NA"
                                with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                    file.write("family rank not found in NCBI taxonomy for query: "+query_name+'\n')
                                    file.close()
                            elif len(sub_fam)==0:
                                res="NA"
                                with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                    file.write("family rank not found in NCBI taxonomy for subject: "+sub_name+'\n')
                                    file.close()
                            else:
                                if query_fam == sub_fam:
                                    res="TRUE"
                                else:
                                    res="FALSE"
                        else:
                            if len(ncbi.get_taxid_translator([int(query_taxid)]))==0:
                                if len(ncbi.get_name_translator([" ".join(query_name.split("_")[0:2])]))==0:
                                    if len(ncbi.get_name_translator([" ".join(query_name.split("_")[0:1])]))==0:
                                        res="NA"
                                        with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                            file.write("taxid not found in NCBI taxonomy for query: "+query_name+'\n')
                                            file.close()
                                    else:
                                        species_name=" ".join(query_name.split("_")[0:1])
                                        query_taxid=int(ncbi.get_name_translator([species_name])[species_name][0])
                                        query_lineages=ncbi.get_lineage(query_taxid)
                                        query_names = ncbi.get_taxid_translator(query_lineages)
                                        query_ranks=ncbi.get_rank(query_names.keys())
                                        query_fam=str(query_names[[key  for (key, value) in query_ranks.items() if value == u'family'][0]])
                                        if len(query_fam)==0:
                                            res="NA"
                                            with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                                file.write("family rank not found in NCBI taxonomy for query: "+query_name+'\n')
                                                file.close()
                                        elif len(sub_fam)==0:
                                            res="NA"
                                            with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                                file.write("family rank not found in NCBI taxonomy for subject: "+sub_name+'\n')
                                                file.close()
                                        else:
                                            if query_fam == sub_fam:
                                                res="TRUE"
                                            else:
                                                res="FALSE"
                                else:
                                    species_name=" ".join(query_name.split("_")[0:2])
                                    query_taxid=int(ncbi.get_name_translator([species_name])[species_name][0])
                                    query_lineages=ncbi.get_lineage(query_taxid)
                                    query_names = ncbi.get_taxid_translator(query_lineages)
                                    query_ranks=ncbi.get_rank(query_names.keys())
                                    query_fam=str(query_names[[key  for (key, value) in query_ranks.items() if value == u'family'][0]])
                                    if len(query_fam)==0:
                                        res="NA"
                                        with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                            file.write("family rank not found in NCBI taxonomy for query: "+query_name+'\n')
                                            file.close()
                                    elif len(sub_fam)==0:
                                        res="NA"
                                        with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                            file.write("family rank not found in NCBI taxonomy for subject: "+sub_name+'\n')
                                            file.close()
                                    else:
                                        if query_fam == sub_fam:
                                            res="TRUE"
                                        else:
                                            res="FALSE"
                            else:
                                query_lineages=ncbi.get_lineage(query_taxid)
                                query_names = ncbi.get_taxid_translator(query_lineages)
                                query_ranks=ncbi.get_rank(query_names.keys())
                                query_fam=str(query_names[[key  for (key, value) in query_ranks.items() if value == u'family'][0]])
                                if len(query_fam)==0:
                                    res="NA"
                                    with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                        file.write("family rank not found in NCBI taxonomy for query: "+query_name+'\n')
                                        file.close()
                                elif len(sub_fam)==0:
                                    res="NA"
                                    with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                        file.write("family rank not found in NCBI taxonomy for subject: "+sub_name+'\n')
                                        file.close()
                                else:
                                    if query_fam == sub_fam:
                                        res="TRUE"
                                    else:
                                        res="FALSE"
                    else:
                        if len(ncbi.get_taxid_translator([int(query_taxid)]))==0:
                            if len(ncbi.get_name_translator([" ".join(query_name.split("_")[0:2])]))==0:
                                if len(ncbi.get_name_translator([" ".join(query_name.split("_")[0:1])]))==0:
                                    res="NA"
                                    with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                        file.write("taxid not found in NCBI taxonomy for query: "+query_name+'\n')
                                        file.close()
                                else:
                                    species_name=" ".join(query_name.split("_")[0:1])
                                    query_taxid=int(ncbi.get_name_translator([species_name])[species_name][0])
                                    query_lineages=ncbi.get_lineage(query_taxid)
                                    query_names = ncbi.get_taxid_translator(query_lineages)
                                    query_ranks=ncbi.get_rank(query_names.keys())
                                    query_fam=str(query_names[[key  for (key, value) in query_ranks.items() if value == u'family'][0]])
                                    if len(query_fam)==0:
                                        res="NA"
                                        with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                            file.write("family rank not found in NCBI taxonomy for query: "+query_name+'\n')
                                            file.close()
                                    elif len(sub_fam)==0:
                                        res="NA"
                                        with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                            file.write("family rank not found in NCBI taxonomy for subject: "+sub_name+'\n')
                                            file.close()
                                    else:
                                        if query_fam == sub_fam:
                                            res="TRUE"
                                        else:
                                            res="FALSE"
                            else:
                                species_name=" ".join(query_name.split("_")[0:2])
                                query_taxid=int(ncbi.get_name_translator([species_name])[species_name][0])
                                query_lineages=ncbi.get_lineage(query_taxid)
                                query_names = ncbi.get_taxid_translator(query_lineages)
                                query_ranks=ncbi.get_rank(query_names.keys())
                                query_fam=str(query_names[[key  for (key, value) in query_ranks.items() if value == u'family'][0]])
                                if len(query_fam)==0:
                                    res="NA"
                                    with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                        file.write("family rank not found in NCBI taxonomy for query: "+query_name+'\n')
                                        file.close()
                                elif len(sub_fam)==0:
                                    res="NA"
                                    with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                        file.write("family rank not found in NCBI taxonomy for subject: "+sub_name+'\n')
                                        file.close()
                                else:
                                    if query_fam == sub_fam:
                                        res="TRUE"
                                    else:
                                        res="FALSE"
                        else:
                            query_lineages=ncbi.get_lineage(query_taxid)
                            query_names = ncbi.get_taxid_translator(query_lineages)
                            query_ranks=ncbi.get_rank(query_names.keys())
                            query_fam=str(query_names[[key  for (key, value) in query_ranks.items() if value == u'family'][0]])
                            if len(query_fam)==0:
                                res="NA"
                                with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                    file.write("family rank not found in NCBI taxonomy for query: "+query_name+'\n')
                                    file.close()
                            elif len(sub_fam)==0:
                                res="NA"
                                with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                                    file.write("family rank not found in NCBI taxonomy for subject: "+sub_name+'\n')
                                    file.close()
                            else:
                                if query_fam == sub_fam:
                                    res="TRUE"
                                else:
                                    res="FALSE"
                else:
                    res="NA"
                    with open(str(outpath+"/"+b+"_errors"+".log"),"a+") as file:
                        file.write("taxid error on subject: "+sub_name+'\n')
                        file.close()
                taxa_dict[query_name][b]=res



for t in taxa_dict:
    toprint=list()
    toprint.append(t)
    for b in list_barcodes:
        toprint.append(taxa_dict[t][b])
    print("\t".join(toprint))
