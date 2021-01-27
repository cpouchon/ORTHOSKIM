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
parser.add_argument("-i","--infile", help="samples file",
                    type=str)
parser.add_argument("-o","--outfile", help="out file path/name",
                    type=str)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()


sample_file = args.infile
out_file=args.outfile


ncbi = NCBITaxa()
with open(out_file, 'w') as outlog:
    outlog.close()

samples_dict={}
with open(sample_file) as f:
    for l in f:
        tab=l.rstrip().split('\t')
        code=tab[0]
        res = [int(i) for i in code.split("_") if i.isdigit()]
        if len(res)>0:
            code_txid=res[0]
        else:
            code_txid='NA'
        code_genus=code.split("_")[0]
        code_species=" ".join(code.split("_")[0:2])
        samples_dict[code]=dict()
        samples_dict[code]["taxid"]=code_txid
        samples_dict[code]["genus"]=code_genus
        samples_dict[code]["species"]=code_species


for taxa in samples_dict.keys():
    if samples_dict[taxa]["taxid"]=="NA":
        if len(ncbi.get_name_translator([samples_dict[taxa]["species"]]))>0:
            query_taxid=int(ncbi.get_name_translator([samples_dict[taxa]["species"]])[samples_dict[taxa]["species"]][0])
            sub_lineages=ncbi.get_lineage(query_taxid)
            sub_names=ncbi.get_taxid_translator(sub_lineages)
            sub_ranks=ncbi.get_rank(sub_names.keys())
            list_fam=[key  for (key, value) in sub_ranks.items() if value == u'family']
            list_order=[key  for (key, value) in sub_ranks.items() if value == u'order']
            list_phyllum=[key  for (key, value) in sub_ranks.items() if value == u'phylum']
            list_kingdom=[key  for (key, value) in sub_ranks.items() if value == u'kingdom']
            if len(list_fam)==0:
                fam="NA"
            else:
                fam=ncbi.get_taxid_translator(list_fam)[list_fam[0]]
            if len(list_order)==0:
                ord="NA"
            else:
                ord=ncbi.get_taxid_translator(list_order)[list_order[0]]
            if len(list_phyllum)==0:
                phyl="NA"
            else:
                phyl=ncbi.get_taxid_translator(list_phyllum)[list_phyllum[0]]
            if len(list_kingdom)==0:
                king="NA"
            else:
                king=ncbi.get_taxid_translator(list_kingdom)[list_kingdom[0]]
        else:
            if len(ncbi.get_name_translator([samples_dict[taxa]["genus"]]))>0:
                query_taxid=int(ncbi.get_name_translator([samples_dict[taxa]["genus"]])[samples_dict[taxa]["genus"]][0])
                sub_lineages=ncbi.get_lineage(query_taxid)
                sub_names=ncbi.get_taxid_translator(sub_lineages)
                sub_ranks=ncbi.get_rank(sub_names.keys())
                list_fam=[key  for (key, value) in sub_ranks.items() if value == u'family']
                list_order=[key  for (key, value) in sub_ranks.items() if value == u'order']
                list_phyllum=[key  for (key, value) in sub_ranks.items() if value == u'phylum']
                list_kingdom=[key  for (key, value) in sub_ranks.items() if value == u'kingdom']
                if len(list_fam)==0:
                    fam="NA"
                else:
                    fam=ncbi.get_taxid_translator(list_fam)[list_fam[0]]
                if len(list_order)==0:
                    ord="NA"
                else:
                    ord=ncbi.get_taxid_translator(list_order)[list_order[0]]
                if len(list_phyllum)==0:
                    phyl="NA"
                else:
                    phyl=ncbi.get_taxid_translator(list_phyllum)[list_phyllum[0]]
                if len(list_kingdom)==0:
                    king="NA"
                else:
                    king=ncbi.get_taxid_translator(list_kingdom)[list_kingdom[0]]
    else:
        if len(ncbi.get_taxid_translator([samples_dict[taxa]["taxid"]]))==0:
            if len(ncbi.get_name_translator([samples_dict[taxa]["species"]]))>0:
                query_taxid=int(ncbi.get_name_translator([samples_dict[taxa]["species"]])[samples_dict[taxa]["species"]][0])
                sub_lineages=ncbi.get_lineage(query_taxid)
                sub_names=ncbi.get_taxid_translator(sub_lineages)
                sub_ranks=ncbi.get_rank(sub_names.keys())
                list_fam=[key  for (key, value) in sub_ranks.items() if value == u'family']
                list_order=[key  for (key, value) in sub_ranks.items() if value == u'order']
                list_phyllum=[key  for (key, value) in sub_ranks.items() if value == u'phylum']
                list_kingdom=[key  for (key, value) in sub_ranks.items() if value == u'kingdom']
                if len(list_fam)==0:
                    fam="NA"
                else:
                    fam=ncbi.get_taxid_translator(list_fam)[list_fam[0]]
                if len(list_order)==0:
                    ord="NA"
                else:
                    ord=ncbi.get_taxid_translator(list_order)[list_order[0]]
                if len(list_phyllum)==0:
                    phyl="NA"
                else:
                    phyl=ncbi.get_taxid_translator(list_phyllum)[list_phyllum[0]]
                if len(list_kingdom)==0:
                    king="NA"
                else:
                    king=ncbi.get_taxid_translator(list_kingdom)[list_kingdom[0]]
            else:
                if len(ncbi.get_name_translator([samples_dict[taxa]["genus"]]))>0:
                    query_taxid=int(ncbi.get_name_translator([samples_dict[taxa]["genus"]])[samples_dict[taxa]["genus"]][0])
                    sub_lineages=ncbi.get_lineage(query_taxid)
                    sub_names=ncbi.get_taxid_translator(sub_lineages)
                    sub_ranks=ncbi.get_rank(sub_names.keys())
                    list_fam=[key  for (key, value) in sub_ranks.items() if value == u'family']
                    list_order=[key  for (key, value) in sub_ranks.items() if value == u'order']
                    list_phyllum=[key  for (key, value) in sub_ranks.items() if value == u'phylum']
                    list_kingdom=[key  for (key, value) in sub_ranks.items() if value == u'kingdom']
                    if len(list_fam)==0:
                        fam="NA"
                    else:
                        fam=ncbi.get_taxid_translator(list_fam)[list_fam[0]]
                    if len(list_order)==0:
                        ord="NA"
                    else:
                        ord=ncbi.get_taxid_translator(list_order)[list_order[0]]
                    if len(list_phyllum)==0:
                        phyl="NA"
                    else:
                        phyl=ncbi.get_taxid_translator(list_phyllum)[list_phyllum[0]]
                    if len(list_kingdom)==0:
                        king="NA"
                    else:
                        king=ncbi.get_taxid_translator(list_kingdom)[list_kingdom[0]]
                else:
                    fam="NA"
                    ord="NA"
                    phyl="NA"
                    king="NA"

        else:
            sub_lineages=ncbi.get_lineage(samples_dict[taxa]["taxid"])
            sub_names=ncbi.get_taxid_translator(sub_lineages)
            sub_ranks=ncbi.get_rank(sub_names.keys())
            list_fam=[key  for (key, value) in sub_ranks.items() if value == u'family']
            list_order=[key  for (key, value) in sub_ranks.items() if value == u'order']
            list_phyllum=[key  for (key, value) in sub_ranks.items() if value == u'phylum']
            list_kingdom=[key  for (key, value) in sub_ranks.items() if value == u'kingdom']
            if len(list_fam)==0:
                fam="NA"
            else:
                fam=ncbi.get_taxid_translator(list_fam)[list_fam[0]]
            if len(list_order)==0:
                ord="NA"
            else:
                ord=ncbi.get_taxid_translator(list_order)[list_order[0]]
            if len(list_phyllum)==0:
                phyl="NA"
            else:
                phyl=ncbi.get_taxid_translator(list_phyllum)[list_phyllum[0]]
            if len(list_kingdom)==0:
                king="NA"
            else:
                king=ncbi.get_taxid_translator(list_kingdom)[list_kingdom[0]]
    toprint=list()
    toprint.append(taxa)
    toprint.append(king)
    toprint.append(phyl)
    toprint.append(ord)
    toprint.append(fam)
    with open(out_file, 'a+') as outlog:
        outlog.write('\t'.join(toprint)+'\n')
        outlog.close()
