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



parser = argparse.ArgumentParser(description='Capture of the besthits from blast table output according to the pident. Script was writen by C. Pouchon (2020).')
parser.add_argument("-o","--outdir", help="out directory path",
                    type=str)
parser.add_argument('-b', '--barcodes', help='list of barcodes', type=str)
parser.add_argument("-t","--taxa", help="list of taxa to include in the final concatenated file",
                    type=str)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

outpath=args.outdir
giventaxa = args.taxa
# list_accessions = ""
# list_barcodes=['matK']
# outpath="./"
# giventaxa = "testdata_taxa"

ncbi = NCBITaxa()

list_barcodes = [str(item) for item in args.barcodes.split(',')]

taxa_dict={}
with open(giventaxa) as f:
    for l in f:
        taxa_dict.setdefault(l.rstrip(), {})
        for b in list_barcodes:
            taxa_dict[l.rstrip()][b]="NA"

for b in list_barcodes:
    bestscore={}
    besthits={}
    with open(str(outpath+"/"+"blast_"+b+".out")) as out:
        for l in out:
            sub_error=0
            tab=l.rstrip().split('\t')
            query_name=tab[0].replace(";","")
            sub_name=tab[1]
            score=tab[2]

            if query_name in taxa_dict.keys():
                if query_name in besthits.keys():
                    if int(score) > bestscore[query_name]:
                        besthits[query_name]=[sub_name]
                        bestscore[query_name]=int(score)
                    elif int(score) == bestscore[query_name]:
                        if sub_name not in besthits[query_name]:
                            besthits[query_name].append(sub_name)
                    elif int(score) > bestscore[query_name]:
                        pass
                else:
                    besthits[query_name]=[sub_name]
                    bestscore[query_name]=int(score)
            else:
                with open(str(outpath+"/"+b+"_besthits_errors"+".log"),"a+") as file:
                    file.write(query_name+' not found in '+giventaxa+'\n')
                    file.close()

    for taxa in besthits:
        query_name=taxa
        sub_name=besthits[taxa][0]
        with open(str(outpath+"/"+b+"_besthits"+".out"),"a+") as file:
            file.write(query_name+'\t'+sub_name+'\n')
            file.close()
