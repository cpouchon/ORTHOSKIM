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

import pandas as pd
import numpy as np
import random
import argparse
import sys
import os, errno



parser = argparse.ArgumentParser(description='Extract closed genes for a query taxa. Script was writen by C. Pouchon (2020).')
parser.add_argument("-in","--infile", help="input reference DB (fasta format)",
                    type=str)
parser.add_argument("-g","--geneslist", help="list of genes to extract",
                    type=str)
parser.add_argument("-q","--query", help="taxa name or taxid if --taxid mode",
                    type=str)
parser.add_argument("-o","--outdir", help="out directory path",
                    type=str)
parser.add_argument("-m","--mode", help="extraction mode",
                    type=str,choices=["taxonomy", "distance"])
parser.add_argument("--taxid", help="query is given by taxid (else it from it's name)",
                    action="store_true")
parser.add_argument("--update", help="update NCBI taxonomy",
                    action="store_true")
parser.add_argument("-d","--distance", help="phylogenetic distance matrix file of genus",
                    type=str)
parser.add_argument("--gene_type", help="type of genes to extract",
                    type=str,choices=["CDS", "rRNA","tRNA"])
parser.add_argument("--compartment", help="cellular compartment",
                    type=str,choices=["chloroplast", "mitochondrion","nucrdna","nucleus"])

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()


def mkdir(path, overwrite=False):
    '''
    function to create a directory for output fasta
    '''
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            if not overwrite:
                print ("path '%s' already exists" % path)   # overwrite == False and we've hit a directory that exists
        else: raise

class ProgressBar:
	'''
	Progress bar
	'''
	def __init__ (self, valmax, maxbar, title):
		if valmax == 0:  valmax = 1
		if maxbar > 200: maxbar = 200
		self.valmax = valmax
		self.maxbar = maxbar
		self.title  = title
	def update(self, val):
		import sys
		perc  = round((float(val) / float(self.valmax)) * 100)
		scale = 100.0 / float(self.maxbar)
		bar   = int(perc / scale)
		out = '\r %20s [%s%s] %3d %% ' % (self.title, '.' * bar, ' ' * (self.maxbar - bar), perc)
		sys.stdout.write(out)
		sys.stdout.flush()


model=args.mode
typeg = args.gene_type
outpath=args.outdir
query_name = args.query
mol=args.compartment
#query_name = "Androsace_pubescens_229451_PHA000540_AXZ_B"

if model=="taxonomy":
    res = [int(i) for i in query_name.split("_") if i.isdigit()]
    query=res[0]
elif model=="distance":
    res=query_name.split("_")[0]
    if res[0].isupper():
        query=res
    else:
        query=res.capitalize()

if model=="taxonomy":
    ncbi = NCBITaxa()
    if args.update:
        ncbi.update_taxonomy_database()
else:
    df = pd.read_csv(args.distance,sep=",", index_col=[0])
    genera = [col for col  in df.columns]
    if query not in genera:
        model="taxonomy"
        ncbi = NCBITaxa()
        if args.update:
            ncbi.update_taxonomy_database()
    else:
        pass



input_file = args.infile
#input_file = "/Users/pouchonc/PhyloAlps/mitoDB/refrences/mit_CDS_unaligned.fa"

input_list_genes = args.geneslist
#input_list_genes = "/Users/pouchonc/Desktop/Scripts/OrthoSkim/ressources/listGenes.mito"
with open(input_list_genes) as f:
    genes = f.readlines()
types=[]
for g in genes:
    g_tab = g.rstrip().split('\t')
    types.append(g_tab[0])

mkdir(outpath)
#model = "distance"

stored={}

'we store all sequences'
cur_genome = SeqIO.parse(input_file, "fasta")
for record in cur_genome:
    seqID=record.id
    sequence=record.seq
    genename=seqID.split("_")[0]

    tax_id=seqID.split("_")[1]
    species="_".join(seqID.split("_")[2:len(seqID.split("_"))])
    genus=seqID.split("_")[2:len(seqID.split("_"))][0]

    if model=="taxonomy":
        if genename not in stored.keys():
            stored[genename]=dict()
            if tax_id not in stored[genename].keys():
                stored[genename][tax_id]=[]
                stored[genename][tax_id].append(sequence)
            else:
                stored[genename][tax_id].append(sequence)
        else:
            if tax_id not in stored[genename].keys():
                stored[genename][tax_id]=[]
                stored[genename][tax_id].append(sequence)
            else:
                stored[genename][tax_id].append(sequence)

    elif model=="distance":
        if genename not in stored.keys():
            stored[genename]=dict()
            if genus not in stored[genename].keys():
                stored[genename][genus]=[]
                stored[genename][genus].append(sequence)
            else:
                stored[genename][genus].append(sequence)
        else:
            if genus not in stored[genename].keys():
                stored[genename][genus]=[]
                stored[genename][genus].append(sequence)
            else:
                stored[genename][genus].append(sequence)

#typeg = "CDS"
fname = "closed_"+str(mol)+"_"+str(typeg)+".fa"
open(os.path.join(outpath, fname), 'w').close()

"we search the closed ref for each gene"
for g in genes:
    g_tab = g.rstrip().split('\t')
    gene_id=g_tab[1]
    gene_type=g_tab[0]

    if gene_type == typeg:
        if model=="taxonomy":
            taxid_list = [int(i) for i in stored[gene_id].keys()]
            taxid_all=taxid_list
            if query in taxid_all:
                # a sequence of the same taxid is already present, we keep the query taxid
                closed_taxa = [query]
            else:
                # we infer the taxonomy tree with the ref and the query taxid
                taxid_all.append(query)
                tree = ncbi.get_topology(taxid_all)
                t2=tree.search_nodes(name=str(query))[0].up
                closed_taxa = [str(l.__dict__['taxid']) for l in t2.get_leaves()]
                closed_taxa.remove(str(query))

        else:
            subdf = df.sort_values(by=query)[query]
            genera_list = [str(i) for i in stored[gene_id].keys()]
            part = 0
            range_list = list()
            closed_list = list()
            for taxa in range(len(subdf.index)):
                if subdf.index[taxa] in genera_list:
                    # we create list for taxa that are at the same phylogenetic distance from the query
                    score=subdf[taxa]
                    if len(range_list)>0:
                        if score in range_list[part]:
                            closed_list[part].extend([subdf.index[taxa]])
                        else:
                            range_list.append([score])
                            closed_list.append([subdf.index[taxa]])
                            part = part+1
                    else:
                        range_list.append([score])
                        closed_list.append([subdf.index[taxa]])
            closed_taxa=closed_list[0]

        sub_best=list()
        sub_name=list()
        for taxa in closed_taxa:
            sub_lmax=list()
            for subtaxa in stored[gene_id][taxa]:
                sub_lmax.append(len(subtaxa))
            sub_orderindx=sorted(range(len(sub_lmax)), key=lambda k: sub_lmax[k])
            sub_orderindx.reverse()
            sub_best.append(stored[gene_id][taxa][sub_orderindx[0]])
            sub_name.append(taxa)

        lmax=list()
        for seq in sub_best:
            lmax.append(len(seq))
        orderindx=sorted(range(len(lmax)), key=lambda k: lmax[k])
        orderindx.reverse()

        out_header=">"+str(gene_id)+"_"+str(sub_name[orderindx[0]])
        out_seq=str(sub_best[orderindx[0]])

        if os.path.isfile(os.path.join(outpath, fname)):
            with open(os.path.join(outpath, fname), 'a+') as file:
                old_headers = []
                end_file=file.tell()
                file.seek(0)
                for line in file:
                    if line.startswith(">"):
                        old_headers.append(line.rstrip())
                if not out_header in old_headers:
                    file.seek(end_file)
                    file.write(out_header+'\n')
                    file.write(str(out_seq)+'\n')
                else:
                    pass
        else :
            with open(os.path.join(outpath, fname), 'a') as out:
                out.write(out_header+'\n')
                out.write(str(out_seq)+'\n')

    # function to now extract the ref from the closed taxa list (either taxid or genus name)
