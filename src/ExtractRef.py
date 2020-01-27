#!/usr/bin/env python

import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.Alphabet import *
from Bio.SeqRecord import *
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA;
from Bio.Seq import *
from Bio.SeqUtils import *
from ete2 import NCBITaxa

import pandas as pd
import numpy as np
import random
import argparse
import sys
import os, errno



parser = argparse.ArgumentParser(description='Extract genes for annotation files. Script was writen by C. Pouchon (2019).')
parser.add_argument("-in","--infile", help="input reference DB (fasta format)",
                    type=str)
parser.add_argument("-g","--geneslist", help="list of genes to extract",
                    type=str)
parser.add_argument("-q","--query", help="taxa name or taxid if --taxid mode",
                    type=str)
parser.add_argument("-o","--outdir", help="out directory path",
                    type=str)
parser.add_argument("-m","--mode", help="extraction mode",
                    type=str,choices=["taxonomy", "tree"])
parser.add_argument("--taxid", help="query is given by taxid (else it from it's name)",
                    action="store_true")
parser.add_argument("--update", help="update NCBI taxonomy",
                    action="store_true")
parser.add_argument("-t","--tree", help="phylogenetic tree of genus",
                    type=str)
parser.add_argument("-tfmt","--treefmt", help="format of tree",
                    type=str,choices=["newick", "nexus"])

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


query_name = "Androsace_pubescens_229451_PHA000540_AXZ_B"
if model=="taxonomy":
    res = [int(i) for i in query_name.split("_") if i.isdigit()]
    query=res[0]
elif model=="tree":
    res=query_name.split("_")[0]
    if res[0].isupper():
        query=res
    else:
        query=res.capitalize()


input_tree="/Users/pouchonc/PhyloAlps/FlorAlp_Genus_tree.tre"
treefmt="newick"

if model=="taxonomy":
    ncbi = NCBITaxa()
    if args.update:
        ncbi.update_taxonomy_database()
else:
    input_tree = args.tree
    treefmt = args.treefmt

    tree = Tree(input_tree,format=1)
    taxa=tree.get_leaf_names()

    if query in taxa:
        print 'tutu'
    else:
        model="taxonomy"
        ncbi = NCBITaxa()
        if args.update:
            ncbi.update_taxonomy_database()

    tree = dendropy.Tree.get(path=input_tree,schema=treefmt)
    pdm = tree.phylogenetic_distance_matrix()
    taxa = [i.label for i in tree.taxon_set]
    colname = taxa
    rowname = taxa
    matrix = np.reshape(range(len(taxa)*len(taxa)),(len(taxa),len(taxa)))
    df = pd.DataFrame(matrix,columns=colname,index=rowname)

    for t1 in taxa:
        for t2 in taxa:
            if t1 == t2:
                df.loc[t1,t2] = 0.0
            else:
                df.loc[t1,t2] = tree.get_distance(t1,t2)

# ENFAIT IL FAUT FAIRE LA MATRICE DE DISTANCE AVANT

input_file = args.infile
input_file = "mit_CDS_unaligned.fa"

input_list_genes = args.geneslist
with open(input_list_genes) as f:
    genes = f.readlines()
types=[]
for g in genes:
    g_tab = g.rstrip().split('\t')
    types.append(g_tab[0])



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

    elif model=="tree":
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



"we search the closed ref for each gene"
for g in genes:
    g_tab = g.rstrip().split('\t')
    gene_id=g_tab[1]
    gene_type=g_tab[0]


    if model=="taxonomy":
        taxid_list = [int(i) for i in stored[gene_id].keys()]
        taxid_all=taxid_list
        taxid_all.append(query)

        tree = ncbi.get_topology(taxid_all)

        t2=tree.search_nodes(name=str(query))[0].up
        closed_taxa = [str(l.__dict__['taxid']) for l in t2.get_leaves()]
        closed_taxa.remove(str(query))

    else:

        subdf = df.sort_values(by=ingenus)[ingenus]
for taxa in range(len(subdf.index)):
    if subdf.index[taxa] in stored.keys():
        closedtaxa_file=stored[subdf.index[taxa]][0]
        break
