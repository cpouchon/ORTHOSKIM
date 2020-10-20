#!/usr/bin/env python

import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
#from Bio.Alphabet import *
from Bio.SeqRecord import *
#from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA;
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
parser.add_argument("-s","--seeds", help="sequence seeds",
                    type=str)
parser.add_argument("--target", help="gene type targeted",
                    type=str,choices=["chloroplast_CDS","chloroplast_tRNA","chloroplast_rRNA", "mitochondrion_CDS","mitochondrion_rRNA","nucrdna","nucleus_aa","nucleus_nt"])

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

ref_seeds=args.seeds
model=args.mode
outpath=args.outdir
query_name = args.query
mol=args.target
input_file = args.infile


# mol="chloroplast"
# #query_name = "Androsace_pubescens_229451_PHA000540_AXZ_B"
# outpath="./"
# typeg = "CDS"
# input_file = "chl_CDS_unaligned.fa"
# input_list_genes = "/Users/pouchonc/Desktop/Scripts/OrthoSkim/ressources/listGenes.chloro"
# ref_seeds="/Users/pouchonc/Desktop/Scripts/OrthoSkim/ressources/chloroCDS.seeds"
# model = "distance"
# # /Users/pouchonc/Desktop/Scripts/OrthoSkim/ressources/FlorAlp_Genus_matrix.csv

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
        res = [int(i) for i in query_name.split("_") if i.isdigit()]
        query=res[0]
        ncbi = NCBITaxa()
        if args.update:
            ncbi.update_taxonomy_database()
    else:
        pass

mkdir(outpath)


fname = "closed_"+str(mol)+".fa"
open(os.path.join(outpath, fname), 'w').close()

# with open(input_list_genes) as f:
#     genes = f.readlines()
# types=[]
# for g in genes:
#     g_tab = g.rstrip().split('\t')
#     types.append(g_tab[0])


seeds_seq={}
# store seeds sequences
seeds_genome = SeqIO.parse(ref_seeds, "fasta")
for record in seeds_genome:
    seqID=record.id
    sequence=record.seq
    genename=seqID.split("_")[0]
    tax_id=seqID.split("_")[1]
    species="_".join(seqID.split("_")[2:len(seqID.split("_"))])
    genus=seqID.split("_")[2:len(seqID.split("_"))][0]
    tokeep={}
    tokeep['seq']=sequence
    tokeep['id']=str(tax_id)+"_"+str(species)
    if genename not in seeds_seq.keys():
        seeds_seq[genename]=dict()
        if tax_id not in seeds_seq[genename].keys():
            seeds_seq[genename][tax_id]=[]
            seeds_seq[genename][tax_id].append(tokeep)
        else:
            seeds_seq[genename][tax_id].append(tokeep)
    else:
        if tax_id not in seeds_seq[genename].keys():
            seeds_seq[genename][tax_id]=[]
            seeds_seq[genename][tax_id].append(tokeep)
        else:
            seeds_seq[genename][tax_id].append(tokeep)


if model=="taxonomy":
    if len(ncbi.get_taxid_translator([int(query)]))==0:
        if len(ncbi.get_name_translator([" ".join(query_name.split("_")[0:2])]))==0:
            if len(ncbi.get_name_translator([" ".join(query_name.split("_")[0:1])]))==0:
                print("WARN: Error on query taxid: %s, the seed sequences will be considered as closest sequences for each gene for this taxa." % str(query_name))
                for g in seeds_seq.keys():
                    if '5.8S' in g:
                        gene_id="rrn5.8S"
                    else:
                        gene_id=g
                    out_header=">"+str(gene_id)+"_"+str(seeds_seq[gene_id][list(seeds_seq[gene_id].keys())[0]][0]['id'])
                    out_seq=str(seeds_seq[gene_id][list(seeds_seq[gene_id].keys())[0]][0]['seq'])
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
                os._exit(0)
            else:
                genus_name=" ".join(query_name.split("_")[0:1])
                query=int(ncbi.get_name_translator([genus_name])[genus_name][0])
        else:
            species_name=" ".join(query_name.split("_")[0:2])
            query=int(ncbi.get_name_translator([species_name])[species_name][0])
    else:
        pass


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

    tokeep={}
    tokeep['seq']=sequence
    tokeep['id']=str(tax_id)+"_"+str(species)

    if tax_id=="NA":
        continue
    else:

        if model=="taxonomy":
            if genename not in stored.keys():
                stored[genename]=dict()
                if tax_id not in stored[genename].keys():
                    stored[genename][tax_id]=[]
                    stored[genename][tax_id].append(tokeep)
                else:
                    stored[genename][tax_id].append(tokeep)
            else:
                if tax_id not in stored[genename].keys():
                    stored[genename][tax_id]=[]
                    stored[genename][tax_id].append(tokeep)
                else:
                    stored[genename][tax_id].append(tokeep)

        elif model=="distance":
            if genename not in stored.keys():
                stored[genename]=dict()
                if genus not in stored[genename].keys():
                    stored[genename][genus]=[]
                    stored[genename][genus].append(tokeep)
                else:
                    stored[genename][genus].append(tokeep)
            else:
                if genus not in stored[genename].keys():
                    stored[genename][genus]=[]
                    stored[genename][genus].append(tokeep)
                else:
                    stored[genename][genus].append(tokeep)

"we search the closed ref for each gene"
for g in seeds_seq.keys():
    if '5.8S' in g:
        gene_id="rrn5.8S"
    else:
        gene_id=g
    if gene_id in stored.keys():
        if model=="taxonomy":
            taxid_list=[]
            if str(query) in stored[gene_id].keys():
                closed_taxa = [str(query)]
                sub_best=list()
                sub_name=list()
                for taxa in closed_taxa:
                    sub_lmax=list()
                    for subtaxa in stored[gene_id][taxa]:
                        sub_lmax.append(len(subtaxa['seq']))
                    sub_orderindx=sorted(range(len(sub_lmax)), key=lambda k: sub_lmax[k])
                    sub_orderindx.reverse()
                    sub_best.append(stored[gene_id][taxa][sub_orderindx[0]])

                lmax=list()
                for seq in sub_best:
                    lmax.append(len(seq['seq']))
                orderindx=sorted(range(len(lmax)), key=lambda k: lmax[k])
                orderindx.reverse()

                out_header=">"+str(gene_id)+"_"+str(sub_best[orderindx[0]]['id'])
                out_seq=str(sub_best[orderindx[0]]['seq'])

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
            else:
                for i in stored[gene_id].keys():
                    if len(ncbi.get_taxid_translator([int(i)]))==0:
                        continue
                    else:
                        taxid_list.append(int(i))

                if len(taxid_list)==0:
                    out_header=">"+str(gene_id)+"_"+str(seeds_seq[gene_id][list(seeds_seq[gene_id].keys())[0]][0]['id'])
                    out_seq=str(seeds_seq[gene_id][list(seeds_seq[gene_id].keys())[0]][0]['seq'])
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
                    print("WARN: Error to extract %s locus for the taxa: %s, the seed sequence will be choosed as closed reference." % (str(gene_id),str(query_name)))
                    continue
                else:
                    taxid_all=taxid_list
                    if query in taxid_all:
                        # a sequence of the same taxid is already present, we keep the query taxid
                        closed_taxa = [str(query)]
                    else:
                        # we infer the taxonomy tree with the ref and the query taxid
                        taxid_all.append(query)
                        tree = ncbi.get_topology(taxid_all)
                        # we check if the query taxid was translated or not
                        if len(tree.search_nodes(name=str(query)))==0:
                            if len(ncbi.get_lineage(query))>0:
                                search_taxid=ncbi.get_lineage(query)[len(ncbi.get_lineage(query))-1]
                                if len(tree.search_nodes(name=str(search_taxid)))>0:
                                    t2=tree.search_nodes(name=str(search_taxid))[0].up
                                    closed_taxa = []
                                    for l in t2.get_leaves():
                                        if l.__dict__['taxid'] in taxid_all:
                                            closed_taxa.append(str(l.__dict__['taxid']))
                                        else:
                                            continue
                                    if str(query) in closed_taxa:
                                        closed_taxa.remove(str(query))
                                    elif str(search_taxid) in closed_taxa:
                                        closed_taxa.remove(str(search_taxid))
                                    else:
                                        pass
                                else:
                                    out_header=">"+str(gene_id)+"_"+str(seeds_seq[gene_id][list(seeds_seq[gene_id].keys())[0]][0]['id'])
                                    out_seq=str(seeds_seq[gene_id][list(seeds_seq[gene_id].keys())[0]][0]['seq'])
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
                                    print("WARN: Error to extract %s locus for the taxa: %s, the seed sequence will be choosed as closed reference." % (str(gene_id),str(query_name)))
                                    continue
                            else:
                                out_header=">"+str(gene_id)+"_"+str(seeds_seq[gene_id][list(seeds_seq[gene_id].keys())[0]][0]['id'])
                                out_seq=str(seeds_seq[gene_id][list(seeds_seq[gene_id].keys())[0]][0]['seq'])
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
                                print("WARN: Error to extract %s locus for the taxa: %s, the seed sequence will be choosed as closed reference." % (str(gene_id),str(query_name)))
                                continue
                        else:
                            t2=tree.search_nodes(name=str(query))[0].up
                            closed_taxa = []
                            for l in t2.get_leaves():
                                if l.__dict__['taxid'] in taxid_all:
                                    closed_taxa.append(str(l.__dict__['taxid']))
                                else:
                                    continue
                            if str(query) in closed_taxa:
                                closed_taxa.remove(str(query))
                            else:
                                pass

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
            if len(closed_list)==0:
                lengthlist=list()
                namelist=list()
                seqlist=list()
                seqlist=[b['seq'] for c in list(stored[gene_id].values()) for b in c]
                lengthlist=[len(b['seq']) for c in list(stored[gene_id].values()) for b in c]
                namelist=[a for a in list(stored[gene_id].keys())]
                suborderindx=sorted(range(len(lengthlist)), key=lambda k: lengthlist[k])
                suborderindx.reverse()
                closed_taxa=[namelist[suborderindx[0]]]
            else:
                closed_taxa=closed_list[0]

        sub_best=list()
        sub_name=list()

        if len(closed_taxa)==0:
            out_header=">"+str(gene_id)+"_"+str(seeds_seq[gene_id][list(seeds_seq[gene_id].keys())[0]][0]['id'])
            out_seq=str(seeds_seq[gene_id][list(seeds_seq[gene_id].keys())[0]][0]['seq'])
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
            print("WARN: Error to extract %s locus for the taxa: %s, the seed sequence will be choosed as closed reference." % (str(gene_id),str(query_name)))
            continue
        else:
            for taxa in closed_taxa:
                sub_lmax=list()
                for subtaxa in stored[gene_id][taxa]:
                    sub_lmax.append(len(subtaxa['seq']))
                sub_orderindx=sorted(range(len(sub_lmax)), key=lambda k: sub_lmax[k])
                sub_orderindx.reverse()
                sub_best.append(stored[gene_id][taxa][sub_orderindx[0]])
                #sub_name.append(taxa)

            lmax=list()
            for seq in sub_best:
                lmax.append(len(seq['seq']))
            orderindx=sorted(range(len(lmax)), key=lambda k: lmax[k])
            orderindx.reverse()

            out_header=">"+str(gene_id)+"_"+str(sub_best[orderindx[0]]['id'])
            out_seq=str(sub_best[orderindx[0]]['seq'])

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
    else:
        out_header=">"+str(gene_id)+"_"+str(seeds_seq[gene_id][list(seeds_seq[gene_id].keys())[0]][0]['id'])
        out_seq=str(seeds_seq[gene_id][list(seeds_seq[gene_id].keys())[0]][0]['seq'])
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
        print("WARN: Error %s locus not found in %s dataset, the seed sequence will be choosed as closed reference." % (str(gene_id),str(input_file)))
        continue
