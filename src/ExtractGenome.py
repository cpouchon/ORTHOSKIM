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


parser = argparse.ArgumentParser(description='Extract closed genomes for a query taxa to assign contigs in genome type. Script was writen by C. Pouchon (2020).')
parser.add_argument("-in","--infile", help="input genomes annotation file",
                    type=str)
parser.add_argument("-fmt","--annotationfmt", help="format of annotated file",choices=["embl", "genbank"],
                    type=str)
parser.add_argument("-q","--query", help="taxa name or taxid if --taxid mode",
                    type=str)
parser.add_argument("-o","--outdir", help="out directory path",
                    type=str)
parser.add_argument("-m","--mode", help="extraction mode",
                    type=str,choices=["taxonomy", "distance"])
parser.add_argument("-d","--distance", help="phylogenetic distance matrix file of genus",
                    type=str)
parser.add_argument("-s","--size", help="minimum size of genome",
                    type=int)
parser.add_argument("--compartment", help="cellular compartment",
                    type=str,choices=["chloroplast", "mitochondrion","nucrdna"])


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


model=args.mode
outpath=args.outdir
query_name = args.query
mol=args.compartment
formatin = args.annotationfmt
cond_size = args.size

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
else:
    df = pd.read_csv(args.distance,sep=",", index_col=[0])
    genera = [col for col  in df.columns]
    if query not in genera:
        model="taxonomy"
        res = [int(i) for i in query_name.split("_") if i.isdigit()]
        query=res[0]
        ncbi = NCBITaxa()
    else:
        pass


mkdir(outpath)

fname = "closed_"+str(mol)+".fa"
open(os.path.join(outpath, fname), 'w').close()


stored={}

f=args.infile
file_path = f
file_name = os.path.basename(file_path)
cur_genome = SeqIO.parse(f, formatin)
for record in cur_genome:
    break_record=0
    for feat in record.features:
        if feat.type == 'source':
            if "/" in feat.qualifiers['organism'][0]:
                break_record=break_record+1 #verify this process
            else:
                org=feat.qualifiers['organism'][0].split(" ")
                for parts in feat.qualifiers['db_xref']:
                    if "taxon" in parts:
                        taxid=parts.split(":")[1]
                taxaname=str(taxid+"_"+"_".join(org))
        else:
            break
    if break_record==0:
        sequence=str(record.seq)
        genus=org[0]
        tokeep={}
        tokeep[taxaname]=sequence

        if taxid=="NA":
            continue
        else:
            if len(sequence) >= cond_size:
                if model=="taxonomy":
                    if taxid not in stored.keys():
                        stored[taxid]={}
                        stored[taxid][taxaname]=sequence
                    else:
                        if taxaname in stored[taxid].keys():
                            pass
                        else:
                            stored[taxid][taxaname]=sequence

                elif model=="distance":
                    if genus not in stored.keys():
                        stored[genus]={}
                        stored[genus][taxaname]=sequence
                    else:
                        if taxaname in stored[genus].keys():
                            pass
                        else:
                            stored[genus][taxaname]=sequence
            else:
                continue

    else:
        continue


if model=="taxonomy":
    if len(ncbi.get_taxid_translator([int(query)]))==0:
        if len(ncbi.get_name_translator([" ".join(query_name.split("_")[0:2])]))==0:
            if len(ncbi.get_name_translator([" ".join(query_name.split("_")[0:1])]))==0:
                closed_taxa = random.sample(stored.keys(), k=5)
                cond_len = 1
                sub_best={}
                for taxa in closed_taxa:
                    if cond_len == 6:
                        break
                    else:
                        for subtaxa in stored[taxa]:
                            if cond_len == 6:
                                break
                            else:
                                sub_best[subtaxa]=stored[taxa][subtaxa]
                                cond_len += 1
                for seq in sub_best:
                    out_header=">"+str(seq)
                    out_seq=str(sub_best[seq])
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
                print("WARN: Error on query taxid: %s, five random genomes will be considered." % str(query_name))
                os._exit(0)
            else:
                genus_name=" ".join(query_name.split("_")[0:1])
                query=int(ncbi.get_name_translator([genus_name])[genus_name][0])
        else:
            species_name=" ".join(query_name.split("_")[0:2])
            query=int(ncbi.get_name_translator([species_name])[species_name][0])
    else:
        pass

if model=="taxonomy":
    taxid_list=[]
    for i in stored.keys():
        if len(ncbi.get_taxid_translator([int(i)]))==0:
            continue
        else:
            taxid_list.append(int(i))

    if len(taxid_list)==0:
        closed_taxa = random.sample(stored.keys(), k=5)
    else:
        taxid_all=taxid_list
        if query in taxid_all:
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
                        print("Error on taxid in NCBI for query taxa: %s; five random genomes will be considered." % str(query_name))
                        closed_taxa = random.sample(stored.keys(), k=5)
                else:
                    print("Error on taxid in NCBI for query taxa: %s; five random genomes will be considered." % str(query_name))
                    closed_taxa = random.sample(stored.keys(), k=5)
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
                        print("Error on taxid in NCBI for query taxa: %s; five random genomes will be considered." % str(query_name))
                        closed_taxa = random.sample(stored.keys(), k=5)
                else:
                    print("Error on taxid in NCBI for query taxa: %s; five random genomes will be considered." % str(query_name))
                    closed_taxa = random.sample(stored.keys(), k=5)
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
    genera_list = [str(i) for i in stored.keys()]
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

    cond_len = 1
    closed_taxa=[]
    for i in closed_list:
        if cond_len == 6:
            break
        else:
            for j in i:
                if cond_len == 6:
                    break
                else:
                    closed_taxa.append(j)
                    cond_len += 1


if len(closed_taxa)==0:
    print("Error on selection of closed genomes for query taxa: %s; five random genomes will be considered." % str(query_name))
    closed_taxa = random.sample(stored.keys(), k=5)
    cond_len = 1
    sub_best={}
    for taxa in closed_taxa:
        if cond_len == 6:
            break
        else:
            for subtaxa in stored[taxa]:
                if cond_len == 6:
                    break
                else:
                    sub_best[subtaxa]=stored[taxa][subtaxa]
                    cond_len += 1

    for seq in sub_best:
        out_header=">"+str(seq)
        out_seq=str(sub_best[seq])
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
    if len(closed_taxa)>=5:
        cond_len = 1
        sub_best={}
        for taxa in closed_taxa:
            if cond_len == 6:
                break
            else:
                for subtaxa in stored[taxa]:
                    if cond_len == 6:
                        break
                    else:
                        sub_best[subtaxa]=stored[taxa][subtaxa]
                        cond_len += 1

        for seq in sub_best:
            out_header=">"+str(seq)
            out_seq=str(sub_best[seq])
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
        n_to_add=5-len(closed_taxa)
        closed_taxa.extend(random.sample(stored.keys(), k=n_to_add))
        cond_len = 1
        sub_best={}
        for taxa in closed_taxa:
            if cond_len == 6:
                break
            else:
                for subtaxa in stored[taxa]:
                    if cond_len == 6:
                        break
                    else:
                        sub_best[subtaxa]=stored[taxa][subtaxa]
                        cond_len += 1

        for seq in sub_best:
            out_header=">"+str(seq)
            out_seq=str(sub_best[seq])
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
