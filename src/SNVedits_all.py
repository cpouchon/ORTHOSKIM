#! /usr/bin/env python

import glob
import sys
import os, errno
import argparse
import fnmatch
#import numpy
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqFeature import FeatureLocation


parser = argparse.ArgumentParser(description='Extract CU/UC RNA editing sites for given taxa pairs')
parser.add_argument("-p","--path", help="searching path of aligned sequences files",
                    type=str)
parser.add_argument("-t","--taxa", help="taxa pairs to compare",
                    type=str)
parser.add_argument("-trim", help="[mode] statistics made on trimmed alignments",
                    action="store_true")
parser.add_argument("-o","--out", help="taxa number in alignment",
                    type=str)
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}


if args.trim:
    'Search all fasta files files'
    path = args.path
    in_files = []
    for r, d, f in os.walk(path):
        for file in f:
            if '.trim' in file:
                in_files.append(os.path.join(r, file))
else:
        'Search all fasta files files'
        path = args.path
        in_files = []
        for r, d, f in os.walk(path):
            for file in f:
                if '.fa' in file:
                    in_files.append(os.path.join(r, file))


outpath=args.out

giventaxa = args.taxa
with open(giventaxa) as f:
    taxal = f.readlines()

taxa_dict={}
taxa_seq={}
taxa_list = list()
for line in taxal:
    l = line.rstrip().split("\t")
    if len(l)>0:
        taxa_dict.setdefault(l[0], {})
        taxa_seq.setdefault(l[0], {})

for file in in_files:
    cf = file.rstrip()
    cfile_path = cf
    cfile_name = os.path.basename(cfile_path)
    cindex_of_dot = cfile_name.index('.')
    gene_name = cfile_name[:cindex_of_dot]
    gene_type = cf.split("/")[len(cf.split("/"))-3]

    cur_genome=to_dict_remove_dups(SeqIO.parse(file, "fasta"))

    for line in taxal:
        l = line.rstrip().split("\t")
        geno=l[0]
        rna=l[1]
        if geno in list(cur_genome.keys()) and rna in list(cur_genome.keys()):
            base_count=0
            at=0
            ac=0
            ag=0
            ta=0
            ca=0
            ga=0
            tc=0
            tg=0
            ct=0
            gt=0
            cg=0
            gc=0
            taxa_dict[geno][gene_name]={}
            for b in range(len(cur_genome[geno].seq)):
                geno_b=cur_genome[geno].seq[b]
                rna_b=cur_genome[rna].seq[b]
                if geno_b =="-" and rna_b=="-":
                    continue
                elif geno_b =="-" and rna_b=="n":
                    continue
                elif geno_b =="n" and rna_b=="n":
                    continue
                elif geno_b =="n" and rna_b=="-":
                    continue
                else:
                    base_count=base_count+1
                    if geno_b==rna_b:
                        continue
                    else:
                        if geno_b=="-" or rna_b=="-":
                            continue
                        else:
                            #count
                            if geno_b=="a" and rna_b=="t":
                                 at=at+1
                            elif geno_b=="a" and rna_b=="c":
                                 ac=ac+1
                            elif geno_b=="a" and rna_b=="g":
                                 ag=ag+1
                            elif geno_b=="t" and rna_b=="a":
                                 ta=ta+1
                            elif geno_b=="c" and rna_b=="a":
                                 ca=ca+1
                            elif geno_b=="g" and rna_b=="a":
                                 ga=ga+1
                            elif geno_b=="t" and rna_b=="c":
                                 tc=tc+1
                                 if b>=4 and b<=(len(cur_genome[geno].seq)-4):
                                     # cond no variants
                                     cond_seq=0
                                     for w in range(b-3,b+3):
                                         if str(cur_genome[geno].seq[w]) != str(cur_genome[rna].seq[w]):
                                             cond_seq=cond_seq+1
                                     if cond_seq <2:
                                         if len(taxa_seq[geno].keys())==0:
                                             taxa_seq[geno]["UC"]={}
                                             taxa_seq[geno]["UC"]["edit"]=[]
                                             taxa_seq[geno]["UC"]["unedit"]=[]
                                             taxa_seq[geno]["CU"]={}
                                             taxa_seq[geno]["CU"]["edit"]=[]
                                             taxa_seq[geno]["CU"]["unedit"]=[]
                                         else:
                                             pass
                                         taxa_seq[geno]["UC"]["unedit"].append(str(cur_genome[geno].seq[b-3:b+4]))
                                         taxa_seq[geno]["UC"]["edit"].append(str(cur_genome[rna].seq[b-3:b+4]))
                            # elif geno_b=="t" and rna_b=="g":
                            #     tg=tg+1
                            elif geno_b=="c" and rna_b=="t":
                                 ct=ct+1
                                 if b>=4 and b<=(len(cur_genome[geno].seq)-4):
                                     # cond no variants
                                     cond_seq=0
                                     for w in range(b-3,b+3):
                                         if str(cur_genome[geno].seq[w]) != str(cur_genome[rna].seq[w]):
                                             cond_seq=cond_seq+1
                                     if cond_seq <2:
                                         if len(taxa_seq[geno].keys())==0:
                                             taxa_seq[geno]["UC"]={}
                                             taxa_seq[geno]["UC"]["edit"]=[]
                                             taxa_seq[geno]["UC"]["unedit"]=[]
                                             taxa_seq[geno]["CU"]={}
                                             taxa_seq[geno]["CU"]["edit"]=[]
                                             taxa_seq[geno]["CU"]["unedit"]=[]
                                         else:
                                             pass
                                         taxa_seq[geno]["CU"]["unedit"].append(str(cur_genome[geno].seq[b-3:b+4]))
                                         taxa_seq[geno]["CU"]["edit"].append(str(cur_genome[rna].seq[b-3:b+4]))
                            elif geno_b=="g" and rna_b=="t":
                                 gt=gt+1
                            elif geno_b=="c" and rna_b=="g":
                                 cg=cg+1
                            elif geno_b=="g" and rna_b=="c":
                                 gc=gc+1

            taxa_dict[geno][gene_name]["count"]=base_count
            taxa_dict[geno][gene_name]["at"]=at
            taxa_dict[geno][gene_name]["ac"]=ac
            taxa_dict[geno][gene_name]["ag"]=ag
            taxa_dict[geno][gene_name]["ta"]=ta
            taxa_dict[geno][gene_name]["ca"]=ca
            taxa_dict[geno][gene_name]["ga"]=ga
            taxa_dict[geno][gene_name]["tc"]=tc
            taxa_dict[geno][gene_name]["tg"]=tg
            taxa_dict[geno][gene_name]["ct"]=ct
            taxa_dict[geno][gene_name]["gt"]=gt
            taxa_dict[geno][gene_name]["cg"]=cg
            taxa_dict[geno][gene_name]["gc"]=gc
            taxa_dict[geno][gene_name]["type"]=gene_type

# output les frequences de CU et UC RNA editing par gene (edits/1000nt)

for line in taxal:
    l = line.rstrip().split("\t")
    geno=l[0]
    for g in list(taxa_dict[geno].keys()):
        fcu=round(float(taxa_dict[geno][g]["ct"])*1000/float(taxa_dict[geno][g]["count"]),2)
        fuc=round(float(taxa_dict[geno][g]["tc"])*1000/float(taxa_dict[geno][g]["count"]),2)
        toprint=list()
        toprint.append(geno)
        toprint.append(g)
        toprint.append(taxa_dict[geno][g]["type"])
        toprint.append(str(fcu))
        toprint.append(str(fuc))
        with open(os.path.join(outpath, "frequencies_RNAediting_CU-UC.log"),"a+") as file:
            file.write("\t".join(toprint)+"\n")
        toprint=list()
        toprint.append(geno)
        toprint.append(g)
        toprint.append(taxa_dict[geno][g]["type"])
        toprint.append(str(taxa_dict[geno][g]["count"]))
        toprint.append(str(taxa_dict[geno][g]["at"]))
        toprint.append(str(taxa_dict[geno][g]["ac"]))
        toprint.append(str(taxa_dict[geno][g]["ag"]))
        toprint.append(str(taxa_dict[geno][g]["ta"]))
        toprint.append(str(taxa_dict[geno][g]["ca"]))
        toprint.append(str(taxa_dict[geno][g]["ga"]))
        toprint.append(str(taxa_dict[geno][g]["tc"]))
        toprint.append(str(taxa_dict[geno][g]["tg"]))
        toprint.append(str(taxa_dict[geno][g]["ct"]))
        toprint.append(str(taxa_dict[geno][g]["gt"]))
        toprint.append(str(taxa_dict[geno][g]["cg"]))
        toprint.append(str(taxa_dict[geno][g]["gc"]))    
        with open(os.path.join(outpath, "count_RNAediting_all.log"),"a+") as file:
            file.write("\t".join(toprint)+"\n")

    if "CU" in list(taxa_seq[geno].keys()):
        for s in taxa_seq[geno]["CU"]["edit"]:
            with open(os.path.join(outpath, str(geno+"_CU.log")),"a+") as file:
                file.write(s+"\n")
    if "UC" in list(taxa_seq[geno].keys()):
        for s in taxa_seq[geno]["UC"]["edit"]:
            with open(os.path.join(outpath, str(geno+"_UC.log")),"a+") as file:
                file.write(s+"\n")
