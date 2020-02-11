#! /usr/bin/env python

import glob
import sys
import os, errno
import argparse
import fnmatch
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

parser = argparse.ArgumentParser(description='Get statistics for matK and rbcL genes')
parser.add_argument("--samplep", help="searching path of samples assemblies",
                    type=str)
parser.add_argument("--genesp", help="searching path of genes files",
                    type=str)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

samppath = args.samplep
samples = []
for r, d, f in os.walk(samppath):
    for file in f:
        if '.fa' in file:
            samples.append(file.replace(".fa",""))

genepath = args.genesp
genesfile = []
for r, d, f in os.walk(genepath):
    for file in f:
        if '.fa' in file:
            genesfile.append(os.path.join(r, file))

stored={}
for file in genesfile:
    cf = file.rstrip()
    cfile_path = cf
    cfile_name = os.path.basename(cfile_path)
    cindex_of_dot = cfile_name.index('.')
    gene_name = cfile_name[:cindex_of_dot]
    cur_genome = SeqIO.parse(file, "fasta")
    for record in cur_genome:
            seqID=record.id.replace(";","")
            if gene_name not in stored.keys():
                stored[gene_name]=dict()
                stored[gene_name]=[]
                stored[gene_name].append(seqID)
            else:
                stored[gene_name].append(seqID)

print("taxon\ttaxid\tspecimen\tlibrary\tchloro_has_matK\tchloro_has_rbcL")
for samp in samples:
    seqID=samp
    num_split=[]
    for indx in range(len(seqID.split("_"))):
        if any(ch.isdigit() for ch in seqID.split("_")[indx]):
            num_split.append(indx)
    species_name=" ".join(seqID.split("_")[0:num_split[0]])
    taxid=seqID.split("_")[num_split[0]]
    code="_".join(seqID.split("_")[num_split[0]+1:num_split[1]+1])
    sequencing=":".join(seqID.split("_")[num_split[1]+1:len(seqID.split("_"))])
    matk=seqID in stored['matK']
    rbcl=seqID in stored['rbcL']
    print ('%s\t%s\t%s\t%s\t%s\t%s' % (str(species_name),str(taxid),str(code),str(sequencing),str(matk).upper(),str(rbcl).upper()))
