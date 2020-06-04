#! /usr/bin/env python

import glob
import sys
import os, errno
import argparse
import fnmatch
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from joblib import Parallel, delayed
import multiprocessing


parser = argparse.ArgumentParser(description='Selection of taxa from fasta files')
parser.add_argument("--inpath", help="searching path of sequences files",
                    type=str)
parser.add_argument("--outpath", help="output path to write sequences",
                    type=str)
parser.add_argument("-t","--taxa", help="list of taxa to include in the final concatenated file",
                    type=str)
parser.add_argument("-e","--extension", help="extension of files to search (i.e fa,fna,fasta or custom)",
                    type=str)
parser.add_argument("--threads", help="number of threads to use",
                    type=int)
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


file_extension=args.extension
'Search all fasta files files'
path = args.inpath
in_files = []
for r, d, f in os.walk(path):
    for file in f:
        if str("."+file_extension) in file:
            in_files.append(os.path.join(r, file))


model=os.path.basename(os.path.normpath(path))
outpath=args.outpath
mkdir(outpath+"/"+model)

num_cores = args.threads

giventaxa = args.taxa
with open(giventaxa) as f:
    taxal = f.readlines()

taxa_dict={}
taxa_list = list()
for line in taxal:
    l = line.rstrip()
    #taxa_list.append(l)
    taxa_dict.setdefault(l, [])

# we need to parallelize the selection --> so to make a function working by gene
# by gene we read the fasta and store ID/seq into a dict --> then extract according to the list

def getTaxa(genefilenumber):
    seqs={}
    genefile = in_files[genefilenumber]
    cur_genome = SeqIO.parse(genefile, "fasta")
    fname=os.path.basename(genefile)
    'we read the fasta and store the info ID+seqs into a dict'
    for record in cur_genome:
        seqID=record.id.replace(";","")
        sequence=record.seq
        seqs.setdefault(seqID, []).append(sequence)
    'we seek taxa of interest from the input list'
    for taxa in taxa_dict:
        if taxa in list(seqs.keys()):
            header=">"+taxa
            newseq=str(seqs[taxa][0])
            if os.path.isfile(os.path.join(outpath+"/"+model, fname)):
                with open(os.path.join(outpath+"/"+model, fname), 'a+') as file:
                    old_headers = []
                    end_file=file.tell()
                    file.seek(0)
                    for line in file:
                        if line.startswith(">"):
                            old_headers.append(line.replace(">","").split(";")[0])
                    if not taxa in old_headers:
                        file.seek(end_file)
                        file.write(header+'\n')
                        file.write(str(newseq)+'\n')
                    else:
                        pass
            else :
                with open(os.path.join(outpath+"/"+model, fname), 'w') as out:
                    out.write(header+'\n')
                    out.write(str(newseq)+'\n')
        else:
            pass


print ("Selection of taxa for %s fasta files" % model)
inputs=range(len(in_files))
Parallel(n_jobs=num_cores)(delayed(getTaxa)(i) for i in inputs)
