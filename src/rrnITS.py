#!/usr/bin/env python

import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import sys
import os, errno
import argparse
from joblib import Parallel, delayed
import multiprocessing

from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna


parser = argparse.ArgumentParser(description='Create ITS probes. Script was writen by C. Pouchon (2019).')
parser.add_argument("-p","--inpath", help="path to clean nucrdna rRNA genes fasta files",
                    type=str)
parser.add_argument("-o","--outpath", help="Out directory path to clean nucrdna misc_RNA",
                    type=str)

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



opath=args.outpath
path = args.inpath
#path = "./references/clean/nucrdna_rRNA/"

in_files = []
for r, d, f in os.walk(path):
    for file in f:
        if '.fa' in file:
            if os.path.join(r, file) not in in_files:
                in_files.append(os.path.join(r,file))

stored={}

mkdir(opath)

for file in in_files:
    cur_genome = SeqIO.parse(file, "fasta")
    for record in cur_genome:
        geneid = record.id.split("_")[0]
        seqID="_".join(record.id.split("_")[1:len(record.id.split("_"))])
        sequence=record.seq
        if seqID not in stored.keys():
            stored[seqID]=dict()
            if geneid not in stored[seqID].keys():
                stored[seqID][geneid]=sequence
            else:
                pass
        else:
            stored[seqID][geneid]=sequence

for taxa in stored.keys():
    if 'rrn18S' in stored[taxa].keys() and 'rrn5.8S' in stored[taxa].keys():
        if len(str(stored[taxa]['rrn18S']))>50:
            amorce_18S=str(stored[taxa]['rrn18S'])[-51:len(str(stored[taxa]['rrn18S']))]
        else:
            amorce_18S=str(stored[taxa]['rrn18S'])
        catseq="".join((str(amorce_18S),str(stored[taxa]['rrn5.8S'])))
        genename = "rrnITS1"
        header=">"+genename+"_"+taxa
        fname=taxa+".fa"
        if os.path.isfile(os.path.join(opath, fname)):
            with open(os.path.join(opath, fname), 'a+') as file:
                old_headers = []
                end_file=file.tell()
                file.seek(0)
                for line in file:
                    if line.startswith(">"):
                        old_headers.append(line.replace(">","").split(";")[0])
                if not str(genename+"_"+taxa) in old_headers:
                    file.seek(end_file)
                    file.write(header+'\n')
                    file.write(str(catseq)+'\n')
                else:
                    pass
        else :
            with open(os.path.join(opath, fname), 'w') as out:
                out.write(header+'\n')
                out.write(str(catseq)+'\n')
    if 'rrn28S' in stored[taxa].keys() and 'rrn5.8S' in stored[taxa].keys():
        if len(str(stored[taxa]['rrn28S']))>50:
            amorce_28S=str(stored[taxa]['rrn28S'])[0:50]
        else:
            amorce_28S=str(stored[taxa]['rrn28S'])
        catseq="".join((str(stored[taxa]['rrn5.8S']),str(amorce_28S)))
        genename = "rrnITS2"
        header=">"+genename+"_"+taxa
        fname=taxa+".fa"
        if os.path.isfile(os.path.join(opath, fname)):
            with open(os.path.join(opath, fname), 'a+') as file:
                old_headers = []
                end_file=file.tell()
                file.seek(0)
                for line in file:
                    if line.startswith(">"):
                        old_headers.append(line.replace(">","").split(";")[0])
                if not str(genename+"_"+taxa) in old_headers:
                    file.seek(end_file)
                    file.write(header+'\n')
                    file.write(str(catseq)+'\n')
                else:
                    pass
        else :
            with open(os.path.join(opath, fname), 'w') as out:
                out.write(header+'\n')
                out.write(str(catseq)+'\n')
    elif 'rrn26S' in stored[taxa].keys() and 'rrn5.8S' in stored[taxa].keys():
        if len(str(stored[taxa]['rrn26S']))>50:
            amorce_28S=str(stored[taxa]['rrn26S'])[0:50]
        else:
            amorce_28S=str(stored[taxa]['rrn26S'])
        catseq="".join((str(stored[taxa]['rrn5.8S']),str(amorce_28S)))
        genename = "rrnITS2"
        header=">"+genename+"_"+taxa
        fname=taxa+".fa"
        if os.path.isfile(os.path.join(opath, fname)):
            with open(os.path.join(opath, fname), 'a+') as file:
                old_headers = []
                end_file=file.tell()
                file.seek(0)
                for line in file:
                    if line.startswith(">"):
                        old_headers.append(line.replace(">","").split(";")[0])
                if not str(genename+"_"+taxa) in old_headers:
                    file.seek(end_file)
                    file.write(header+'\n')
                    file.write(str(catseq)+'\n')
                else:
                    pass
        else :
            with open(os.path.join(opath, fname), 'w') as out:
                out.write(header+'\n')
                out.write(str(catseq)+'\n')
