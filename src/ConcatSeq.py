#! /usr/bin/env python

import glob
import sys
import os, errno
import argparse
import fnmatch
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation




parser = argparse.ArgumentParser(description='Concatenation of sequences alignments')
parser.add_argument("-p","--path", help="searching path of sequences files",
                    type=str)
parser.add_argument("-i","--infile", help="input file list of sequence files, sequences have to be aligned before",
                    type=str)
parser.add_argument("-pfind","--pathfind", help="[mode] search all sequences files in given path (-p) otherwise parse a given list (-i)",
                    action="store_true")
parser.add_argument("-o","--outfile", help="output file list of concatenated sequences",
                    type=str)
parser.add_argument("-t","--taxa", help="list of taxa to include in the final concatenated file",
                    type=str)
parser.add_argument("--target", help="list of targeted gene directory",
                    type=str)
parser.add_argument("-g","--genes", help="list of gene to include in the final concatenated file",
                    type=str)
parser.add_argument("--trim", help="[mode] if alignments files were trimmed, genes are collected in trim/ folder in rigth target. Otherwise, they are collected in files/",
                    action="store_true")
parser.add_argument("--missingfract", help="maximal missing data threshold allowed to consider final sequence (e.g. 0.5 meaning that final sequence has fewer than 0.5 of missing data)",
                     type=float,default=None)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()


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



if args.pathfind:
    list_target = [str(item) for item in args.target.split(',')]
    'Search all fasta files files'
    path = args.path
    in_files = []
    for p in list_target:
        if args.trim:
            path_to_seek=os.walk(os.path.join(path,p,"trim"))
            for r, d, f in path_to_seek:
                for file in f:
                    if str("."+"trim") in file:
                        in_files.append(os.path.join(r, file))
        else:
            path_to_seek=os.walk(os.path.join(path,p,"files"))
            for r, d, f in path_to_seek:
                for file in f:
                    if str("."+"fa") in file:
                        in_files.append(os.path.join(r, file))
else:
    'We parse a list of files'
    input_list = args.infile
    with open(input_list) as f:
        in_files = f.readlines()

#cat *.fa | grep '>' | sort | uniq | perl -pe 's/>//' > list_taxa

fname=args.outfile+".fa"
Nthreshold=args.missingfract

giventaxa = args.taxa
with open(giventaxa) as f:
    taxal = f.readlines()

taxa_dict={}
taxa_list = list()
for line in taxal:
    l = line.rstrip()
    if len(l)>0:
        taxa_dict.setdefault(l, {})

givengenes = args.genes
with open(givengenes) as f:
    genel = f.readlines()

gene_list = list()
for gene in genel:
    l = gene.rstrip()
    if len(l)>0:
        gene_list.append(l)

end=0
partnumber=0

for file in in_files:
    seqs={}
    if args.trim:
        genename=os.path.basename(file).replace(str("."+"trim"),"")
    else:
        genename=os.path.basename(file).replace(str("."+"fa"),"")

    if genename in gene_list:
        cur_genome = SeqIO.parse(file, "fasta")
        for record in cur_genome:
            if "_R_" in record.id.split(";")[0]:
                seqID=record.id.split(";")[0].replace("_R_","")
            else:
                seqID=record.id.split(";")[0]
            sequence=record.seq.upper()

            if seqID not in list(seqs.keys()):
                seqs[seqID]=dict()
                if genename not in list(seqs[seqID].keys()):
                    seqs[seqID][genename]=sequence
                else:
                    pass
            else:
                if genename not in list(seqs[seqID].keys()):
                    seqs[seqID][genename]=sequence
                else:
                    pass

        start=end+1
        end=(start-1)+len(str(seqs[list(seqs.keys())[0]][genename]))
        partnumber=partnumber+1

        splitfile=os.path.split(os.path.abspath(fname))
        partfile=str(splitfile[0]+'/'+splitfile[1].split(".")[0]+".partitions")
        partinfo=str(splitfile[0]+'/'+splitfile[1].split(".")[0]+".info")

        if os.path.isfile(partfile):
            with open(partfile, 'a+') as file:
                file.write(str("DNA, part"+str(partnumber)+" = "+str(start)+"-"+str(end)+"\n"))
        else :
            with open(partfile, 'w') as out:
                out.write(str("DNA, part"+str(partnumber)+" = "+str(start)+"-"+str(end)+"\n"))

        if os.path.isfile(partinfo):
            with open(partinfo, 'a+') as file:
                file.write(str(str(start)+"\t"+str(end)+"\t"+str(genename)+"\t"+"part"+str(partnumber)+"\n"))
        else :
            with open(partinfo, 'w') as out:
                out.write(str(str(start)+"\t"+str(end)+"\t"+str(genename)+"\t"+"part"+str(partnumber)+"\n"))

        for taxa in list(taxa_dict.keys()):
            if taxa in list(seqs.keys()):
                taxa_dict[taxa][genename]=seqs[taxa][genename]
            else:
                'we create a null sequence for missing taxa per gene'
                lenseq = len(str(seqs[list(seqs.keys())[0]][genename]))
                nullseq = "N"*lenseq
                taxa_dict[taxa][genename]=nullseq
    else:
        continue


'concatenation of sequences in dictionary and output sequences'
for taxa in taxa_dict:
    header= ">"+str(taxa)
    concatenate_seq=""
    for file in in_files:
        if args.trim:
            genename=os.path.basename(file).replace(str("."+"trim"),"")
        else:
            genename=os.path.basename(file).replace(str("."+"fa"),"")
        if genename in gene_list:
            concatenate_seq=concatenate_seq+str(taxa_dict[taxa][genename])
    Ncount=int(concatenate_seq.count("N"))
    gappcount=int(concatenate_seq.count("-"))
    SeqLength=len(concatenate_seq)
    Nratio=float(float(Ncount)/float(SeqLength))
    AmbCount=Ncount+gappcount
    Ambratio=float(float(AmbCount)/float(SeqLength))
    'we add condition to check if a taxa is completly missing'
    seqtestN='N'*len(concatenate_seq)
    splitfile=os.path.split(os.path.abspath(fname))
    missfile=str(splitfile[0]+'/'+splitfile[1].split(".")[0]+".missingdata")
    with open(missfile, 'a+') as file:
        file.write(str(str(taxa)+"\t"+str(Ambratio)+"\n"))
    if seqtestN == concatenate_seq:
        print ("WARN: %s is missing" % (taxa))
    elif Nratio > Nthreshold:
        print ("WARN: %s has too missing data, not passed filter fraction of %s" % (taxa,Nthreshold))
        if os.path.isfile("missing_seqs.fasta"):
            'we stored only the missing sequences in a file'
            with open("missing_seqs.fasta", 'a+') as file:
                old_headers = []
                end_file=file.tell()
                file.seek(0)
                for line in file:
                    if line.startswith(">"):
                        old_headers.append(line.replace(">","").split(";")[0])
                if not taxa in old_headers:
                    file.seek(end_file)
                    file.write(header+'\n')
                    file.write(str(concatenate_seq)+'\n')
                else:
                    pass
        else :
            with open("missing_seqs.fasta", 'w') as out:
                out.write(header+'\n')
                out.write(str(concatenate_seq)+'\n')
    else:
        if os.path.isfile(fname):
            with open(fname, 'a+') as file:
                old_headers = []
                end_file=file.tell()
                file.seek(0)
                for line in file:
                    if line.startswith(">"):
                        old_headers.append(line.replace(">","").split(";")[0])
                if not taxa in old_headers:
                    file.seek(end_file)
                    file.write(header+'\n')
                    file.write(str(concatenate_seq)+'\n')
                else:
                    pass
        else :
            with open(fname, 'w') as out:
                out.write(header+'\n')
                out.write(str(concatenate_seq)+'\n')
