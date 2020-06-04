#!/usr/bin/env python

import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.Alphabet import *
from Bio.SeqRecord import *
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA;
from Bio.Seq import *
from Bio.SeqUtils import *

import argparse
import sys
import os, errno


parser = argparse.ArgumentParser(description='Extract genes for annotation files. Script was writen by C. Pouchon (2019).')
parser.add_argument("-in","--infile", help="input assembled annotated files / list",
                    type=str)
parser.add_argument("-m","--model", help="molecular type position",
                    type=str,choices=["chloroplast", "mitochondrion","nucrdna"])
parser.add_argument("-fmt","--annotationfmt", help="format of annotated file",choices=["embl", "genbank"],
                    type=str)
parser.add_argument("-o","--outdir", help="out directory path",
                    type=str)
parser.add_argument("--single", help="[mode] parse a single annotated file instead of a list",
                    action="store_true")

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


'#########################################################'
'parse all input files'
'#########################################################'

model=args.model


modelextens=''
if model == "chloroplast":
    modelextens = '.chloro'
if model == "mitochondrion":
    modelextens = '.mito'
if model == "nucrdna":
    modelextens = '.rdnanuc'

formatin = args.annotationfmt
formatextens=''
if formatin == "embl":
    formatextens = '.embl'
elif formatin == "genbank":
    formatextens = '.gb'

outdirect = args.outdir

stored = {}

'#########################################################'
'extraction of genes according to regions and model'
'#########################################################'

mkdir(outdirect)
mkdir(str(outdirect)+"/"+str(model)+"_misc_RNA"+"/")
mkdir(str(outdirect)+"/"+str(model)+"_rRNA"+"/")

if args.single:
    f=args.infile
    file_path = f
    file_name = os.path.basename(file_path)
    cur_genome = SeqIO.parse(f, formatin)

    for record in cur_genome:
        Bar = ProgressBar(len(record.features), 60, '\t Extraction of sequences')
        barp=0
        for feat in record.features:

            if feat.type == 'source':
                if "/" in feat.qualifiers['organism'][0]:
                    break #verify this process
                else:
                    org=feat.qualifiers['organism'][0].split(" ")
                    taxid=feat.qualifiers['db_xref'][0].split(":")[1]
                    taxaname=str(taxid+"_"+"_".join(org))

            if "misc_" in feat.type:
                if 'gene' in feat.qualifiers:
                    gene = feat.qualifiers['gene'][0].replace(" ","")
                    barp=barp+1
                    if feat.location_operator=="join":
                        length=0
                        out_seq=""
                        joinlocation=str(feat.location).replace("join{","").replace("}","").replace("(-)","").replace("(+)","").replace("]","").split(",")
                        for parts in joinlocation:
                            clean=parts.split("[")[1]
                            if ">" in clean or "<" in clean:
                                pass
                            else:
                                spart = int(float(clean.split(":")[0]))
                                epart = int(float(clean.split(":")[1]))
                                strand = feat.location.strand
                                flanked = FeatureLocation(spart, epart, strand)
                                seq = flanked.extract(record.seq)
                                length = length+len(str(seq))
                                out_seq=out_seq+seq
                    else:
                        s, e, strand = feat.location.start, feat.location.end, feat.location.strand
                        length = e-s
                        flanked = FeatureLocation(s, e, strand)
                        out_seq = flanked.extract(record.seq)

                    header = ">"+str(gene)+"_"+str(taxid)+"_"+"_".join(org)
                    fname = taxaname+".fa"
                    subdir=str(model)+"_misc_RNA"

                    'check for nucleotidic sequence if existing or not'
                    if os.path.isfile(os.path.join(str(outdirect),str(subdir), fname)):
                        with open(os.path.join(str(outdirect),str(subdir), fname), 'a+') as file:
                            old_headers = []
                            end_file=file.tell()
                            file.seek(0)
                            for line in file:
                                if line.startswith(">"):
                                    old_headers.append(line.rstrip())

                            if not header in old_headers:
                                Bar.update(barp)
                                file.seek(end_file)
                                file.write(header+'\n')
                                file.write(str(out_seq)+'\n')
                            else:
                                pass
                    else :
                        with open(os.path.join(str(outdirect),str(subdir), fname), 'a') as out:
                            Bar.update(barp)
                            out.write(header+'\n')
                            out.write(str(out_seq)+'\n')
                else:
                    pass

            elif feat.type == "rRNA":
                if 'gene' in feat.qualifiers:
                    gene = feat.qualifiers['gene'][0].replace(" ","")
                    barp=barp+1
                    if feat.location_operator=="join":
                        length=0
                        out_seq=""
                        joinlocation=str(feat.location).replace("join{","").replace("}","").replace("(-)","").replace("(+)","").replace("]","").split(",")
                        for parts in joinlocation:
                            clean=parts.split("[")[1]
                            if ">" in clean or "<" in clean:
                                pass
                            else:
                                spart = int(float(clean.split(":")[0]))
                                epart = int(float(clean.split(":")[1]))
                                strand = feat.location.strand
                                flanked = FeatureLocation(spart, epart, strand)
                                seq = flanked.extract(record.seq)
                                length = length+len(str(seq))
                                out_seq=out_seq+seq
                    else:
                        s, e, strand = feat.location.start, feat.location.end, feat.location.strand
                        length = e-s
                        flanked = FeatureLocation(s, e, strand)
                        out_seq = flanked.extract(record.seq)

                    header = ">"+str(gene)+"_"+str(taxid)+"_"+"_".join(org)
                    fname = taxaname+".fa"
                    subdir=str(model)+"_"+str("rRNA")

                    'check for nucleotidic sequence if existing or not'
                    if os.path.isfile(os.path.join(str(outdirect),str(subdir), fname)):
                        with open(os.path.join(str(outdirect),str(subdir), fname), 'a+') as file:
                            old_headers = []
                            end_file=file.tell()
                            file.seek(0)
                            for line in file:
                                if line.startswith(">"):
                                    old_headers.append(line.rstrip())

                            if not header in old_headers:
                                Bar.update(barp)
                                file.seek(end_file)
                                file.write(header+'\n')
                                file.write(str(out_seq)+'\n')
                            else:
                                pass
                    else :
                        with open(os.path.join(str(outdirect),str(subdir), fname), 'a') as out:
                            Bar.update(barp)
                            out.write(header+'\n')
                            out.write(str(out_seq)+'\n')
                else:
                    pass

else:
    input_list_taxa = args.infile
    with open(input_list_taxa) as f:
        in_files = f.readlines()

    for files in in_files:
        f = files.rstrip()
        file_path = f
        file_name = os.path.basename(file_path)

        cur_genome = SeqIO.parse(f, formatin)

        for record in cur_genome:
            Bar = ProgressBar(len(record.features), 60, '\t Extraction of sequences')
            barp=0
            for feat in record.features:
                # we seek type "CDS" in the annotation file

                if feat.type == 'source':
                    org=feat.qualifiers['organism'][0].split(" ")
                    taxid=feat.qualifiers['db_xref'][0].split(":")[1]
                    taxaname=str(taxid+"_"+"_".join(org))

                if feat.type == "CDS":
                    if 'gene' in feat.qualifiers:
                        gene = feat.qualifiers['gene'][0]
                        barp=barp+1
                        header = ">"+str(gene)+"_"+str(taxid)+"_"+"_".join(org)
                        if 'translation' in feat.qualifiers:
                            out_seq = feat.qualifiers['translation'][0]
                            fname = taxaname+".fa"
                            subdir=str(model)+"_"+str("CDS")

                            'check for nucleotidic sequence if existing or not'
                            if os.path.isfile(os.path.join(str(outdirect),str(subdir), fname)):
                                with open(os.path.join(str(outdirect),str(subdir), fname), 'a+') as file:
                                    old_headers = []
                                    end_file=file.tell()
                                    file.seek(0)
                                    for line in file:
                                        if line.startswith(">"):
                                            old_headers.append(line.rstrip())

                                    if not header in old_headers:
                                        Bar.update(barp)
                                        file.seek(end_file)
                                        file.write(header+'\n')
                                        file.write(str(out_seq)+'\n')
                                    else:
                                        pass
                            else :
                                with open(os.path.join(str(outdirect),str(subdir), fname), 'a') as out:
                                    Bar.update(barp)
                                    out.write(header+'\n')
                                    out.write(str(out_seq)+'\n')
                        else:
                            pass

                    else:
                        pass

                elif feat.type == "tRNA":
                    if 'gene' in feat.qualifiers:
                        gene = feat.qualifiers['gene'][0]
                        barp=barp+1
                        if "anticodon" in feat.qualifiers:
                            code=feat.qualifiers['anticodon'][0]
                            gene_codon=gene+"-"+code
                        else:
                            gene_codon=gene
                        if feat.location_operator=="join":
                            length=0
                            out_seq=""
                            joinlocation=str(feat.location).replace("join{","").replace("}","").replace("(-)","").replace("(+)","").replace("]","").split(",")
                            for parts in joinlocation:
                                clean=parts.split("[")[1]
                                if ">" in clean or "<" in clean:
                                    pass
                                else:
                                    spart = int(float(clean.split(":")[0]))
                                    epart = int(float(clean.split(":")[1]))
                                    strand = feat.location.strand
                                    flanked = FeatureLocation(spart, epart, strand)
                                    seq = flanked.extract(record.seq)
                                    length = length+len(str(seq))
                                    out_seq=out_seq+seq
                        else:
                            s, e, strand = feat.location.start, feat.location.end, feat.location.strand
                            length = e-s
                            flanked = FeatureLocation(s, e, strand)
                            out_seq = flanked.extract(record.seq)

                        header = ">"+str(gene_codon)+"_"+str(taxid)+"_"+"_".join(org)
                        fname = taxaname+".fa"
                        subdir=str(model)+"_"+str("tRNA")

                        'check for nucleotidic sequence if existing or not'
                        if os.path.isfile(os.path.join(str(outdirect),str(subdir), fname)):
                            with open(os.path.join(str(outdirect),str(subdir), fname), 'a+') as file:
                                old_headers = []
                                end_file=file.tell()
                                file.seek(0)
                                for line in file:
                                    if line.startswith(">"):
                                        old_headers.append(line.rstrip())

                                if not header in old_headers:
                                    Bar.update(barp)
                                    file.seek(end_file)
                                    file.write(header+'\n')
                                    file.write(str(out_seq)+'\n')
                                else:
                                    pass
                        else :
                            with open(os.path.join(str(outdirect),str(subdir), fname), 'a') as out:
                                Bar.update(barp)
                                out.write(header+'\n')
                                out.write(str(out_seq)+'\n')
                    else:
                        pass


                elif feat.type == "rRNA":
                    if 'gene' in feat.qualifiers:
                        gene = feat.qualifiers['gene'][0]
                        barp=barp+1
                        if feat.location_operator=="join":
                            length=0
                            out_seq=""
                            joinlocation=str(feat.location).replace("join{","").replace("}","").replace("(-)","").replace("(+)","").replace("]","").split(",")
                            for parts in joinlocation:
                                clean=parts.split("[")[1]
                                if ">" in clean or "<" in clean:
                                    pass
                                else:
                                    spart = int(float(clean.split(":")[0]))
                                    epart = int(float(clean.split(":")[1]))
                                    strand = feat.location.strand
                                    flanked = FeatureLocation(spart, epart, strand)
                                    seq = flanked.extract(record.seq)
                                    length = length+len(str(seq))
                                    out_seq=out_seq+seq
                        else:
                            s, e, strand = feat.location.start, feat.location.end, feat.location.strand
                            length = e-s
                            flanked = FeatureLocation(s, e, strand)
                            out_seq = flanked.extract(record.seq)

                        header = ">"+str(gene)+"_"+str(taxid)+"_"+"_".join(org)
                        fname = taxaname+".fa"
                        subdir=str(model)+"_"+str("rRNA")

                        'check for nucleotidic sequence if existing or not'
                        if os.path.isfile(os.path.join(str(outdirect),str(subdir), fname)):
                            with open(os.path.join(str(outdirect),str(subdir), fname), 'a+') as file:
                                old_headers = []
                                end_file=file.tell()
                                file.seek(0)
                                for line in file:
                                    if line.startswith(">"):
                                        old_headers.append(line.rstrip())

                                if not header in old_headers:
                                    Bar.update(barp)
                                    file.seek(end_file)
                                    file.write(header+'\n')
                                    file.write(str(out_seq)+'\n')
                                else:
                                    pass
                        else :
                            with open(os.path.join(str(outdirect),str(subdir), fname), 'a') as out:
                                Bar.update(barp)
                                out.write(header+'\n')
                                out.write(str(out_seq)+'\n')
                    else:
                        pass
