#!/usr/bin/env python

import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
#from Bio.Alphabet import *
from Bio.SeqRecord import *
#from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA;
from Bio.Seq import *
from Bio.SeqUtils import *
from joblib import Parallel, delayed
import multiprocessing

import argparse
import sys
import os, errno


parser = argparse.ArgumentParser(description='Extract genes for an annotated file. Script was writen by C. Pouchon (2019).')
parser.add_argument("-t","--intaxa", help="input assembled annotated file",
                    type=str)
parser.add_argument("-n","--taxaname", help="name of sample",
                    type=str)
parser.add_argument("-m","--model", help="molecular type position",
                    type=str,choices=["chloroplast", "mitochondrion","nucrdna"])
parser.add_argument("-g","--geneslist", help="list of genes to extract",
                    type=str)
parser.add_argument("--threads", help="number of threads to use",
                    type=int)
parser.add_argument("-fmt","--annotationfmt", help="format of annotated file",choices=["embl", "genbank"],
                    type=str)
parser.add_argument("-o","--outdir", help="out directory path",
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


def GeneExtraction(genenumber):
    g=genes[genenumber]
    g_tab = g.rstrip().split('\t')
    cur_genome = SeqIO.parse(f, formatin)
    for record in cur_genome:
        for feat in record.features:
                if feat.type == g_tab[0]:
                    if 'gene' in feat.qualifiers:
                        gene = feat.qualifiers['gene'][0]
                        if (gene == g_tab[1]) or (gene == g_tab[1]+"_1") or (gene in g_tab[1]):
                            'condition: if gene was the same as in list'
                            'we seek if gene is a tRNA type to extract the codon OR a CDS to keep also proteic sequence'
                            'AND we check for the curent record if we have only one copy of this gene'
                            if feat.type == "tRNA" :
                                code=feat.qualifiers['anticodon'][0]
                                gene_codon=gene+"-"+code
                                if gene_codon == g_tab[1]:
                                    if feat.location_operator=="join":
                                        length=0
                                        out_seq=""
                                        joinlocation=str(feat.location).replace("join{","").replace("}","").replace("(-)","").replace("(+)","").replace("]","").split(",")
                                        spart=0
                                        epart=0
                                        for parts in joinlocation:
                                            clean=parts.split("[")[1]
                                            if ">" in clean or "<" in clean:
                                                pass
                                            else:
                                                mp = int(float(clean.split(":")[0]))
                                                Mp = int(float(clean.split(":")[1]))
                                                if spart==0 and epart==0:
                                                    spart = mp
                                                    epart = Mp
                                                else:
                                                    if mp > spart:
                                                        spart=spart
                                                        epart=Mp
                                                    else:
                                                        spart=mp
                                                        epart=epart

                                        strand = feat.location.strand
                                        flanked = FeatureLocation(spart, epart, strand)
                                        out_seq = flanked.extract(record.seq)
                                        length = epart-spart

                                        header = ">"+str(name)+"; "+"gene="+str(gene_codon)+"; "+"type="+str(model)+"; "+"length="+str(length)+"; "+"taxa="+taxaname+"; "+"taxid="+taxid_num+"; "+"count="+str(dict_genes[gene_codon])
                                        fname = gene_codon+".fa"
                                        typeofgene = g_tab[0]

                                        if os.path.isfile(os.path.join(outpath+"/"+model+'_'+typeofgene, fname)):
                                            with open(os.path.join(outpath+"/"+model+'_'+typeofgene, fname), 'a+') as file:
                                                old_headers = []
                                                end_file=file.tell()
                                                file.seek(0)
                                                for line in file:
                                                    if line.startswith(">"):
                                                        old_headers.append(line.replace(">","").split(";")[0])
                                                if not name in old_headers:
                                                    file.seek(end_file)
                                                    file.write(header+'\n')
                                                    file.write(str(out_seq)+'\n')
                                                else:
                                                    pass
                                        else :
                                            with open(os.path.join(outpath+"/"+model+'_'+typeofgene, fname), 'a') as out:
                                                out.write(header+'\n')
                                                out.write(str(out_seq)+'\n')
                                    else:
                                        continue
                                else:
                                    pass

                            elif feat.type == "CDS" :
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

                                header = ">"+str(name)+"; "+"gene="+g_tab[1]+"; "+"type="+str(model)+"; "+"length="+str(length)+"; "+"taxa="+taxaname+"; "+"taxid="+taxid_num+"; "+"count="+str(dict_genes[gene])

                                fname = g_tab[1]+".fa"
                                typeofgene = g_tab[0]


                                'check for nucleotidic sequence if existing or not'
                                if os.path.isfile(os.path.join(outpath+"/"+model+'_'+typeofgene, fname)):
                                    with open(os.path.join(outpath+"/"+model+'_'+typeofgene, fname), 'a+') as file:
                                        old_headers = []
                                        end_file=file.tell()
                                        file.seek(0)
                                        for line in file:
                                            if line.startswith(">"):
                                                old_headers.append(line.replace(">","").split(";")[0])
                                        if not name in old_headers:
                                            file.seek(end_file)
                                            file.write(header+'\n')
                                            file.write(str(out_seq)+'\n')
                                        else:
                                            pass
                                else :
                                    with open(os.path.join(outpath+"/"+model+'_'+typeofgene, fname), 'a') as out:
                                        out.write(header+'\n')
                                        out.write(str(out_seq)+'\n')

                            else:
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

                                if g_tab[1]=="5.8S rRNA":
                                    header = ">"+str(name)+"; "+"gene="+"rrn5.8S"+"; "+"type="+str(model)+"; "+"length="+str(length)+"; "+"taxa="+taxaname+"; "+"taxid="+taxid_num+"; "+"count="+str(dict_genes[gene])
                                    fname = "rrn5.8S"+".fa"
                                elif g_tab[1]=="rrnITS1":
                                    header = ">"+str(name)+"; "+"gene="+"ITS1"+"; "+"type="+str(model)+"; "+"length="+str(length)+"; "+"taxa="+taxaname+"; "+"taxid="+taxid_num+"; "+"count="+str(dict_genes[gene])
                                    fname = "ITS1"+".fa"
                                elif g_tab[1]=="rrnITS2":
                                    header = ">"+str(name)+"; "+"gene="+"ITS2"+"; "+"type="+str(model)+"; "+"length="+str(length)+"; "+"taxa="+taxaname+"; "+"taxid="+taxid_num+"; "+"count="+str(dict_genes[gene])
                                    fname = "ITS2"+".fa"
                                else:
                                    header = ">"+str(name)+"; "+"gene="+g_tab[1]+"; "+"type="+str(model)+"; "+"length="+str(length)+"; "+"taxa="+taxaname+"; "+"taxid="+taxid_num+"; "+"count="+str(dict_genes[gene])
                                    fname = g_tab[1]+".fa"

                                typeofgene = g_tab[0]

                                if model=="nucrdna":
                                    if os.path.isfile(os.path.join(outpath+"/"+model, fname)):
                                        with open(os.path.join(outpath+"/"+model, fname), 'a+') as file:
                                            old_headers = []
                                            end_file=file.tell()
                                            file.seek(0)
                                            for line in file:
                                                if line.startswith(">"):
                                                    old_headers.append(line.replace(">","").split(";")[0])
                                            if not name in old_headers:
                                                file.seek(end_file)
                                                file.write(header+'\n')
                                                file.write(str(out_seq)+'\n')
                                            else:
                                                pass
                                    else :
                                        with open(os.path.join(outpath+"/"+model, fname), 'w') as out:
                                            out.write(header+'\n')
                                            out.write(str(out_seq)+'\n')
                                else:
                                    if os.path.isfile(os.path.join(outpath+"/"+model+'_'+typeofgene, fname)):
                                        with open(os.path.join(outpath+"/"+model+'_'+typeofgene, fname), 'a+') as file:
                                            old_headers = []
                                            end_file=file.tell()
                                            file.seek(0)
                                            for line in file:
                                                if line.startswith(">"):
                                                    old_headers.append(line.replace(">","").split(";")[0])
                                            if not name in old_headers:
                                                file.seek(end_file)
                                                file.write(header+'\n')
                                                file.write(str(out_seq)+'\n')
                                            else:
                                                pass
                                    else :
                                        with open(os.path.join(outpath+"/"+model+'_'+typeofgene, fname), 'w') as out:
                                            out.write(header+'\n')
                                            out.write(str(out_seq)+'\n')

                    else:
                        pass
                else:
                    pass

'#########################################################'
'parse all input files'
'#########################################################'

input = args.intaxa
input_list_genes = args.geneslist

with open(input_list_genes) as f:
    genes = f.readlines()


name = args.taxaname
model=args.model

outpath=args.outdir

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

stored = {}

num_cores = args.threads

'#########################################################'
'extraction of genes according to regions and model'
'#########################################################'

if model in ("chloroplast","nucrdna","mitochondrion"):
    '''
    mkdir for differents regions of a genes list for a given model
    '''
    types=[]
    for g in genes:
        g_tab = g.rstrip().split('\t')
        types.append(g_tab[0])
    if model=="nucrdna":
        mkdir(str(outpath+"/"+model))
    else:
        for type in set(types):
            mkdir(str(outpath+"/"+model+'_'+type))
    '''
    process the extraction
    '''
    f = input.rstrip()
    file_path = f
    file_name = os.path.basename(file_path)
    index_of_dot = file_name.index('.')
    file_name_without_extension = file_name[:index_of_dot]
    taxaname = file_name_without_extension.split(":")[0]
    if len(file_name_without_extension.split(":")) > 1:
        taxid_num = file_name_without_extension.split(":")[1]
    else:
        taxid_num = "NA"

    dict_genes = {}
    for g in genes:
        g_tab = g.rstrip().split('\t')
        cur_genome = SeqIO.parse(f, formatin)
        for record in cur_genome:
            for feat in record.features:
                if feat.type == g_tab[0]:
                    if 'gene' in feat.qualifiers:
                        gene = feat.qualifiers['gene'][0]
                        if (gene == g_tab[1]) or (gene == g_tab[1]+"_1") or (gene in g_tab[1]):
                            if feat.type == "tRNA" :
                                code=feat.qualifiers['anticodon'][0]
                                gene_codon=gene+"-"+code
                                if gene_codon==g_tab[1]:
                                    if gene_codon in dict_genes.keys():
                                        dict_genes[gene_codon] = dict_genes[gene_codon]+1
                                    else:
                                        dict_genes[gene_codon] = 1
                            else:
                                if gene in dict_genes.keys():
                                    dict_genes[gene] = dict_genes[gene]+1
                                else:
                                    dict_genes[gene] = 1


    print ("Extraction of sequences for %s" % file_name_without_extension)
    inputs=range(len(genes))
    Parallel(n_jobs=num_cores)(delayed(GeneExtraction)(i) for i in inputs)

else:
    print ("model not recognized - Please, re-run with model=chloroplast,mitochondrion,nucrdna")
    pass
