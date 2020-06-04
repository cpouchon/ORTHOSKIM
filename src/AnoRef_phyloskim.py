#!/usr/bin/env python

import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.Alphabet import *
from Bio.SeqRecord import *
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA;
from Bio.Seq import *
from Bio.SeqUtils import *
from joblib import Parallel, delayed
import multiprocessing

import argparse
import sys
import os, errno


parser = argparse.ArgumentParser(description='Extract genes for an annotated file. Script was writen by C. Pouchon (2019).')
parser.add_argument("-i","--infile", help="input assembled annotated file",
                    type=str)
parser.add_argument("-m","--model", help="molecular type position",
                    type=str,choices=["chloroplast", "mitochondrion","nucrdna"])
parser.add_argument("-g","--geneslist", help="list of genes to extract",
                    type=str)
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


'#########################################################'
'parse all input files'
'#########################################################'

input = args.infile
input_list_genes = args.geneslist

with open(input_list_genes) as f:
    genes = f.readlines()

dict_genes={}
for g in genes:
    g_tab = g.rstrip().split('\t')
    if g_tab[0] not in dict_genes.keys():
        dict_genes[g_tab[0]]=[]
        dict_genes[g_tab[0]].append(g_tab[1])
    else:
        dict_genes[g_tab[0]].append(g_tab[1])

# acess all genes --> [item for sublist in dict_genes.values() for item in sublist]

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


'#########################################################'
'extraction of genes according to regions and model'
'#########################################################'

if model in ("chloroplast","nucrdna","mitochondrion"):
    '''
    mkdir for differents regions of a genes list for a given model
    '''

    for type in dict_genes.keys():
        mkdir(str(outpath+"/"+model+'_'+type))
    '''
    process the extraction
    '''
    f = input.rstrip()
    file_path = f
    file_name = os.path.basename(file_path)

    to_annot_file={}
    to_annot_file['seq']=[]

    cur_genome = SeqIO.parse(f, formatin)
    for record in cur_genome:

        break_record=0

        for feat in record.features:
            if feat.type == 'source':
                if "/" in feat.qualifiers['organism'][0]:
                    break_record=break_record+1 #verify this process
                else:
                    org=feat.qualifiers['organism'][0].split(" ")
                    taxid=feat.qualifiers['db_xref'][0].split(":")[1]
                    taxaname=str(taxid+"_"+"_".join(org))
                    orga_info=feat.qualifiers['organism']
                    db_info=feat.qualifiers['db_xref']

            if break_record>0:
                break
            else:
                if feat.type in dict_genes.keys():
                    gene_type=feat.type
                    if 'gene' in feat.qualifiers:
                        if "_1" in feat.qualifiers['gene'][0]:
                            gene = feat.qualifiers['gene'][0].replace("_1","")
                        elif "_2" in feat.qualifiers['gene'][0]:
                            gene = feat.qualifiers['gene'][0].replace("_2","")
                        elif feat.qualifiers['gene'][0]=="5.8S rRNA":
                            gene = "rrn5.8S"
                        else:
                            gene = feat.qualifiers['gene'][0]
                        if feat.type == "tRNA" :
                            if (gene == dict_genes["tRNA"][0].split("-")[0]):
                                if "anticodon" in feat.qualifiers:
                                    code=feat.qualifiers['anticodon'][0]
                                    if code==dict_genes["tRNA"][0].split("-")[1]:
                                        gene_codon=gene+"-"+code

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
                                        fname = gene_codon+".fa"
                                        typeofgene = "tRNA"
                                        if os.path.isfile(os.path.join(outpath+"/"+model+'_'+typeofgene, fname)):
                                            with open(os.path.join(outpath+"/"+model+'_'+typeofgene, fname), 'a+') as file:
                                                old_headers = []
                                                end_file=file.tell()
                                                file.seek(0)
                                                for line in file:
                                                    if line.startswith(">"):
                                                        old_headers.append(line.rstrip())
                                                if not header in old_headers:
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
                                        pass

                        elif feat.type == "CDS" :
                            if gene in [item for sublist in dict_genes.values() for item in sublist]:
                                header = ">"+str(gene)+"_"+str(taxid)+"_"+"_".join(org)
                                if 'translation' in feat.qualifiers:
                                    out_seq = feat.qualifiers['translation'][0]
                                    fname = gene+".fa"
                                    typeofgene = "CDS"

                                    'check for nucleotidic sequence if existing or not'
                                    if os.path.isfile(os.path.join(outpath+"/"+model+'_'+typeofgene, fname)):
                                        with open(os.path.join(outpath+"/"+model+'_'+typeofgene, fname), 'a+') as file:
                                            old_headers = []
                                            end_file=file.tell()
                                            file.seek(0)
                                            for line in file:
                                                if line.startswith(">"):
                                                    old_headers.append(line.rstrip())
                                            if not header in old_headers:
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
                                    pass

                        else:
                            if gene in [item for sublist in dict_genes.values() for item in sublist]:
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
                                fname = gene+".fa"
                                typeofgene=feat.type

                                if 'ITS' in gene:
                                    pass
                                else:
                                    if os.path.isfile(os.path.join(outpath+"/"+model+'_'+typeofgene, fname)):
                                        with open(os.path.join(outpath+"/"+model+'_'+typeofgene, fname), 'a+') as file:
                                            old_headers = []
                                            end_file=file.tell()
                                            file.seek(0)
                                            for line in file:
                                                if line.startswith(">"):
                                                    old_headers.append(line.rstrip())
                                            if not header in old_headers:
                                                file.seek(end_file)
                                                file.write(header+'\n')
                                                file.write(str(out_seq)+'\n')
                                            else:
                                                pass
                                    else :
                                        with open(os.path.join(outpath+"/"+model+'_'+typeofgene, fname), 'w') as out:
                                            out.write(header+'\n')
                                            out.write(str(out_seq)+'\n')


        if break_record==0:
            sequence=str(record.seq)
            to_annot_file['seq'].append(sequence)
            to_annot_file['org']=orga_info
            to_annot_file['db']=db_info
            to_annot_file['id']=record.id
            to_annot_file['name']=record.name
            to_annot_file['description']=record.description
        else:
            continue


    concat_seq="NNNNNNNNNNNNNNNNNNNN".join(to_annot_file['seq'])

    if model=="chloroplast":
        cond_size=140000
    elif model=="nucrdna":
        cond_size=3000
    else:
        cond_size=50000

    if len("".join(to_annot_file['seq'])) >= cond_size:
        newrecord=SeqRecord(seq=Seq(concat_seq), id=to_annot_file['id'], name=to_annot_file['name'],
        description=to_annot_file['description'],dbxrefs=[])
        from Bio.Alphabet import generic_dna
        newrecord.seq.alphabet = generic_dna

        from Bio import SeqFeature
        my_start_pos = SeqFeature.ExactPosition(0)
        my_end_pos = SeqFeature.ExactPosition(len(concat_seq))

        from Bio.SeqFeature import FeatureLocation
        my_feature_location = FeatureLocation(my_start_pos,my_end_pos)
        my_feature_type = "source"

        from Bio.SeqFeature import SeqFeature
        my_feature = SeqFeature(my_feature_location,type=my_feature_type)
        my_feature.location.strand=1
        my_feature.qualifiers['organism']=orga_info
        my_feature.qualifiers['db_xref']=db_info
        newrecord.features.append(my_feature)

        mkdir(str(outpath))

        with open(str(str(outpath)+"/"+model+"_genomes"+formatextens), "a+") as output_handle:
            SeqIO.write(newrecord, output_handle, formatin)
    else:
        pass


    print ("Extraction of sequences for %s" % file_name)

else:
    print ("model not recognized - Please, re-run with model=chloroplast,mitochondrion,nucrdna")
    pass
