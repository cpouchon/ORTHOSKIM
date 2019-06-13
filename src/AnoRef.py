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
parser.add_argument("-in","--taxalist", help="input assembled annotated files",
                    type=str)
parser.add_argument("-m","--model", help="molecular type position",
                    type=str,choices=["chloroplast", "mitochondrion","nucrdna"])
parser.add_argument("-type","--typeseq", help="type of sequences for extraction",
                    type=str)
parser.add_argument("-fmt","--annotationfmt", help="format of annotated file",choices=["embl", "genbank"],
                    type=str)
parser.add_argument("-o","--outdir", help="out directory path",
                    type=str)
parser.add_argument("-protonly","--only_protein", help="[mode] extraction is made only on CDS type by keeping translated sequences",
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
                print "path '%s' already exists" % path   # overwrite == False and we've hit a directory that exists
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


input_list_taxa = args.taxalist
with open(input_list_taxa) as f:
    in_files = f.readlines()

# input_list_genes = args.geneslist
# with open(input_list_genes) as f:
#     genes = f.readlines()

model=args.model

feattype=args.typeseq

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

if model in ("chloroplast","nucrdna","mitochondrion"):
    '''
    mkdir for differents regions of a genes list for a given model
    '''
    # types=[]
    # for g in genes:
    #     g_tab = g.rstrip().split('\t')
    #     types.append(g_tab[0])

    mkdir(outdirect)
    for files in in_files:
        f = files.rstrip()
        '''
        process the extraction
        '''
        file_path = f
        file_name = os.path.basename(file_path)
        index_of_dot = file_name.index('.')
        file_name_without_extension = file_name[:index_of_dot]
        taxaname = file_name_without_extension.split(":")[0]
        if len(file_name_without_extension.split(":")) > 1:
            taxid_num = file_name_without_extension.split(":")[1]
        else:
            taxid_num = "NA"

        # dict_genes = {}
        # for g in genes:
        #     g_tab = g.rstrip().split('\t')
        #     cur_genome = SeqIO.parse(f, formatin)
        #     for record in cur_genome:
        #         for feat in record.features:
        #             if feat.type == g_tab[0]:
        #                 if 'gene' in feat.qualifiers:
        #                     gene = feat.qualifiers['gene'][0]
        #                     if (gene == g_tab[1]) or (gene == g_tab[1]+"_1"):
        #                         if feat.type == "tRNA" :
        #                             code=feat.qualifiers['anticodon'][0]
        #                             gene_codon=gene+"_"+code
        #                             if gene_codon in dict_genes.keys():
        #                                 dict_genes[gene_codon] = dict_genes[gene_codon]+1
        #                             else:
        #                                 dict_genes[gene_codon] = 1
        #                         else:
        #                             if gene in dict_genes.keys():
        #                                 dict_genes[gene] = dict_genes[gene]+1
        #                             else:
        #                                 dict_genes[gene] = 1

        cur_genome = SeqIO.parse(f, formatin)
        for record in cur_genome:
            Bar = ProgressBar(len(record.features), 60, '\t Extraction of sequences for %s' % file_name_without_extension)
            barp=0
            for feat in record.features:
                # we seek type "CDS" in the annotation file

                if feat.type == 'source':
                    org=feat.qualifiers['organism'][0].split(" ")
                    taxid=feat.qualifiers['db_xref'][0].split(":")[1]
                    taxaname=str(org[0]+"_"+org[1]+"_"+taxid+"_"+formatin)

                if feat.type == feattype:
                    if 'gene' in feat.qualifiers:
                        gene = feat.qualifiers['gene'][0]
                        if args.only_protein:
                            if feat.type == "CDS" :
                                barp=barp+1
                                'gene could be in multiple part if it is partial in multiple contigs'
                                if feat.location_operator=="join":
                                    joinlocation = str(feat.location).replace("join{","").replace("}","").replace("(-)","").replace("(+)","").replace("[","").replace("]","").replace(" ","").split(",")
                                    part1, part2 = joinlocation[0], joinlocation[1]
                                    spart1, epart1 = int(float(part1.split(":")[0])), int(float(part1.split(":")[1]))
                                    spart2, epart2 = int(float(part2.split(":")[0])), int(float(part2.split(":")[1]))
                                    length = (epart1 - spart1)+(epart2-spart2)
                                    strand = feat.location.strand
                                    flanked1 = FeatureLocation(spart1, epart1, strand)
                                    flanked2 = FeatureLocation(spart2, epart2, strand)
                                    out_seq1 = flanked1.extract(record.seq)
                                    out_seq2 = flanked2.extract(record.seq)
                                    out_seq = out_seq1+out_seq2

                                else:
                                    s, e, strand = feat.location.start, feat.location.end, feat.location.strand
                                    length = e-s
                                    flanked = FeatureLocation(s, e, strand)
                                    out_seq = flanked.extract(record.seq)

                                header = ">"+str(taxaname)+"_"+str(gene)+"; "+"gene="+str(gene)+"; "+"type="+str(model)+"; "+"length="+str(length)+"; "+"taxa="+str(org[0]+"_"+org[1])+"; "+"taxid="+taxid_num
                                fname = taxaname+".fna"

                                if 'translation' in feat.qualifiers:
                                    pseq = feat.qualifiers['translation'][0]

                                'check for proteic sequence if existing or not'
                                if 'translation' in feat.qualifiers:
                                    if os.path.isfile(os.path.join(outdirect, fname)):
                                        with open(os.path.join(outdirect, fname), 'a+') as file:
                                            old_headers = []
                                            end_file=file.tell()
                                            file.seek(0)
                                            for line in file:
                                                if line.startswith(">"):
                                                    old_headers.append(line.rstrip())

                                            if not header in old_headers:
                                                file.seek(end_file)
                                                file.write(header+'\n')
                                                file.write(str(pseq)+'\n')
                                            else:
                                                pass
                                    else :
                                        with open(os.path.join(outdirect, fname), 'a') as out:
                                            out.write(header+'\n')
                                            out.write(str(pseq)+'\n')
                            else:
                                continue

                        else:
                            if feat.type == "tRNA" :
                                code=feat.qualifiers['anticodon'][0]
                                gene_codon=gene+"_"+code
                                barp=barp+1
                                if feat.location_operator=="join":
                                    joinlocation = str(feat.location).replace("join{","").replace("}","").replace("(-)","").replace("(+)","").replace("[","").replace("]","").replace(" ","").split(",")
                                    part1, part2 = joinlocation[0], joinlocation[1]
                                    spart1, epart1 = int(float(part1.split(":")[0])), int(float(part1.split(":")[1]))
                                    spart2, epart2 = int(float(part2.split(":")[0])), int(float(part2.split(":")[1]))
                                    length = (epart1 - spart1)+(epart2-spart2)
                                    strand = feat.location.strand
                                    flanked1 = FeatureLocation(spart1, epart1, strand)
                                    flanked2 = FeatureLocation(spart2, epart2, strand)
                                    out_seq1 = flanked1.extract(record.seq)
                                    out_seq2 = flanked2.extract(record.seq)
                                    out_seq = out_seq1+out_seq2

                                else:
                                    s, e, strand = feat.location.start, feat.location.end, feat.location.strand
                                    length = e-s
                                    flanked = FeatureLocation(s, e, strand)
                                    out_seq = flanked.extract(record.seq)

                                header = ">"+str(taxaname)+"_"+str(gene)+"; "+"gene="+str(gene_codon)+"; "+"type="+str(model)+"; "+"length="+str(length)+"; "+"taxa="+str(org[0]+"_"+org[1])+"; "+"taxid="+taxid_num
                                fname = taxaname+".fna"

                                if os.path.isfile(os.path.join(outdirect, fname)):
                                    with open(os.path.join(outdirect, fname), 'a+') as file:
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
                                    with open(os.path.join(outdirect, fname), 'a') as out:
                                        Bar.update(barp)
                                        out.write(header+'\n')
                                        out.write(str(out_seq)+'\n')


                            elif feat.type == "CDS" :
                                barp=barp+1
                                'gene could be in multiple part if it is partial in multiple contigs'
                                if feat.location_operator=="join":
                                    joinlocation = str(feat.location).replace("join{","").replace("}","").replace("(-)","").replace("(+)","").replace("[","").replace("]","").replace(" ","").split(",")
                                    part1, part2 = joinlocation[0], joinlocation[1]
                                    spart1, epart1 = int(float(part1.split(":")[0])), int(float(part1.split(":")[1]))
                                    spart2, epart2 = int(float(part2.split(":")[0])), int(float(part2.split(":")[1]))
                                    length = (epart1 - spart1)+(epart2-spart2)
                                    strand = feat.location.strand
                                    flanked1 = FeatureLocation(spart1, epart1, strand)
                                    flanked2 = FeatureLocation(spart2, epart2, strand)
                                    out_seq1 = flanked1.extract(record.seq)
                                    out_seq2 = flanked2.extract(record.seq)
                                    out_seq = out_seq1+out_seq2

                                else:
                                    s, e, strand = feat.location.start, feat.location.end, feat.location.strand
                                    length = e-s
                                    flanked = FeatureLocation(s, e, strand)
                                    out_seq = flanked.extract(record.seq)

                                header = ">"+str(taxaname)+"_"+str(gene)+"; "+"gene="+str(gene)+"; "+"type="+str(model)+"; "+"length="+str(length)+"; "+"taxa="+str(org[0]+"_"+org[1])+"; "+"taxid="+taxid_num
                                fname = taxaname+".fna"

                                'check for nucleotidic sequence if existing or not'
                                if os.path.isfile(os.path.join(outdirect, fname)):
                                    with open(os.path.join(outdirect, fname), 'a+') as file:
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
                                    with open(os.path.join(outdirect, fname), 'a') as out:
                                        Bar.update(barp)
                                        out.write(header+'\n')
                                        out.write(str(out_seq)+'\n')

                            else:
                                barp=barp+1
                                'gene could be in multiple part if it is partial in multiple contigs'
                                if feat.location_operator=="join":
                                    joinlocation = str(feat.location).replace("join{","").replace("}","").replace("(-)","").replace("(+)","").replace("[","").replace("]","").replace(" ","").split(",")
                                    part1, part2 = joinlocation[0], joinlocation[1]
                                    spart1, epart1 = int(float(part1.split(":")[0])), int(float(part1.split(":")[1]))
                                    spart2, epart2 = int(float(part2.split(":")[0])), int(float(part2.split(":")[1]))
                                    length = (epart1 - spart1)+(epart2-spart2)
                                    strand = feat.location.strand
                                    flanked1 = FeatureLocation(spart1, epart1, strand)
                                    flanked2 = FeatureLocation(spart2, epart2, strand)
                                    out_seq1 = flanked1.extract(record.seq)
                                    out_seq2 = flanked2.extract(record.seq)
                                    out_seq = out_seq1+out_seq2

                                else:
                                    s, e, strand = feat.location.start, feat.location.end, feat.location.strand
                                    length = e-s
                                    flanked = FeatureLocation(s, e, strand)
                                    out_seq = flanked.extract(record.seq)

                                header = ">"+str(taxaname)+"_"+str(gene)+"; "+"gene="+str(gene)+"; "+"type="+str(model)+"; "+"length="+str(length)+"; "+"taxa="+str(org[0]+"_"+org[1])+"; "+"taxid="+taxid_num
                                fname = taxaname+".fna"

                                if os.path.isfile(os.path.join(outdirect, fname)):
                                    with open(os.path.join(outdirect, fname), 'a+') as file:
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
                                    with open(os.path.join(outdirect, fname), 'w') as out:
                                        Bar.update(barp)
                                        out.write(header+'\n')
                                        out.write(str(out_seq)+'\n')




else:
    print "model not recognized - Please, re-run with model=chloroplast,mitochondrion,nucrdna"
    pass
