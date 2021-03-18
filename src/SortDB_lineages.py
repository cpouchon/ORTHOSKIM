#!/usr/bin/env python

import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
#from Bio.Alphabet import *
from Bio.SeqRecord import *
#from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA;
from Bio.Seq import *
from Bio.SeqUtils import *
from ete3 import NCBITaxa
from joblib import Parallel, delayed
import joblib
import multiprocessing

import argparse
import sys
import os, errno
import time

start_time = time.time()

parser = argparse.ArgumentParser(description='Selection of database sequences according to query rank/phylum. Script was writen by C. Pouchon (2020).')
parser.add_argument("-i","--infile", help="input genes/genomes file",
                    type=str)
parser.add_argument("-o","--outname", help="out path/name file (path is required) (fasta for genes, embl for genomes)",
                    type=str)
parser.add_argument("-p","--phylum", help="query phylum (e.g. Primulaceae, Ericales)",
                    type=str)
parser.add_argument("-r","--rank", help="query rank (e.g. family, order)",
                    type=str)


if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}


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


outname=args.outname
input_file = args.infile
qrank=args.rank
qphylum=args.phylum

open(outname, 'w').close()

cur_genome=to_dict_remove_dups(SeqIO.parse(input_file, "fasta"))

stored_families={}
ncbi = NCBITaxa()

print("rank queried: %s" % (str(qrank)))
print("phylum queried: %s" % (str(qphylum)))


seeds_genome=[k.split("_")[0] for k in cur_genome.keys()]
inputs=range(len(list(set(seeds_genome))))
Bar = ProgressBar(len(seeds_genome), 60, '\t parsing sequences')
barp=0
print("total sequence(s) found: %s" % (str(len(seeds_genome))))
print("total gene(s) found: %s" % (str(len(inputs))))
for genenumber in inputs:
    geneid=list(set(seeds_genome))[genenumber]
    subgenome={k:v for k,v in cur_genome.items() if k.startswith(str(geneid+"_"))}
    'we store all sequences'
    for record in subgenome:
        barp=barp+1
        Bar.update(barp)
        seqID=record
        sequence=subgenome[record].seq
        genename=seqID.split("_")[0]

        if genename==geneid:
            tax_id=seqID.split("_")[1]
            species="_".join(seqID.split("_")[2:len(seqID.split("_"))])
            genus=seqID.split("_")[2:len(seqID.split("_"))][0]

            tokeep={}
            tokeep['seq']=sequence
            tokeep['id']=str(tax_id)+"_"+str(species)
            taxaname=str(tax_id)+"_"+str(species)

            # faire condition pour que le taxid soit sur l'espece s'il n'existe pas
            if tax_id == "NA":
                continue
            else:
                if len(ncbi.get_taxid_translator([tax_id]))>0:
                    sub_lineages=ncbi.get_lineage(int(tax_id))
                    sub_names=ncbi.get_taxid_translator(sub_lineages)
                    sub_ranks=ncbi.get_rank(sub_names.keys())
                    list_fam=[key  for (key, value) in sub_ranks.items() if value == str(qrank)]
                    if len(list_fam)==0:
                        continue
                    else:
                        fam=ncbi.get_taxid_translator(list_fam)[list_fam[0]]
                        fam_taxid=int(list_fam[0])
                        if qphylum.lower() in fam.lower():
                            if genename not in list(stored_families.keys()):
                                stored_families[genename]={}
                                stored_families[genename][fam_taxid]={}
                                stored_families[genename][fam_taxid][int(tax_id)]={}
                                stored_families[genename][fam_taxid][int(tax_id)][taxaname]=sequence
                            else:
                                if fam_taxid not in list(stored_families[genename].keys()):
                                    stored_families[genename][fam_taxid]={}
                                    stored_families[genename][fam_taxid][int(tax_id)]={}
                                    stored_families[genename][fam_taxid][int(tax_id)][taxaname]=sequence
                                else:
                                    if int(tax_id) in list(stored_families[genename][fam_taxid].keys()):
                                        stored_families[genename][fam_taxid][int(tax_id)][taxaname]=sequence
                                    else:
                                        stored_families[genename][fam_taxid][int(tax_id)]={}
                                        stored_families[genename][fam_taxid][int(tax_id)][taxaname]=sequence
                        else:
                            continue
                else:
                    spname=" ".join(seqID.split("_")[2:4])
                    if len(ncbi.get_name_translator([spname]))>0:
                        tax_id=int(ncbi.get_name_translator([spname])[spname][0])
                        sub_lineages=ncbi.get_lineage(tax_id)
                        sub_names=ncbi.get_taxid_translator(sub_lineages)
                        sub_ranks=ncbi.get_rank(sub_names.keys())
                        list_fam=[key  for (key, value) in sub_ranks.items() if value == str(qrank)]
                        if len(list_fam)==0:
                            continue
                        else:
                            fam=ncbi.get_taxid_translator(list_fam)[list_fam[0]]
                            fam_taxid=int(list_fam[0])
                            if qphylum.lower() in fam.lower():
                                if genename not in list(stored_families.keys()):
                                    stored_families[genename]={}
                                    stored_families[genename][fam_taxid]={}
                                    stored_families[genename][fam_taxid][int(tax_id)]={}
                                    stored_families[genename][fam_taxid][int(tax_id)][taxaname]=sequence
                                else:
                                    if fam_taxid not in list(stored_families[genename].keys()):
                                        stored_families[genename][fam_taxid]={}
                                        stored_families[genename][fam_taxid][int(tax_id)]={}
                                        stored_families[genename][fam_taxid][int(tax_id)][taxaname]=sequence
                                    else:
                                        if int(tax_id) in list(stored_families[genename][fam_taxid].keys()):
                                            stored_families[genename][fam_taxid][int(tax_id)][taxaname]=sequence
                                        else:
                                            stored_families[genename][fam_taxid][int(tax_id)]={}
                                            stored_families[genename][fam_taxid][int(tax_id)][taxaname]=sequence
                            else:
                                continue
                    else:
                        continue
        else:
            continue


selected={}
taxa=[]
for g in list(stored_families.keys()):
    geneid=g
    for fam in list(stored_families[geneid].keys()):
        for sub in list(stored_families[geneid][fam].keys()):
            for all_sub in list(stored_families[geneid][fam][sub].keys()):
                fname=str(geneid+"_"+all_sub.replace(".",""))
                if all_sub in taxa:
                    pass
                else:
                    taxa.append(all_sub)
                if fname in list(selected.keys()):
                    pass
                else:
                    selected[fname]=str(stored_families[geneid][fam][sub][all_sub])

for seq in selected:
    out_header=">"+str(seq)
    out_seq=str(selected[seq])
    if os.path.isfile(outname):
        with open(outname, 'a+') as file:
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
        with open(outname, 'w') as out:
            out.write(out_header+'\n')
            out.write(str(out_seq)+'\n')

print("")
print("%s taxa found, %s sequences extracted" % (str(len(taxa)),str(len(selected))))
print("--- %s seconds ---" % (time.time() - start_time))
