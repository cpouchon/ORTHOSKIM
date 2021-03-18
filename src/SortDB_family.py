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

try:
    from Bio.Alphabet import generic_dna
except ImportError:
    generic_dna is None


import argparse
import sys
import os, errno
import time

start_time = time.time()

parser = argparse.ArgumentParser(description='Selection of database sequences by family. Script was writen by C. Pouchon (2020).')
parser.add_argument("-i","--infile", help="input genes/genomes file",
                    type=str)
parser.add_argument("-f","--format", help="input format of sequences file",choices=["embl", "genbank","fasta"],
                    type=str)
parser.add_argument("-o","--outname", help="out path/name file (path is required) (fasta for genes, embl for genomes)",
                    type=str)
parser.add_argument("-l","--lineages", help="number of lineages required by family",
                    type=int)
parser.add_argument("-m","--mode", help="selection mode",
                    type=str,choices=["gene", "genome"])


if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}



model=args.mode
outname=args.outname
input_file = args.infile
taxanumb=args.lineages
formatin=args.format


open(outname, 'w').close()

cur_genome=to_dict_remove_dups(SeqIO.parse(input_file, formatin))

stored_families={}
ncbi = NCBITaxa()

if model=="gene":
    seeds_genome=[k.split("_")[0] for k in cur_genome.keys()]
    inputs=range(len(list(set(seeds_genome))))
    for genenumber in inputs:
        geneid=list(set(seeds_genome))[genenumber]
        subgenome={k:v for k,v in cur_genome.items() if k.startswith(str(geneid+"_"))}
        'we store all sequences'
        for record in subgenome:
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
                        list_fam=[key  for (key, value) in sub_ranks.items() if value == u'family']
                        if len(list_fam)==0:
                            continue
                        else:
                            fam=ncbi.get_taxid_translator(list_fam)[list_fam[0]]
                            fam_taxid=int(list_fam[0])

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
                        spname=" ".join(seqID.split("_")[2:4])
                        if len(ncbi.get_name_translator([spname]))>0:
                            tax_id=int(ncbi.get_name_translator([spname])[spname][0])
                            sub_lineages=ncbi.get_lineage(tax_id)
                            sub_names=ncbi.get_taxid_translator(sub_lineages)
                            sub_ranks=ncbi.get_rank(sub_names.keys())
                            list_fam=[key  for (key, value) in sub_ranks.items() if value == u'family']
                            if len(list_fam)==0:
                                continue
                            else:
                                fam=ncbi.get_taxid_translator(list_fam)[list_fam[0]]
                                fam_taxid=int(list_fam[0])

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


    selected={}
    for g in list(stored_families.keys()):
        geneid=g
        select_taxon={}
        for fam in list(stored_families[geneid].keys()):
            stored_families[geneid][fam]
            dict_close={k:v for x in list(stored_families[geneid][fam].keys()) for k,v in stored_families[geneid][fam][int(x)].items()}
            genera_in=list()
            value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
            cond_lineages=0
            for v in value_sort:
                genus=v.split("_")[1]
                if genus in genera_in:
                    continue
                else:
                    if cond_lineages >= taxanumb:
                        break
                    else:
                        select_taxon[v]=dict_close[v]
                        cond_lineages=cond_lineages+1
                        genera_in.append(genus)
        for s in list(select_taxon.keys()):
            fname=str(geneid+"_"+s)
            selected[fname]=select_taxon[s]

elif model == "genome":
    for rec in cur_genome:
        record=cur_genome[rec]
        break_record=0
        for feat in record.features:
            if feat.type == 'source':
                if "/" in feat.qualifiers['organism'][0]:
                    break_record=break_record+1 #verify this process
                else:
                    org=feat.qualifiers['organism'][0].split(" ")
                    for parts in feat.qualifiers['db_xref']:
                        if "taxon" in parts:
                            taxid=parts.split(":")[1]
                    taxaname=str(taxid+"_"+"_".join(org))
            else:
                break
        if break_record==0:
            sequence=str(record.seq)

            if taxid =="NA":
                continue
            else:
                if len(ncbi.get_taxid_translator([taxid]))>0:
                    sub_lineages=ncbi.get_lineage(int(taxid))
                    sub_names=ncbi.get_taxid_translator(sub_lineages)
                    sub_ranks=ncbi.get_rank(sub_names.keys())
                    list_genus=[key  for (key, value) in sub_ranks.items() if value == u'genus']
                    if len(list_genus)==0:
                        genus=org[0]
                        if len(ncbi.get_name_translator([genus]))>0:
                            genus_taxid=int(ncbi.get_name_translator([genus])[genus][0])
                        else:
                            taxid="NA"
                            continue
                    else:
                        genus=ncbi.get_taxid_translator(list_genus)[list_genus[0]]
                        genus_taxid=list_genus[0]
                    list_fam=[key  for (key, value) in sub_ranks.items() if value == u'family']
                    if len(list_fam)==0:
                        fam="NA"
                    else:
                        fam=ncbi.get_taxid_translator(list_fam)[list_fam[0]]
                    if fam=="NA":
                        continue
                    else:
                        if fam not in stored_families.keys():
                            stored_families[fam]={}
                            stored_families[fam][taxaname]=sequence
                        else:
                            if taxaname in stored_families[fam].keys():
                                pass
                            else:
                                stored_families[fam][taxaname]=sequence
                else:
                    continue
        else:
            continue

    selected={}
    for fam in list(stored_families.keys()):
        select_taxon={}
        dict_close={k:v for k,v in stored_families[fam].items()}
        value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
        cond_lineages=0
        for v in value_sort:
            if cond_lineages >= taxanumb:
                break
            else:
                selected[v]=dict_close[v]
                cond_lineages=cond_lineages+1


if model=="gene":
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

elif model=="genome":
    for seq in selected:
        taxid=seq.split("_")[0]
        names=seq.split("_")[1:len(seq.split("_"))]
        id=" ".join(names)

        if generic_dna:
            newrecord=SeqRecord(seq=Seq(selected[seq], generic_dna), id=seq, name=id,description='',dbxrefs=[])
        else:
            newrecord=SeqRecord(seq=Seq(selected[seq]), id=id, name=seq,description='',dbxrefs=[])
        newrecord.annotations["molecule_type"] = "DNA"
        #from Bio.Alphabet import generic_dna
        #newrecord.seq.alphabet = generic_dna

        from Bio import SeqFeature
        my_start_pos = SeqFeature.ExactPosition(0)
        my_end_pos = SeqFeature.ExactPosition(len(selected[seq]))

        from Bio.SeqFeature import FeatureLocation
        my_feature_location = FeatureLocation(my_start_pos,my_end_pos)
        my_feature_type = "source"

        from Bio.SeqFeature import SeqFeature
        my_feature = SeqFeature(my_feature_location,type=my_feature_type)
        my_feature.location.strand=1
        my_feature.qualifiers['organism']=id
        my_feature.qualifiers['db_xref']=str("taxon:"+taxid)
        newrecord.features.append(my_feature)

        with open(outname, 'a+') as file:
            SeqIO.write(newrecord, file, "embl")

print("--- %s seconds ---" % (time.time() - start_time))
