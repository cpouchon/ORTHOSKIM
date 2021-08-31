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

import argparse
import sys
import os, errno
import time



parser = argparse.ArgumentParser(description='Extract closed genomes for a query taxa to assign contigs in genome type. Script was writen by C. Pouchon (2020).')
parser.add_argument("-in","--infile", help="input genomes annotation file",
                    type=str)
parser.add_argument("-fmt","--annotationfmt", help="format of annotated file",choices=["embl", "genbank"],
                    type=str)
parser.add_argument("-q","--query", help="taxa name or taxid if --taxid mode",
                    type=str)
parser.add_argument("-o","--outdir", help="out directory path",
                    type=str)
parser.add_argument("-m","--min", help="minimum size of genome",
                    type=int)
parser.add_argument("-M","--max", help="maximum size of genome",
                    type=int)
parser.add_argument("--compartment", help="cellular compartment",
                    type=str,choices=["chloroplast", "mitochondrion","nucrdna"])


if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

start_time = time.time()


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

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}



outpath=args.outdir
query_name = args.query
mol=args.compartment
formatin = args.annotationfmt
cond_size_min = args.min
cond_size_max = args.max
f=args.infile

#f="/Users/pouchonc/PhyloAlps/mitoDB/mitochondria_plants.gb"
#outpath="./test_close_ref/"
#query_name = "Pouteria_pubescens_1671090_TOU007152_BGN_GER"
#mol="mitochondrion"
#formatin = "genbank"
#cond_size_min = 200000
#cond_size_max = 1000000



res = [int(i) for i in query_name.split("_") if i.isdigit()]

ncbi = NCBITaxa()
if len(res)>0:
    if len(ncbi.get_taxid_translator([int(res[0])]))>0:
        query=res[0]
    else:
        query_genus=query_name.split("_")[0]
        query_species=" ".join(query_name.split("_")[0:2])
        if len(ncbi.get_name_translator([query_species]))>0:
            query=int(ncbi.get_name_translator([query_species])[query_species][0])
        else:
            if len(ncbi.get_name_translator([query_genus]))>0:
                query=int(ncbi.get_name_translator([query_genus])[query_genus][0])
            else:
                print("WARN: no taxid found for %s or %s" % (str(query_name),str(query_genus)))
                query=000000
else:
    query_genus=query_name.split("_")[0]
    query_species=" ".join(query_name.split("_")[0:2])
    if len(ncbi.get_name_translator([query_species]))>0:
        query=int(ncbi.get_name_translator([query_species])[query_species][0])
    else:
        if len(ncbi.get_name_translator([query_genus]))>0:
            query=int(ncbi.get_name_translator([query_genus])[query_genus][0])
        else:
            print("WARN: no taxid found for %s or %s" % (str(query_name),str(query_genus)))
            query=000000

mkdir(outpath)

fname = "close_"+str(mol)+".fa"
open(os.path.join(outpath, fname), 'w').close()



stored={}
stored_families={}
stat_size={}
families_taxid={}


file_path = f
file_name = os.path.basename(file_path)

cur_genome=to_dict_remove_dups(SeqIO.parse(f, formatin))
#cur_genome = SeqIO.parse(f, formatin)
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
                if len(sequence) >= cond_size_min and len(sequence) <= cond_size_max:
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
                        if list_fam[0] in list(families_taxid.keys()):
                            pass
                        else:
                            families_taxid[list_fam[0]]=fam
                    tokeep={}
                    tokeep[taxaname]=sequence

                    if taxid=="NA":
                        continue
                    else:
                        stat_size[taxaname]=len(sequence)
                        if taxaname not in stored.keys():
                            stored[taxaname]=sequence
                        else:
                            pass
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
    else:
        continue


if len(ncbi.get_taxid_translator([int(query)]))==0:
    if len(ncbi.get_name_translator([" ".join(query_name.split("_")[0:2])]))==0:
        if len(ncbi.get_name_translator([" ".join(query_name.split("_")[0:1])]))==0:
            sort_size=[k for k,v in sorted(stat_size.items(), key=lambda item: item[1])]
            sort_size.reverse()
            if len(sort_size)>=5:
                closed_taxa=sort_size[0:5]
            else:
                closed_taxa=sort_size
            for seq in closed_taxa:
                out_header=">"+str(seq)
                out_seq=str(stored[seq])
                if os.path.isfile(os.path.join(outpath, fname)):
                    with open(os.path.join(outpath, fname), 'a+') as file:
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
                    with open(os.path.join(outpath, fname), 'a') as out:
                        out.write(out_header+'\n')
                        out.write(str(out_seq)+'\n')
            print("WARN: Error on query taxid: %s, the five longest genomes were considered." % str(query_name))
            os._exit(0)
        else:
            genus_name=" ".join(query_name.split("_")[0:1])
            query=int(ncbi.get_name_translator([genus_name])[genus_name][0])
    else:
        species_name=" ".join(query_name.split("_")[0:2])
        query=int(ncbi.get_name_translator([species_name])[species_name][0])
else:
    pass

sub_lineages=ncbi.get_lineage(int(query))
sub_names=ncbi.get_taxid_translator(sub_lineages)
sub_ranks=ncbi.get_rank(sub_names.keys())
list_fam=[key  for (key, value) in sub_ranks.items() if value == u'family']
if len(list_fam)==0:
    sort_size=[k for k,v in sorted(stat_size.items(), key=lambda item: item[1])]
    sort_size.reverse()
    if len(sort_size)>=5:
        closed_taxa=sort_size[0:5]
    else:
        closed_taxa=sort_size
    for seq in closed_taxa:
        out_header=">"+str(seq)
        out_seq=str(stored[seq])
        if os.path.isfile(os.path.join(outpath, fname)):
            with open(os.path.join(outpath, fname), 'a+') as file:
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
            with open(os.path.join(outpath, fname), 'a') as out:
                out.write(out_header+'\n')
                out.write(str(out_seq)+'\n')
    print("WARN: Error on query taxid: %s, the five longest genomes were considered." % str(query_name))
    os._exit(0)
else:
    query_fam=ncbi.get_taxid_translator(list_fam)[list_fam[0]]
    query_fam_taxid=list_fam[0]



taxid_fam_list=[]
for i in families_taxid.keys():
    if len(ncbi.get_taxid_translator([int(i)]))==0:
        continue
    else:
        taxid_fam_list.append(int(i))

cond_loop_len_fam=len(taxid_fam_list)

if len(taxid_fam_list)>0:
    if len(taxid_fam_list)>1:
        if query_fam_taxid in taxid_fam_list:
            taxid_all_fam=taxid_fam_list
            tree_fam = ncbi.get_topology(taxid_all_fam)
            # faire comme en dessous mais ne pas supprimé le query_taxid_fam car présent
            list_loop_fam=[]
            closed_list=[]

            if len(tree_fam.search_nodes(name=str(query_fam_taxid)))==0:
                if len(ncbi.get_lineage(query_fam_taxid))>0:
                    search_taxid_fam=ncbi.get_lineage(query_fam_taxid)[len(ncbi.get_lineage(query_fam_taxid))-1]
                    if len(tree_fam.search_nodes(name=str(search_taxid_fam)))>0:
                        if len(list_loop_fam)==0:
                            t2fam=tree_fam.search_nodes(name=str(search_taxid_fam))[0]
                            closed = []
                            for l in t2fam.get_leaves():
                                if l.__dict__['taxid'] in taxid_all_fam:
                                    closed.append(str(l.__dict__['taxid']))
                                    if str(l.__dict__['taxid']) in list_loop_fam:
                                        pass
                                    else:
                                        list_loop_fam.append(str(l.__dict__['taxid']))
                                else:
                                    continue
                            closed_list.append(closed)
                        while len(list_loop_fam)<cond_loop_len_fam:
                            t2fam=t2fam.up
                            closed = []
                            for l in t2fam.get_leaves():
                                if l.__dict__['taxid'] in taxid_all_fam:
                                    closed.append(str(l.__dict__['taxid']))
                                    if str(l.__dict__['taxid']) in list_loop_fam:
                                        pass
                                    else:
                                        list_loop_fam.append(str(l.__dict__['taxid']))
                                else:
                                    continue
                            closed_list.append(closed)
                    else:
                        pass
                else:
                    pass
            else:
                if len(list_loop_fam)==0:
                    t2fam=tree_fam.search_nodes(name=str(query_fam_taxid))[0]
                    closed = []
                    for l in t2fam.get_leaves():
                        if l.__dict__['taxid'] in taxid_all_fam:
                            closed.append(str(l.__dict__['taxid']))
                            if str(l.__dict__['taxid']) in list_loop_fam:
                                pass
                            else:
                                list_loop_fam.append(str(l.__dict__['taxid']))
                        else:
                            continue
                    closed_list.append(closed)

                while len(list_loop_fam)<cond_loop_len_fam:
                    t2fam=t2fam.up
                    closed = []
                    for l in t2fam.get_leaves():
                        if l.__dict__['taxid'] in taxid_all_fam:
                            closed.append(str(l.__dict__['taxid']))
                            if str(l.__dict__['taxid']) in list_loop_fam:
                                pass
                            else:
                                list_loop_fam.append(str(l.__dict__['taxid']))
                        else:
                            continue
                    closed_list.append(closed)

        else:
            taxid_all_fam=taxid_fam_list
            taxid_all_fam.append(query_fam_taxid)
            tree_fam = ncbi.get_topology(taxid_all_fam)

            list_loop_fam=[]
            closed_list=[]

            if len(tree_fam.search_nodes(name=str(query_fam_taxid)))==0:
                if len(ncbi.get_lineage(query_fam_taxid))>0:
                    search_taxid_fam=ncbi.get_lineage(query_fam_taxid)[len(ncbi.get_lineage(query_fam_taxid))-1]
                    if len(tree_fam.search_nodes(name=str(search_taxid_fam)))>0:
                        if len(list_loop_fam)==0:
                            tcond=tree_fam.search_nodes(name=str(search_taxid_fam))[0].up
                            if tcond is None:
                                t2fam=tree_fam.search_nodes(name=str(search_taxid_fam))[0]
                            else:
                                t2fam=tree_fam.search_nodes(name=str(search_taxid_fam))[0].up
                            closed = []
                            for l in t2fam.get_leaves():
                                if l.__dict__['taxid'] in taxid_all_fam:
                                    closed.append(str(l.__dict__['taxid']))
                                    if str(l.__dict__['taxid']) in list_loop_fam:
                                        pass
                                    else:
                                        list_loop_fam.append(str(l.__dict__['taxid']))
                                else:
                                    continue
                            if str(search_taxid_fam) in closed:
                                closed.remove(str(search_taxid_fam))
                            else:
                                pass
                            if str(search_taxid_fam) in list_loop_fam:
                                list_loop_fam.remove(str(search_taxid_fam))
                            else:
                                pass
                            closed_list.append(closed)
                        while len(list_loop_fam)<cond_loop_len_fam:
                            t2fam=t2fam.up
                            closed = []
                            for l in t2fam.get_leaves():
                                if l.__dict__['taxid'] in taxid_all_fam:
                                    closed.append(str(l.__dict__['taxid']))
                                    if str(l.__dict__['taxid']) in list_loop_fam:
                                        pass
                                    else:
                                        list_loop_fam.append(str(l.__dict__['taxid']))
                                else:
                                    continue
                            if str(search_taxid_fam) in closed:
                                closed.remove(str(search_taxid_fam))
                            else:
                                pass
                            if str(search_taxid_fam) in list_loop_fam:
                                list_loop_fam.remove(str(search_taxid_fam))
                            else:
                                pass
                            closed_list.append(closed)
                    else:
                        pass
                else:
                    pass
            else:
                if len(list_loop_fam)==0:
                    tcond=tree_fam.search_nodes(name=str(query_fam_taxid))[0].up
                    if tcond is None:
                        t2fam=tree_fam.search_nodes(name=str(query_fam_taxid))[0]
                    else:
                        t2fam=tree_fam.search_nodes(name=str(query_fam_taxid))[0].up
                    closed = []
                    for l in t2fam.get_leaves():
                        if l.__dict__['taxid'] in taxid_all_fam:
                            closed.append(str(l.__dict__['taxid']))
                            if str(l.__dict__['taxid']) in list_loop_fam:
                                pass
                            else:
                                list_loop_fam.append(str(l.__dict__['taxid']))
                        else:
                            continue
                    if str(query_fam_taxid) in closed:
                        closed.remove(str(query_fam_taxid))
                    else:
                        pass
                    if str(query_fam_taxid) in list_loop_fam:
                        list_loop_fam.remove(str(query_fam_taxid))
                    else:
                        pass
                    closed_list.append(closed)

                while len(list_loop_fam)<cond_loop_len_fam:
                    t2fam=t2fam.up
                    closed = []
                    for l in t2fam.get_leaves():
                        if l.__dict__['taxid'] in taxid_all_fam:
                            closed.append(str(l.__dict__['taxid']))
                            if str(l.__dict__['taxid']) in list_loop_fam:
                                pass
                            else:
                                list_loop_fam.append(str(l.__dict__['taxid']))
                        else:
                            continue
                    if str(query_fam_taxid) in closed:
                        closed.remove(str(query_fam_taxid))
                    else:
                        pass
                    if str(query_fam_taxid) in list_loop_fam:
                        list_loop_fam.remove(str(query_fam_taxid))
                    else:
                        pass
                    closed_list.append(closed)

    else:
        closed_list=[]
        closed=[]
        closed.append(str(taxid_fam_list[0]))
        closed_list.append(closed)

    cond_len = 1
    closed_taxa=[]
    if len(closed_list)>0:
        # on choisit par la taille dans chaque close fam
        for i in closed_list:
            seqs=list()
            for j in i:
                sublist=list(stored_families[families_taxid[int(j)]].keys())
                for y in sublist:
                    seqs.append(y)
            lenseq=[stat_size[x] for x in seqs]
            sortseqs = [x for _,x in sorted(zip(lenseq,seqs))]
            sortseqs.reverse()

            for s in sortseqs:
                if cond_len == 6:
                    break
                else:
                    if s in closed_taxa:
                        pass
                    else:
                        closed_taxa.append(s)
                        cond_len += 1

else:
    closed_taxa=[]
    # pas de close car pas de sequences de familles donc --> à la fin prendre sur la taille si >5 --> 5 sinon tout


if len(closed_taxa)==0:
    # prendre 5 plus grand si plus de 5 sequences sinon tout prendre
    print("Error on selection of closed genomes for query taxa: %s; the five longest genomes were considered." % str(query_name))
    sort_size=[k for k,v in sorted(stat_size.items(), key=lambda item: item[1])]
    sort_size.reverse()
    if len(sort_size)>=5:
        closed_taxa=sort_size[0:5]
    else:
        closed_taxa=sort_size
    for seq in closed_taxa:
        out_header=">"+str(seq)
        out_seq=str(stored[seq])
        if os.path.isfile(os.path.join(outpath, fname)):
            with open(os.path.join(outpath, fname), 'a+') as file:
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
            with open(os.path.join(outpath, fname), 'a') as out:
                out.write(out_header+'\n')
                out.write(str(out_seq)+'\n')
else:
    if len(closed_taxa)>=5:
        for seq in closed_taxa:
            out_header=">"+str(seq)
            out_seq=str(stored[seq])
            if os.path.isfile(os.path.join(outpath, fname)):
                with open(os.path.join(outpath, fname), 'a+') as file:
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
                with open(os.path.join(outpath, fname), 'a') as out:
                    out.write(out_header+'\n')
                    out.write(str(out_seq)+'\n')
    else:
        #prendre les plus grand si taille >5 sinon tout prendre en complément si taxon pas deja present
        n_to_add=5-len(closed_taxa)
        sort_size=[k for k,v in sorted(stat_size.items(), key=lambda item: item[1])]
        sort_size.reverse()
        if len(sort_size)>=5:
            cond_len=1
            for tax in sort_size:
                if cond_len == 6:
                    break
                else:
                    if tax in closed_taxa:
                        pass
                    else:
                        closed_taxa.append(tax)
                        cond_len+=1
        else:
            for tax in sort_size:
                if tax in closed_taxa:
                    pass
                else:
                    closed_taxa.append(tax)
        for seq in closed_taxa:
            out_header=">"+str(seq)
            out_seq=str(stored[seq])
            if os.path.isfile(os.path.join(outpath, fname)):
                with open(os.path.join(outpath, fname), 'a+') as file:
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
                with open(os.path.join(outpath, fname), 'a') as out:
                    out.write(out_header+'\n')
                    out.write(str(out_seq)+'\n')

print("--- %s seconds ---" % (time.time() - start_time))
