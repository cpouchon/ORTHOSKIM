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


parser = argparse.ArgumentParser(description='Extract closed genes for a query taxa. Script was writen by C. Pouchon (2020).')
parser.add_argument("-in","--infile", help="input reference DB (fasta format)",
                    type=str)
parser.add_argument("-q","--query", help="taxa name",
                    type=str)
parser.add_argument("-o","--outdir", help="out directory path",
                    type=str)
parser.add_argument("--target", help="gene type targeted",
                    type=str,choices=["chloroplast_CDS","chloroplast_tRNA","chloroplast_rRNA", "mitochondrion_CDS","mitochondrion_rRNA","nucrdna","nucleus_aa","nucleus_nt"])


if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)


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

def ret_stored(genenumber,genome,seeds_genome,query_fam_taxid,query,query_name):
    geneid=list(set(seeds_genome))[genenumber]
    ncbi = NCBITaxa()
    stored_families={}
    stored={}
    subgenome={k:v for k,v in genome.items() if k.startswith(str(geneid+"_"))}
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

                    if genename not in list(stored.keys()):
                        stored[genename]={}
                        stored[genename]={}
                        stored[genename][int(tax_id)]={}
                        stored[genename][int(tax_id)][taxaname]=sequence
                    else:
                        if int(tax_id) in list(stored[genename].keys()):
                            stored[genename][int(tax_id)][taxaname]=sequence
                        else:
                            stored[genename][int(tax_id)]={}
                            stored[genename][int(tax_id)][taxaname]=sequence
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

                        if genename not in list(stored.keys()):
                            stored[genename]={}
                            stored[genename]={}
                            stored[genename][int(tax_id)]={}
                            stored[genename][int(tax_id)][taxaname]=sequence
                        else:
                            if int(tax_id) in list(stored[genename].keys()):
                                stored[genename][int(tax_id)][taxaname]=sequence
                            else:
                                stored[genename][int(tax_id)]={}
                                stored[genename][int(tax_id)][taxaname]=sequence
                    else:
                        continue
        else:
            continue

    if '5.8S' in geneid:
        gene_id="rrn5.8S"
    else:
        gene_id=geneid

    # condition du query ou non pour remplacer les seeds
    # if len(ncbi.get_taxid_translator([int(query)]))==0:
    #     if len(ncbi.get_name_translator([" ".join(query_name.split("_")[0:2])]))==0:
    #         if len(ncbi.get_name_translator([" ".join(query_name.split("_")[0:1])]))==0:
    #             print("WARN: Error on query taxid: %s, the seed sequences will be considered as closest sequences for each gene for this taxa." % str(query_name))
    #             closed_taxa=list(stored[genename].keys())
    #             dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
    #             close_taxon={}
    #             value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
    #             close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
    #         else:
    #             genus_name=" ".join(query_name.split("_")[0:1])
    #             query=int(ncbi.get_name_translator([genus_name])[genus_name][0])
    #     else:
    #         species_name=" ".join(query_name.split("_")[0:2])
    #         query=int(ncbi.get_name_translator([species_name])[species_name][0])
    # else:
    #     pass

    if query!=0:
        if query_fam_taxid!=0:
            if geneid in list(stored_families.keys()):
                if query_fam_taxid in list(stored_families[geneid].keys()):
                    if query in list(stored_families[geneid][query_fam_taxid].keys()):
                        close_taxon={}
                        value_sort = sorted(stored_families[geneid][query_fam_taxid][query], key = lambda key: len(stored_families[geneid][query_fam_taxid][query][key]),reverse=True)
                        close_taxon[str(gene_id+"_"+value_sort[0])]=str(stored_families[geneid][query_fam_taxid][query][value_sort[0]])
                    else:
                        taxid_all=list(stored_families[geneid][query_fam_taxid].keys())
                        taxid_all.append(query)
                        tree = ncbi.get_topology(taxid_all)
                        # we check if the query taxid was translated or not
                        if len(tree.search_nodes(name=str(query)))==0:
                            if len(ncbi.get_lineage(query))>0:
                                search_taxid=ncbi.get_lineage(query)[len(ncbi.get_lineage(query))-1]
                                if len(tree.search_nodes(name=str(search_taxid)))>0:
                                    tcond=tree.search_nodes(name=str(search_taxid))[0].up
                                    if tcond is None:
                                        t2=tree.search_nodes(name=str(search_taxid))[0]
                                    else:
                                        t2=tree.search_nodes(name=str(search_taxid))[0].up
                                    closed_taxa = []
                                    for l in t2.get_leaves():
                                        if l.__dict__['taxid'] in taxid_all:
                                            closed_taxa.append(str(l.__dict__['taxid']))
                                        else:
                                            continue
                                    if str(query) in closed_taxa:
                                        closed_taxa.remove(str(query))
                                    elif str(search_taxid) in closed_taxa:
                                        closed_taxa.remove(str(search_taxid))
                                    else:
                                        pass

                                    if len(closed_taxa)>0:
                                        dict_close={k:v for x in closed_taxa for k,v in stored_families[geneid][query_fam_taxid][int(x)].items()}
                                        close_taxon={}
                                        value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                        close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                                    else:
                                        closed_taxa=list(stored[genename].keys())
                                        dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                        close_taxon={}
                                        value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                        close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                                else:
                                    closed_taxa=list(stored[genename].keys())
                                    dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                    close_taxon={}
                                    value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                    close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                            else:
                                #error donc prendre seeds
                                closed_taxa=list(stored[genename].keys())
                                dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                close_taxon={}
                                value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                        else:
                            tcond=tree.search_nodes(name=str(query))[0].up
                            if tcond is None:
                                t2=tree.search_nodes(name=str(query))[0]
                            else:
                                t2=tree.search_nodes(name=str(query))[0].up
                            closed_taxa = []
                            for l in t2.get_leaves():
                                if l.__dict__['taxid'] in taxid_all:
                                    closed_taxa.append(str(l.__dict__['taxid']))
                                else:
                                    continue
                            if str(query) in closed_taxa:
                                closed_taxa.remove(str(query))
                            else:
                                pass

                            if len(closed_taxa)>0:
                                dict_close={k:v for x in closed_taxa for k,v in stored_families[geneid][query_fam_taxid][int(x)].items()}
                                close_taxon={}
                                value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                            else:
                                closed_taxa=list(stored[genename].keys())
                                dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                close_taxon={}
                                value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                else:
                    #cond2 faire les close_fam

                    taxid_all_fam=list(stored_families[geneid].keys())
                    taxid_all_fam.append(query_fam_taxid)
                    tree_fam = ncbi.get_topology(taxid_all_fam)

                    if len(tree_fam.search_nodes(name=str(query_fam_taxid)))==0:
                        if len(ncbi.get_lineage(query_fam_taxid))>0:
                            search_taxid_fam=ncbi.get_lineage(query_fam_taxid)[len(ncbi.get_lineage(query_fam_taxid))-1]
                            if len(tree_fam.search_nodes(name=str(search_taxid_fam)))>0:
                                tcond=tree_fam.search_nodes(name=str(search_taxid_fam))[0].up
                                if tcond is None:
                                    t2fam=tree_fam.search_nodes(name=str(search_taxid_fam))[0]
                                else:
                                    t2fam=tree_fam.search_nodes(name=str(search_taxid_fam))[0].up
                                closed_fam = []
                                for l in t2fam.get_leaves():
                                    if l.__dict__['taxid'] in taxid_all_fam:
                                        closed_fam.append(str(l.__dict__['taxid']))
                                    else:
                                        continue
                                if str(search_taxid_fam) in closed_fam:
                                    closed_fam.remove(str(search_taxid_fam))
                                else:
                                    pass

                                if len(closed_fam)>0:
                                    taxid_all=list(y for x in closed_fam for y in list(stored_families[geneid][int(x)].keys()))
                                    taxid_all.append(query)
                                    tree = ncbi.get_topology(taxid_all)
                                    # we check if the query taxid was translated or not
                                    if len(tree.search_nodes(name=str(query)))==0:
                                        if len(ncbi.get_lineage(query))>0:
                                            search_taxid=ncbi.get_lineage(query)[len(ncbi.get_lineage(query))-1]
                                            if len(tree.search_nodes(name=str(search_taxid)))>0:
                                                tcond=tree.search_nodes(name=str(search_taxid))[0].up
                                                if tcond is None:
                                                    t2=tree.search_nodes(name=str(search_taxid))[0]
                                                else:
                                                    t2=tree.search_nodes(name=str(search_taxid))[0].up
                                                closed_taxa = []
                                                for l in t2.get_leaves():
                                                    if l.__dict__['taxid'] in taxid_all:
                                                        closed_taxa.append(str(l.__dict__['taxid']))
                                                    else:
                                                        continue
                                                if str(query) in closed_taxa:
                                                    closed_taxa.remove(str(query))
                                                elif str(search_taxid) in closed_taxa:
                                                    closed_taxa.remove(str(search_taxid))
                                                else:
                                                    pass

                                                if len(closed_taxa)>0:
                                                    dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                                    close_taxon={}
                                                    value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                                    close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                                                else:
                                                    closed_taxa=list(stored[genename].keys())
                                                    dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                                    close_taxon={}
                                                    value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                                    close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                                            else:
                                                closed_taxa=list(stored[genename].keys())
                                                dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                                close_taxon={}
                                                value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                                close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                                        else:
                                            #error donc prendre seeds
                                            closed_taxa=list(stored[genename].keys())
                                            dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                            close_taxon={}
                                            value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                            close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                                    else:
                                        tcond=tree.search_nodes(name=str(query))[0].up
                                        if tcond is None:
                                            t2=tree.search_nodes(name=str(query))[0]
                                        else:
                                            t2=tree.search_nodes(name=str(query))[0].up
                                        closed_taxa = []
                                        for l in t2.get_leaves():
                                            if l.__dict__['taxid'] in taxid_all:
                                                closed_taxa.append(str(l.__dict__['taxid']))
                                            else:
                                                continue
                                        if str(query) in closed_taxa:
                                            closed_taxa.remove(str(query))
                                        else:
                                            pass

                                        if len(closed_taxa)>0:
                                            dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                            close_taxon={}
                                            value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                            close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                                        else:
                                            closed_taxa=list(stored[genename].keys())
                                            dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                            close_taxon={}
                                            value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                            close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                                else:
                                    closed_taxa=list(stored[genename].keys())
                                    dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                    close_taxon={}
                                    value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                    close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                            else:
                                #error donc prendre seeds
                                closed_taxa=list(stored[genename].keys())
                                dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                close_taxon={}
                                value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                        else:
                            #error donc prendre seeds
                            closed_taxa=list(stored[genename].keys())
                            dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                            close_taxon={}
                            value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                            close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                    else:
                        tcond=tree_fam.search_nodes(name=str(query_fam_taxid))[0].up
                        if tcond is None:
                            t2fam=tree_fam.search_nodes(name=str(query_fam_taxid))[0]
                        else:
                            t2fam=tree_fam.search_nodes(name=str(query_fam_taxid))[0].up
                        closed_fam = []
                        for l in t2fam.get_leaves():
                            if l.__dict__['taxid'] in taxid_all_fam:
                                closed_fam.append(str(l.__dict__['taxid']))
                            else:
                                continue
                        if str(query_fam_taxid) in closed_fam:
                            closed_fam.remove(str(query_fam_taxid))
                        else:
                            pass

                        if len(closed_fam)>0:
                            taxid_all=list(y for x in closed_fam for y in list(stored_families[geneid][int(x)].keys()))
                            taxid_all.append(query)
                            tree = ncbi.get_topology(taxid_all)
                            # we check if the query taxid was translated or not
                            if len(tree.search_nodes(name=str(query)))==0:
                                if len(ncbi.get_lineage(query))>0:
                                    search_taxid=ncbi.get_lineage(query)[len(ncbi.get_lineage(query))-1]
                                    if len(tree.search_nodes(name=str(search_taxid)))>0:
                                        tcond=tree.search_nodes(name=str(search_taxid))[0].up
                                        if tcond is None:
                                            t2=tree.search_nodes(name=str(search_taxid))[0]
                                        else:
                                            t2=tree.search_nodes(name=str(search_taxid))[0].up
                                        closed_taxa = []
                                        for l in t2.get_leaves():
                                            if l.__dict__['taxid'] in taxid_all:
                                                closed_taxa.append(str(l.__dict__['taxid']))
                                            else:
                                                continue
                                        if str(query) in closed_taxa:
                                            closed_taxa.remove(str(query))
                                        elif str(search_taxid) in closed_taxa:
                                            closed_taxa.remove(str(search_taxid))
                                        else:
                                            pass

                                        if len(closed_taxa)>0:
                                            dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                            close_taxon={}
                                            value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                            close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                                        else:
                                            closed_taxa=list(stored[genename].keys())
                                            dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                            close_taxon={}
                                            value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                            close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                                    else:
                                        closed_taxa=list(stored[genename].keys())
                                        dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                        close_taxon={}
                                        value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                        close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                                else:
                                    #error donc prendre seeds
                                    closed_taxa=list(stored[genename].keys())
                                    dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                    close_taxon={}
                                    value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                    close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                            else:
                                tcond=tree.search_nodes(name=str(query))[0].up
                                if tcond is None:
                                    t2=tree.search_nodes(name=str(query))[0]
                                else:
                                    t2=tree.search_nodes(name=str(query))[0].up
                                closed_taxa = []
                                for l in t2.get_leaves():
                                    if l.__dict__['taxid'] in taxid_all:
                                        closed_taxa.append(str(l.__dict__['taxid']))
                                    else:
                                        continue
                                if str(query) in closed_taxa:
                                    closed_taxa.remove(str(query))
                                else:
                                    pass

                                if len(closed_taxa)>0:
                                    dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                    close_taxon={}
                                    value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                    close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                                else:
                                    closed_taxa=list(stored[genename].keys())
                                    dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                                    close_taxon={}
                                    value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                                    close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
                        else:
                            closed_taxa=list(stored[genename].keys())
                            dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                            close_taxon={}
                            value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                            close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
            else:
                closed_taxa=list(stored[genename].keys())
                dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
                close_taxon={}
                value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
                close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
        else:
            closed_taxa=list(stored[genename].keys())
            dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
            close_taxon={}
            value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
            close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])
    else:
        closed_taxa=list(stored[genename].keys())
        dict_close={k:v for x in closed_taxa for k,v in stored[geneid][int(x)].items()}
        close_taxon={}
        value_sort = sorted(dict_close, key = lambda key: len(dict_close[key]),reverse=True)
        close_taxon[str(gene_id+"_"+value_sort[0])]=str(dict_close[value_sort[0]])

    return close_taxon


if __name__ == "__main__":


    args = parser.parse_args()

    model=args.target
    outpath=args.outdir
    query_name = args.query
    input_file = args.infile

    if joblib.cpu_count()>2:
        num_cores=int(joblib.cpu_count()-1)
    else:
        num_cores=int(joblib.cpu_count())


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
                    print("the longest sequence for each gene was considered.")
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
                print("the longest sequence for each gene was considered.")
                query=000000

    mkdir(outpath)

    fname = "close_"+str(model)+".fa"
    open(os.path.join(outpath, fname), 'w').close()

    cur_genome=to_dict_remove_dups(SeqIO.parse(input_file, "fasta"))
    seeds_genome=[k.split("_")[0] for k in cur_genome.keys()]

    # query fam taxid
    if len(ncbi.get_taxid_translator([query]))>0:
        sub_lineages=ncbi.get_lineage(query)
        sub_names=ncbi.get_taxid_translator(sub_lineages)
        sub_ranks=ncbi.get_rank(sub_names.keys())
        list_fam=[key  for (key, value) in sub_ranks.items() if value == u'family']
        if len(list_fam)==0:
            query_fam_taxid=00000
        else:
            query_fam_taxid=list_fam[0]
    else:
        query_fam_taxid=00000

    inputs=range(len(list(set(seeds_genome))))

    #x=Parallel(n_jobs=num_cores,verbose=1)(delayed(ret_stored)(i,cur_genome) for i in inputs)
    mylist=[(x, cur_genome,seeds_genome,query_fam_taxid,query,query_name) for x in inputs]
    nprocs = num_cores
    pool = multiprocessing.Pool(nprocs)
    x=pool.starmap(ret_stored, mylist)

    close_genes={k:v for i in x for k,v in i.items()}

    # maintenant il faut ecrire le fichier
    for s in close_genes:
        out_header=">"+str(s)
        out_seq=str(close_genes[s])

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
