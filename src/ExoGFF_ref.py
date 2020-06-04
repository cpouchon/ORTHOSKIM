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
from Bio.Alphabet import generic_dna


parser = argparse.ArgumentParser(description='Export regions of query database from exonerate alignment. Filtering is made according to rawscore. Script was writen by C. Pouchon (2019).')
parser.add_argument("-i","--infile", help="Fasta input contigs file",
                    type=str)
parser.add_argument("-m","--model", help="molecular type compartment",
                    type=str,choices=["mitochondrion_CDS","mitochondrion_rRNA","mitochondrion_tRNA","chloroplast_CDS","chloroplast_RNA","nucrdna_rRNA"])
parser.add_argument("-o","--outdir", help="Out directory path",
                    type=str)
parser.add_argument("-g","--gfftab", help="input gff table from exonerate alignment",
                    type=str)
parser.add_argument("-s","--seeds", help="input reference seeds",
                    type=str)
parser.add_argument("--threads", help="number of threads to use",
                    type=int)
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
    contpart=[]
    geneid=list(besthits.keys())[genenumber]

    if 'rrnITS' not in geneid:
        if len(besthits[geneid])==1:
            seqid = besthits[geneid][0][1]
            dna=str(seqs[seqid][0])
            newid = str(geneid)+"_"+"_".join(seqid.split("_")[1:len(seqid.split("_"))])
        else:
            hitparts={}
            hitlen={}
            for hits in besthits[geneid]:
                if hits[1] in hitparts.keys():
                    hitparts[hits[1]]=hitparts[hits[1]]+1
                else:
                    hitparts[hits[1]]=1
                if hits[1] in hitlen.keys():
                    lenhit=int(hits[4])-int(hits[3])
                    hitlen[hits[1]]=hitlen[hits[1]]+lenhit
                else:
                    lenhit=int(hits[4])-int(hits[3])
                    hitlen[hits[1]]=lenhit

            if len(hitparts.keys()) ==1:
                'sequence unique donc on recupere le bon id'
                seqid = besthits[geneid][0][1]
                dna=str(seqs[seqid][0])
                newid = str(geneid)+"_"+"_".join(seqid.split("_")[1:len(seqid.split("_"))])

            else:
                lmax=list()
                for uniqhit in hitparts.keys():
                    lmax.append(hitlen[uniqhit])
                orderindx=sorted(range(len(lmax)), key=lambda k: lmax[k])
                orderindx.reverse()
                'on garde le plus grand en alignement (reverse du min)'
                seqid = hitparts.keys()[orderindx[0]]
                dna=str(seqs[seqid][0])
                newid = str(geneid)+"_"+"_".join(seqid.split("_")[1:len(seqid.split("_"))])

        header = ">"+str(newid)
        taxa="_".join(seqid.split("_")[1:len(seqid.split("_"))])

        if os.path.isfile(os.path.join(outpath+"/"+model, taxa+".fa")):
            with open(os.path.join(outpath+"/"+model, taxa+".fa"), 'a+') as file:
                old_headers = []
                end_file=file.tell()
                file.seek(0)
                for line in file:
                    if line.startswith(">"):
                        old_headers.append(line.rstrip().replace(">",""))
                if not str(newid) in old_headers:
                    file.seek(end_file)
                    file.write(header+'\n')
                    file.write(str(dna)+'\n')
                else:
                    pass
        else :
            with open(os.path.join(outpath+"/"+model, taxa+".fa"), 'w') as out:
                out.write(header+'\n')
                out.write(str(dna)+'\n')
    else:
        pass




file=args.gfftab
input_file = args.infile
model = args.model
outpath=args.outdir
num_cores = args.threads
tabseeds=args.seeds

stored={}
dicscore={}
besthits={}
groups={}
seqs={}
counthits={}

if model =="nucrdna_rRNA":
    mkdir(str(outpath+"/"+model))
    mkdir(str(outpath+"/"+model.replace("_rRNA","_misc_RNA")))
else:
    mkdir(str(outpath+"/"+model))

'we store the gff tab'
with open(file) as f:
    tab = f.readlines()

with open(tabseeds) as f:
    reftab = f.readlines()

refnames=list()
for line in reftab:
    if line.startswith(">"):
        l=line.rstrip()
        genename = l.replace(">","")
        if genename not in refnames:
            refnames.append(genename)
        else:
            pass

'we store all sequences'
cur_genome = SeqIO.parse(input_file, "fasta")
for record in cur_genome:
    seqID=record.id
    sequence=record.seq
    seqs.setdefault(seqID, []).append(sequence)


'find each gene of ref and keep seqID,score,length,ref'
for line in tab:
    l = line.rstrip().split("\t")
    'intialisation of the scan'
    if l[0] in refnames:

        genename = l[0].split("_")[0]

        refid = l[0]
        seqid = l[8].split("; ")[1].split(" ")[1]
        rmin = int(l[3])
        rmax = int(l[4])
        rscore = int(l[5])

        refposmin=rmin
        refposmax=rmax



        'we make groups for gene with borne values'
        if genename not in groups.keys():
            dic=[]
            dic.append(rmin)
            dic.append(rmax)
            groups.setdefault(genename, []).append(dic)
        else:
            'if existed, we scan each subgrp to check if the hits overlapp others hits on reference sequence'
            loverlap=list()
            lnoverlap=list()

            for i in range(len(groups[genename])):
                m = min(groups[genename][i])
                M = max(groups[genename][i])
                'we create all condition of overlapping'
                if rmin >= m:
                    if M >= rmin:
                        if M <= rmax:
                            '0'
                            rmin = m
                            loverlap.append(i)
                            continue
                        else:
                            '4'
                            rmin=m
                            rmax=M
                            loverlap.append(i)
                            continue
                    else:
                        lnoverlap.append(i)
                        continue

                if rmin <= m:
                    if rmax > m:
                        if M >= rmin:
                            if M >= rmax:
                                '2'
                                rmax=M
                                loverlap.append(i)
                                continue
                            else:
                                '1'
                                loverlap.append(i)
                                continue
                    else:
                        lnoverlap.append(i)
                        continue

            'on modifie le grp pour la partie overlapp qui remplacera les grp existent / ou on fait un nouveau grp'
            overgp=[]
            overgp.append(rmin)
            overgp.append(rmax)

            novergp=[]
            for i in lnoverlap:
                novergp.append(groups[genename][i])
            novergp.append(overgp)
            groups[genename]=novergp


        contname=[]
        contname.append(seqid)
        if genename not in counthits.keys():
             counthits.setdefault(genename, []).append(contname)
        else:
            newcont=[]

            if len(lnoverlap)>0:
                for i in lnoverlap:
                     newcont.append(counthits[genename][i])

            overcont=[]
            'il faut modifier pour avoir que les element est pas une liste dans une liste'
            if len(loverlap)>0:
                for i in loverlap:
                    overcont.extend(counthits[genename][i])
                if seqid not in overcont:
                    overcont.extend(contname)
                newcont.append(overcont)
            else:
                newcont.append(contname)

            counthits[genename]=newcont

    else:
        if l[2]=="similarity":

            score = l[5]
            frame = l[6]
            id=(genename,seqid,refid,refposmin,refposmax)
            if genename not in dicscore.keys():
                dicscore.setdefault(genename, []).append(int(score))
                #besthits.setdefault(genename, []).append(dicinfo)
                besthits.setdefault(genename, []).append(id)
                #counthits.setdefault(genename, []).append([seqid])
            else:
                'mettre condition pour regarder dans loverlapp index'

                'condition si nouveau group'

                newscore=[]
                newhits=[]
                overscore=[]

                for i in lnoverlap:
                    'ceux pour qui pas overlapp on les gardes dans le meme ordre (new/modif in first)'
                    newscore.append(dicscore[genename][i])
                    newhits.append(besthits[genename][i])

                if len(loverlap)>0:
                    for i in loverlap:
                        overscore.append(dicscore[genename][i])
                    'extract best score on new overlapping group'
                    highscore_index = loverlap[overscore.index(max(overscore))]
                    highscore = dicscore[genename][highscore_index]
                    'compare hit score with the bestscore for this overlapping group'
                    if int(score) > highscore:
                        highscore = int(score)
                        #highhits = dicinfo
                        highhits = id
                    else:
                        highhits = besthits[genename][highscore_index]

                    'maintenant faut faire le nouveau gp'
                    newscore.append(highscore)
                    newhits.append(highhits)
                else:
                    'pas overlapp donc new'
                    newscore.append(int(score))
                    #newhits.append(dicinfo)
                    newhits.append(id)


                dicscore[genename]=newscore
                besthits[genename]=newhits


            'we store all combination genename/seqid/refid + frame'
            #id = (genename,seqid,refid,refposmin,refposmax)
            dicframe=dict()
            dicframe['frame']=str(frame)
            stored.setdefault(id, []).append(dicframe)

            dicinfo = dict()
            minpos = l[3]
            maxpos = l[4]
            dicinfo['min']=str(minpos)
            dicinfo['max']=str(maxpos)
            stored.setdefault(id, []).append(dicinfo)


print ("clean sequences for %s" % input_file)
inputs=range(len(list(besthits.keys())))
Parallel(n_jobs=num_cores)(delayed(GeneExtraction)(i) for i in inputs)
