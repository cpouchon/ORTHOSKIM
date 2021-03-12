#!/usr/bin/env python

import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import sys
import os, errno
import argparse
from joblib import Parallel, delayed
import multiprocessing
import time

from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna
start_time = time.time()

parser = argparse.ArgumentParser(description='Export regions of query database from exonerate alignment. Filtering is made according to rawscore. Script was writen by C. Pouchon (2019).')
parser.add_argument("-i","--infile", help="Fasta input contigs file",
                    type=str)
parser.add_argument("-m","--model", help="molecular type compartment",
                    type=str,choices=["chloroplast", "mitochondrion_CDS","mitochondrion_rRNA","mitochondrion_tRNA","nucleus","chloroplast_CDS","chloroplast_rRNA","chloroplast_tRNA","nucrdna"])
parser.add_argument("-o","--outdir", help="Out directory path",
                    type=str)
parser.add_argument("-n","--namesample", help="sample name used in sequence header for output",
                    type=str)
parser.add_argument("-g","--gfftab", help="input gff table from exonerate alignment",
                    type=str)
parser.add_argument("-t","--typeofseq", help="type of sequence to extract",
                    type=str,choices=["all", "exon","intron"])
parser.add_argument("-l","--minlength", help="minimal length of alignment into reference to consider output",
                    type=int)
parser.add_argument("-c","--mincov", help="minimal coverage of contigs allowed for genomic scan",
                    type=int)
parser.add_argument("-cl","--mincontlen", help="minimal length of contifs allowed for genomic scan",
                    type=int)
parser.add_argument("-rp","--refpct", help="minimal percent of reference covered allowed",
                    type=float)
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

def GeneExtraction(genenumber):
    contpart=[]
    num_exon=[]
    num_intron=[]
    geneid=list(besthits_filtered2.keys())[genenumber]
    extract_contigs=list()
    if "rrnITS" in geneid:
        if len(besthits_filtered2[geneid])==1:
            contpart.append(str(len(set(counthits[geneid][0]))))
            contigcount=len(set(counthits[geneid][0]))

            seqid = besthits_filtered2[geneid][0][1]
            extract_contigs.append(seqid)
            #fr = stored[besthits[geneid][0]][0]['frame']
            dna=str(seqs[seqid][0])
            concat=[]
            list_combo=list()
            pos_to_extract=list()
            exon_ref_cov=0
            for hits in besthits_filtered2[geneid]:
                exon_ref_cov=exon_ref_cov+(hits[4]-hits[3]+1)
            gene_percent=round(float(exon_ref_cov)/float(reflength[geneid]),2)
            if len(stored[besthits_filtered2[geneid][0]])==0:
                pass
            else:
                for parts in stored[besthits_filtered2[geneid][0]]:
                    if 'frame' in parts.keys():
                        fr=parts['frame']
                    else:
                        start=int(parts['min'])
                        end=int(parts['max'])
                        combo=str(start)+"_"+str(end)+"_"+fr

                        if len(list_combo)==0:
                            list_combo.append(combo)
                            pos_to_extract.append(combo)
                        else:
                            num_overlapp=0
                            for elem in list_combo:
                                m=int(elem.split("_")[0])
                                M=int(elem.split("_")[1])
                                if start <= m:
                                    if end >= m:
                                        num_overlapp=num_overlapp+1
                                    elif end == M:
                                        num_overlapp=num_overlapp+1
                                    else:
                                        pass
                                else:
                                    if start <= M:
                                        num_overlapp=num_overlapp+1
                                    else:
                                        pass
                            if num_overlapp==0:
                                pos_to_extract.append(combo)
                                list_combo=pos_to_extract
                            else:
                                pass

                sublmin=list()
                framelist=list()
                for i in pos_to_extract:
                    sublmin.append(int(i.split("_")[0]))
                    if i.split("_")[2] in framelist:
                        pass
                    else:
                        framelist.append(i.split("_")[2])
                suborderindx=sorted(range(len(sublmin)), key=lambda k: sublmin[k])

                if len(framelist)==1:
                    if framelist[0]=="-":
                        suborderindx.reverse()
                    else:
                        pass
                else:
                    pass

                for i in suborderindx:
                    posmin=int(pos_to_extract[i].split("_")[0])
                    posmax=int(pos_to_extract[i].split("_")[1])
                    posfr=pos_to_extract[i].split("_")[2]
                    if posfr == "-":
                        extract=str(Seq(dna[posmin-1:posmax]).reverse_complement())
                    else:
                        extract=str(Seq(dna[posmin-1:posmax]))
                    concat.append(extract)

                #newseq="".join(concat)
        else:
            'condition pour voir quelle position il faut extraire en premier'
            lmin=list()

            exon_ref_cov=0
            for hits in besthits_filtered2[geneid]:
                lmin.append(hits[3])
                exon_ref_cov=exon_ref_cov+(hits[4]-hits[3]+1)
            orderindx=sorted(range(len(lmin)), key=lambda k: lmin[k])
            gene_percent=round(float(exon_ref_cov)/float(reflength[geneid]),2)
            concat=[]
            for idx in orderindx:
                contpart.append(str(len(set(counthits[geneid][idx]))))
                seqid = besthits_filtered2[geneid][idx][1]
                extract_contigs.append(seqid)
                #fr = stored[besthits[geneid][idx]][0]['frame']
                dna=str(seqs[seqid][0])
                list_combo=list()
                pos_to_extract=list()
                if len(stored[besthits_filtered2[geneid][idx]])==0:
                    pass
                else:
                    for parts in stored[besthits_filtered2[geneid][idx]]:
                        if 'frame' in parts.keys():
                            fr=parts['frame']
                        else:
                            start=int(parts['min'])
                            end=int(parts['max'])
                            combo=str(start)+"_"+str(end)+"_"+fr

                            if len(list_combo)==0:
                                list_combo.append(combo)
                                pos_to_extract.append(combo)
                            else:
                                num_overlapp=0
                                for elem in list_combo:
                                    m=int(elem.split("_")[0])
                                    M=int(elem.split("_")[1])
                                    if start <= m:
                                        if end >= m:
                                            num_overlapp=num_overlapp+1
                                        elif end == M:
                                            num_overlapp=num_overlapp+1
                                        else:
                                            pass
                                    else:
                                        if start <= M:
                                            num_overlapp=num_overlapp+1
                                        else:
                                            pass
                                if num_overlapp==0:
                                    pos_to_extract.append(combo)
                                    list_combo=pos_to_extract
                                else:
                                    pass

                    sublmin=list()
                    framelist=list()
                    for i in pos_to_extract:
                        sublmin.append(int(i.split("_")[0]))
                        if i.split("_")[2] in framelist:
                            pass
                        else:
                            framelist.append(i.split("_")[2])
                    suborderindx=sorted(range(len(sublmin)), key=lambda k: sublmin[k])

                    if len(framelist)==1:
                        if framelist[0]=="-":
                            suborderindx.reverse()
                        else:
                            pass
                    else:
                        pass

                    for i in suborderindx:
                        posmin=int(pos_to_extract[i].split("_")[0])
                        posmax=int(pos_to_extract[i].split("_")[1])
                        posfr=pos_to_extract[i].split("_")[2]
                        if posfr == "-":
                            extract=str(Seq(dna[posmin-1:posmax]).reverse_complement())
                        else:
                            extract=str(Seq(dna[posmin-1:posmax]))
                        concat.append(extract)
    else:
        if len(besthits_filtered2[geneid])==1:
            refcoverage_dict={}
            refcoverage_pass=[]

            seqid = besthits_filtered2[geneid][0][1]
            concat=[]

            for hits in besthits_filtered2[geneid]:
                seqid=hits[1]
                if seqid in refcoverage_dict:
                    exon_cov=hits[4]-hits[3]+1
                    refcoverage_dict[seqid]=refcoverage_dict[seqid]+int(exon_cov)
                else:
                    refcoverage_dict[seqid]=hits[4]-hits[3]+1

            for cont in refcoverage_dict.keys():
                genepercent=round(float(refcoverage_dict[cont])/float(reflength[geneid]),2)
                if genepercent>=refcond:
                    refcoverage_pass.append(cont)
                else:
                    pass

            if seqid in refcoverage_pass:
                contpart.append(str(len(set(counthits[geneid][0]))))
                contigcount=len(set(counthits[geneid][0]))
                extract_contigs.append(seqid)
                #fr = stored[besthits[geneid][0]][0]['frame']
                dna=str(seqs[seqid][0])

                list_combo=list()
                pos_to_extract=list()

                if len(stored[besthits_filtered2[geneid][0]])==0:
                    pass
                else:
                    for parts in stored[besthits_filtered2[geneid][0]]:
                        if 'frame' in parts.keys():
                            fr=parts['frame']
                        else:
                            start=int(parts['min'])
                            end=int(parts['max'])
                            combo=str(start)+"_"+str(end)+"_"+fr

                            if len(list_combo)==0:
                                list_combo.append(combo)
                                pos_to_extract.append(combo)
                            else:
                                num_overlapp=0
                                for elem in list_combo:
                                    m=int(elem.split("_")[0])
                                    M=int(elem.split("_")[1])
                                    if start <= m:
                                        if end >= m:
                                            num_overlapp=num_overlapp+1
                                        elif end == M:
                                            num_overlapp=num_overlapp+1
                                        else:
                                            pass
                                    else:
                                        if start <= M:
                                            num_overlapp=num_overlapp+1
                                        else:
                                            pass
                                if num_overlapp==0:
                                    pos_to_extract.append(combo)
                                    list_combo=pos_to_extract
                                else:
                                    pass

                    sublmin=list()
                    framelist=list()
                    for i in pos_to_extract:
                        sublmin.append(int(i.split("_")[0]))
                        if i.split("_")[2] in framelist:
                            pass
                        else:
                            framelist.append(i.split("_")[2])
                    suborderindx=sorted(range(len(sublmin)), key=lambda k: sublmin[k])

                    if len(framelist)==1:
                        if framelist[0]=="-":
                            suborderindx.reverse()
                        else:
                            pass
                    else:
                        pass

                    for i in suborderindx:
                        posmin=int(pos_to_extract[i].split("_")[0])
                        posmax=int(pos_to_extract[i].split("_")[1])
                        posfr=pos_to_extract[i].split("_")[2]
                        if posfr == "-":
                            extract=str(Seq(dna[posmin-1:posmax]).reverse_complement())
                        else:
                            extract=str(Seq(dna[posmin-1:posmax]))
                        concat.append(extract)

                    #newseq="".join(concat)

                    num_exon.append(str(stat_stored[besthits_filtered2[geneid][0]]["exon"]))
                    num_intron.append(str(stat_stored[besthits_filtered2[geneid][0]]["intron"]))

        else:
            'condition pour voir quelle position il faut extraire en premier'
            lmin=list()

            refcoverage_dict={}
            refcoverage_pass=[]

            for hits in besthits_filtered2[geneid]:
                seqid=hits[1]
                lmin.append(hits[3])
                if seqid in refcoverage_dict:
                    exon_cov=hits[4]-hits[3]+1
                    refcoverage_dict[seqid]=refcoverage_dict[seqid]+int(exon_cov)
                else:
                    refcoverage_dict[seqid]=hits[4]-hits[3]+1

            for cont in refcoverage_dict.keys():
                genepercent=round(float(refcoverage_dict[cont])/float(reflength[geneid]),2)
                if genepercent>=refcond:
                    refcoverage_pass.append(cont)
                else:
                    pass

            orderindx=sorted(range(len(lmin)), key=lambda k: lmin[k])
            concat=[]
            for idx in orderindx:
                seqid = besthits_filtered2[geneid][idx][1]
                if seqid in refcoverage_pass:
                    contpart.append(str(len(set(counthits[geneid][idx]))))
                    extract_contigs.append(seqid)
                    #fr = stored[besthits[geneid][idx]][0]['frame']
                    dna=str(seqs[seqid][0])
                    list_combo=list()
                    pos_to_extract=list()
                    if len(stored[besthits_filtered2[geneid][idx]])==0:
                        pass
                    else:
                        for parts in stored[besthits_filtered2[geneid][idx]]:
                            if 'frame' in parts.keys():
                                fr=parts['frame']
                            else:
                                start=int(parts['min'])
                                end=int(parts['max'])
                                combo=str(start)+"_"+str(end)+"_"+fr

                                if len(list_combo)==0:
                                    list_combo.append(combo)
                                    pos_to_extract.append(combo)
                                else:
                                    num_overlapp=0
                                    for elem in list_combo:
                                        m=int(elem.split("_")[0])
                                        M=int(elem.split("_")[1])
                                        if start <= m:
                                            if end >= m:
                                                num_overlapp=num_overlapp+1
                                            elif end == M:
                                                num_overlapp=num_overlapp+1
                                            else:
                                                pass
                                        else:
                                            if start <= M:
                                                num_overlapp=num_overlapp+1
                                            else:
                                                pass
                                    if num_overlapp==0:
                                        pos_to_extract.append(combo)
                                        list_combo=pos_to_extract
                                    else:
                                        pass

                        sublmin=list()
                        framelist=list()
                        for i in pos_to_extract:
                            sublmin.append(int(i.split("_")[0]))
                            if i.split("_")[2] in framelist:
                                pass
                            else:
                                framelist.append(i.split("_")[2])
                        suborderindx=sorted(range(len(sublmin)), key=lambda k: sublmin[k])

                        if len(framelist)==1:
                            if framelist[0]=="-":
                                suborderindx.reverse()
                            else:
                                pass
                        else:
                            pass

                        for i in suborderindx:
                            posmin=int(pos_to_extract[i].split("_")[0])
                            posmax=int(pos_to_extract[i].split("_")[1])
                            posfr=pos_to_extract[i].split("_")[2]
                            if posfr == "-":
                                extract=str(Seq(dna[posmin-1:posmax]).reverse_complement())
                            else:
                                extract=str(Seq(dna[posmin-1:posmax]))
                            concat.append(extract)

                    num_exon.append(str(stat_stored[besthits_filtered2[geneid][idx]]["exon"]))
                    num_intron.append(str(stat_stored[besthits_filtered2[geneid][idx]]["intron"]))
                else:
                    continue

    newseq="".join(concat).replace("N","")
    newlength = len(newseq)

    if 'rrnITS' in geneid:
        cond_rrna_contigs=0
        for cont in extract_contigs:
            if cont in rrna_contigs:
                cond_rrna_contigs= cond_rrna_contigs+1
            else:
                pass
        if cond_rrna_contigs>0:
            header = ">"+str(nameofsample)+"; "+"gene="+str(geneid.replace("rrn",""))+"; "+"info=intergenic spacer"+"; "+"type="+str(model)+"; "+"length="+str(newlength)+"; "+"match_contigs="+str("-".join(contpart))+"; "+"ref_percent="+str(gene_percent)
            fname = str(geneid.replace("rrn",""))+".fa"
            cond_min_length=0
            if newlength>=min_length:
                cond_min_length=cond_min_length+1
            else:
                pass
            if cond_min_length!=0:
                if os.path.isfile(os.path.join(outpath+"/"+model, fname)):
                    with open(os.path.join(outpath+"/"+model, fname), 'a+') as file:
                        old_headers = []
                        end_file=file.tell()
                        file.seek(0)
                        for line in file:
                            if line.startswith(">"):
                                old_headers.append(line.rstrip().replace(">","").split(";")[0])
                        if not nameofsample in old_headers:
                            file.seek(end_file)
                            file.write(header+'\n')
                            file.write(str(newseq)+'\n')
                        else:
                            pass
                else :
                    with open(os.path.join(outpath+"/"+model, fname), 'w') as out:
                        out.write(header+'\n')
                        out.write(str(newseq)+'\n')
                for elem in extract_contigs:
                    return elem
            else:
                pass
        else:
            pass
    else:
        if typeseq=="exon":
            genepart=typeseq
        elif typeseq=="intron":
            genepart=typeseq
        elif typeseq=="all":
            genepart="exon+intron"

        exon_l_recovered=0
        for p in refcoverage_pass:
            exon_l_recovered=exon_l_recovered+int(refcoverage_dict[p])
        gpct=round(float(exon_l_recovered)/float(reflength[geneid]),2)

        fname = geneid+".fa"
        cond_min_length=0
        if newlength>=min_length:
            cond_min_length=cond_min_length+1
        else:
            pass

        if geneid=="rrn5.8S":
            cond_rrna_contigs=0
            for cont in extract_contigs:
                if cont in rrna_contigs:
                    cond_rrna_contigs= cond_rrna_contigs+1
                else:
                    pass
            if cond_rrna_contigs>0:
                if cond_min_length!=0:
                    if len(contpart)>1 or len(num_exon)>1 or len(num_intron)>1:
                        pass
                    else:
                        header = ">"+str(nameofsample)+"; "+"gene="+str(geneid)+"; "+"info="+str(genepart)+"; "+"type="+str(model)+"; "+"length="+str(newlength)+"; "+"match_contigs="+str("-".join(contpart))+"; "+"ref_percent="+str(gpct)+"; "+"n_exons="+str("-".join(num_exon))+"; "+"n_introns="+str("-".join(num_intron))
                        if os.path.isfile(os.path.join(outpath+"/"+model, fname)):
                            with open(os.path.join(outpath+"/"+model, fname), 'a+') as file:
                                old_headers = []
                                end_file=file.tell()
                                file.seek(0)
                                for line in file:
                                    if line.startswith(">"):
                                        old_headers.append(line.replace(">","").split(";")[0])
                                if not nameofsample in old_headers:
                                    file.seek(end_file)
                                    file.write(header+'\n')
                                    file.write(str(newseq)+'\n')
                                else:
                                    pass
                        else :
                            with open(os.path.join(outpath+"/"+model, fname), 'w') as out:
                                out.write(header+'\n')
                                out.write(str(newseq)+'\n')
                        for elem in extract_contigs:
                            return elem
                else:
                    pass
            else:
                pass
        else:
            if cond_min_length!=0:
                if len(contpart)>1 or len(num_exon)>1 or len(num_intron)>1:
                    pass
                else:
                    header = ">"+str(nameofsample)+"; "+"gene="+str(geneid)+"; "+"info="+str(genepart)+"; "+"type="+str(model)+"; "+"length="+str(newlength)+"; "+"match_contigs="+str("-".join(contpart))+"; "+"ref_percent="+str(gpct)+"; "+"n_exons="+str("-".join(num_exon))+"; "+"n_introns="+str("-".join(num_intron))
                    if os.path.isfile(os.path.join(outpath+"/"+model, fname)):
                        with open(os.path.join(outpath+"/"+model, fname), 'a+') as file:
                            old_headers = []
                            end_file=file.tell()
                            file.seek(0)
                            for line in file:
                                if line.startswith(">"):
                                    old_headers.append(line.replace(">","").split(";")[0])
                            if not nameofsample in old_headers:
                                file.seek(end_file)
                                file.write(header+'\n')
                                file.write(str(newseq)+'\n')
                            else:
                                pass
                    else :
                        with open(os.path.join(outpath+"/"+model, fname), 'w') as out:
                            out.write(header+'\n')
                            out.write(str(newseq)+'\n')
                    for elem in extract_contigs:
                        return elem
            else:
                pass

file=args.gfftab
typeseq = args.typeofseq
input_file = args.infile
nameofsample = args.namesample
min_length = args.minlength
model = args.model
outpath=args.outdir
cov=args.mincov
clen=args.mincontlen
num_cores = args.threads
tabseeds=args.seeds
refcond=args.refpct

stored={}
stat_stored={}
dicscore={}
besthits={}
groups={}
seqs={}
counthits={}
dict_aln={}


mkdir(str(outpath+"/"+model))

'we store the gff tab'
with open(file) as f:
    tab = f.readlines()

if len(tab)<2:
    print(str("WARN: "+model+" - no contigs mapped on references for "+nameofsample))
else:
    pass

with open(tabseeds) as f:
    reftab = f.readlines()

refnames=list()
reflength=dict()

for line in reftab:
    if line.startswith(">"):
        l=line.rstrip()
        genename = l.replace(">","")
        if genename not in refnames:
            refnames.append(genename)
        else:
            pass
        if "UCE" in l:
            genename = l.replace(">","").split("_")[0]+"_UCE"
        elif "BUSCO" in l:
            genename = l.replace(">","").split("_")[0]+"_BUSCO"
        else:
            genename = l.replace(">","").split("_")[0]
    else:
        l=line.rstrip()
        reflength[genename]=len(l)

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

        if "UCE" in l[0]:
            genename = l[0].split("_")[0]+"_UCE"
        elif "BUSCO" in l[0]:
            genename = l[0].split("_")[0]+"_BUSCO"
        else:
            genename = l[0].split("_")[0]

        refid = l[0]
        seqid = l[8].split("; ")[1].split(" ")[1]
        rmin = int(l[3])
        rmax = int(l[4])
        rscore = int(l[5])

        num_exon=0
        num_intron=0

        refposmin=rmin
        refposmax=rmax

        aln_info=(str(genename),str(seqid),str(refid),str(refposmin),str(refposmax),str(rscore))
        dict_aln.setdefault(genename, []).append(aln_info)

        contiglength=int(seqid.split("_")[3])
        contigcov=float(seqid.split("_")[5])
    else:
        if 'ITS' in genename:
            if l[2]=="gene":
                if contigcov<cov or contiglength<clen:
                    continue
                else:
                    pass
                score = l[5]
                frame = l[6]
                id=(genename,seqid,refid,refposmin,refposmax,score,frame)

                if id in stat_stored.keys():
                    pass
                else:
                    stat_stored[id]={}
                    stat_stored[id]["exon"]=0
                    stat_stored[id]["intron"]=0

            else:
                if l[2]=='exon':
                    if contigcov<cov or contiglength<clen:
                        continue
                    else:
                        pass
                    old_num=stat_stored[id]["exon"]
                    stat_stored[id]["exon"]=old_num+1
                if l[2]=='intron':
                    if contigcov<cov or contiglength<clen:
                        continue
                    else:
                        pass

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


                    old_num=stat_stored[id]["intron"]
                    stat_stored[id]["intron"]=old_num+1

                    if contigcov<float(cov) or contiglength<clen:
                        continue
                    else:
                        pass

                    'we store all combination genename/seqid/refid + frame'
                    #id = (genename,seqid,refid,refposmin,refposmax)
                    if dict_aln[genename].count(aln_info)==1:
                        dicframe=dict()
                        dicframe['frame']=str(frame)
                        stored.setdefault(id, []).append(dicframe)
                        dicinfo = dict()
                        minpos = l[3]
                        maxpos = l[4]
                        dicinfo['min']=str(minpos)
                        dicinfo['max']=str(maxpos)
                        stored.setdefault(id, []).append(dicinfo)
                    else:
                        pass
                else:
                    pass
        else:

            if contigcov<cov or contiglength<clen:
                continue
            else:
                pass

            ###____1_____
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

            if l[2]=="gene":

                if contigcov<cov or contiglength<clen:
                    continue
                else:
                    pass

                score = l[5]
                frame = l[6]
                id=(genename,seqid,refid,refposmin,refposmax,score,frame)

                if id in stat_stored.keys():
                    pass
                else:
                    stat_stored[id]={}
                    stat_stored[id]["exon"]=0
                    stat_stored[id]["intron"]=0

                'dictionary to identify the best score with combination seqID/refID'
                'we have to keep index in the same order of overlapp and create new dict'
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

                if typeseq=="all":
                    'we store all combination genename/seqid/refid + frame'
                    #id = (genename,seqid,refid,refposmin,refposmax)
                    if dict_aln[genename].count(aln_info)==1:
                        dicframe=dict()
                        dicframe['frame']=str(frame)
                        stored.setdefault(id, []).append(dicframe)
                        dicinfo = dict()
                        minpos = l[3]
                        maxpos = l[4]
                        dicinfo['min']=str(minpos)
                        dicinfo['max']=str(maxpos)
                        stored.setdefault(id, []).append(dicinfo)
                    else:
                        pass

                    if l[2]=='exon':
                        old_num=stat_stored[id]["exon"]
                        stat_stored[id]["exon"]=old_num+1
                    if l[2]=="intron":
                        old_num=stat_stored[id]["intron"]
                        stat_stored[id]["intron"]=old_num+1

            else:
                if typeseq=="exon":
                    if contigcov<cov or contiglength<clen:
                        continue
                    else:
                        pass

                    if dict_aln[genename].count(aln_info)==1:
                        dicframe=dict()
                        dicframe['frame']=str(frame)
                        stored.setdefault(id, []).append(dicframe)
                    else:
                        pass

                    if l[2]=="exon":
                        if contigcov<float(cov) or contiglength<clen:
                            continue
                        else:
                            pass

                        old_num=stat_stored[id]["exon"]
                        stat_stored[id]["exon"]=old_num+1

                        'we keep all cds position for each combination'
                        if dict_aln[genename].count(aln_info)==1:
                            dicinfo = dict()
                            minpos = l[3]
                            maxpos = l[4]
                            dicinfo['min']=str(minpos)
                            dicinfo['max']=str(maxpos)
                            stored.setdefault(id, []).append(dicinfo)
                        else:
                            pass

                    elif l[2]=="intron":
                        old_num=stat_stored[id]["intron"]
                        stat_stored[id]["intron"]=old_num+1

                elif typeseq=="intron":
                    if contigcov<cov or contiglength<clen:
                        continue
                    else:
                        pass

                    if dict_aln[genename].count(aln_info)==1:
                        dicframe=dict()
                        dicframe['frame']=str(frame)
                        stored.setdefault(id, []).append(dicframe)
                    else:
                        pass

                    if l[2]=="intron":
                        if contigcov<float(cov) or contiglength<clen:
                            continue
                        else:
                            pass

                        old_num=stat_stored[id]["intron"]
                        stat_stored[id]["intron"]=old_num+1

                        'we keep all cds position for each combination'
                        if dict_aln[genename].count(aln_info)==1:
                            dicinfo = dict()
                            minpos = l[3]
                            maxpos = l[4]
                            dicinfo['min']=str(minpos)
                            dicinfo['max']=str(maxpos)
                            stored.setdefault(id, []).append(dicinfo)
                        else:
                            pass

                    elif l[2]=="exon":
                        old_num=stat_stored[id]["exon"]
                        stat_stored[id]["exon"]=old_num+1


besthits_filtered=dict()
for g in list(besthits.keys()):
    besthits_filtered[g]=[]
    bhit_partcount=0
    bhitfilt_parcount=0
    for bhits in besthits[g]:
        if 'ITS' in g:
            besthits_filtered[g].append(bhits)
        else:
            bhit_partcount=bhit_partcount+1
            frame_bhits=[i["frame"] for i in stored[bhits] if "frame" in list(i.keys())][0]
            for hit in stored.keys():
                if hit[0]==bhits[0] and hit[1]==bhits[1] and hit[2]==bhits[2]:
                    if hit!=bhits:
                        #cond overlapp ratio 80%
                        minhit=hit[3]
                        maxhit=hit[4]
                        shit=int(hit[5])
                        sbhit=int(bhits[5])
                        rscore=float(shit)/float(sbhit)
                        rhit=set(range(minhit,maxhit))
                        rbhit=set(range(bhits[3],bhits[4]))
                        intersection=rhit.intersection(rbhit)
                        ratio=float(len(intersection))/float(len(rbhit))
                        frame_hit=[i["frame"] for i in stored[hit] if "frame" in list(i.keys())][0]
                        if ratio>0.8 and rscore>0.8:
                            if frame_bhits!=frame_hit:
                                if stat_stored[hit]['intron']<stat_stored[bhits]['intron']:
                                    besthits_filtered[g].append(hit)
                                    bhitfilt_parcount=bhitfilt_parcount+1
                                else:
                                    pass
                            else:
                                pass
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
            if bhitfilt_parcount<bhit_partcount:
                if bhits in besthits_filtered[g]:
                    pass
                else:
                    besthits_filtered[g].append(bhits)
                    bhitfilt_parcount=bhitfilt_parcount+1
            else:
                pass

## UNCOMMENT PART
# we add a new condition to choose no the best for the two rRNA but the contigs which maximise the capture of the all complex
besthits_filtered2=dict()
if "rrn5.8S" in besthits_filtered.keys():
    cont_rrn5=besthits_filtered["rrn5.8S"][0][1]
if "rrnITS1" in besthits_filtered.keys():
    cont_rnITS1=besthits_filtered["rrnITS1"][0][1]
if "rrnITS2" in besthits_filtered.keys():
    cont_rnITS2=besthits_filtered["rrnITS2"][0][1]

for g in list(besthits_filtered.keys()):
    if g=="rrn18S" or g=="rrn16S":
        contid=besthits_filtered[g][0][1]
        if "rrn5.8S" in besthits_filtered.keys() and "rrnITS1" in besthits_filtered.keys():
            if contid==cont_rrn5 and contid==cont_rnITS1:
                besthits_filtered2[g]=besthits_filtered[g]
            else:
                if cont_rrn5==cont_rnITS1:
                    hit_found=list()
                    for hit in stored.keys():
                        if hit[0]==g and hit[1]==cont_rrn5:
                            hit_found.append(hit)
                        else:
                            pass
                    if len(hit_found)>0:
                        score_max=0
                        hit_score=[]
                        for sub in hit_found:
                            sub_score=int(sub[5])
                            if sub_score>=score_max:
                                score_max=sub_score
                                hit_score=[sub]
                        besthits_filtered2[g]=hit_score
                    else:
                        pass
                else:
                    if contid==cont_rnITS1:
                        besthits_filtered2[g]=besthits_filtered[g]
                    else:
                        pass
        else:
            pass
    elif g=="rrn23S" or g=="rrn26S" or g=="rrn28S":
        contid=besthits_filtered[g][0][1]
        if "rrn5.8S" in besthits_filtered.keys() and "rrnITS2" in besthits_filtered.keys():
            if contid==cont_rrn5 and contid==cont_rnITS2:
                besthits_filtered2[g]=besthits_filtered[g]
            else:
                if cont_rrn5==cont_rnITS2:
                    hit_found=list()
                    for hit in stored.keys():
                        if hit[0]==g and hit[1]==cont_rrn5:
                            hit_found.append(hit)
                        else:
                            pass
                    if len(hit_found)>0:
                        score_max=0
                        hit_score=[]
                        for sub in hit_found:
                            sub_score=int(sub[5])
                            if sub_score>=score_max:
                                score_max=sub_score
                                hit_score=[sub]
                        besthits_filtered2[g]=hit_score
                    else:
                        pass
                else:
                    if contid==cont_rnITS2:
                        besthits_filtered2[g]=besthits_filtered[g]
                    else:
                        pass
        else:
            pass
    else:
        besthits_filtered2[g]=besthits_filtered[g]

## COMMENT PART
# # we add a new condition to choose no the best for the two rRNA but the contigs which maximise the capture of the all complex
# besthits_filtered2=dict()
# if "rrn5.8S" in besthits_filtered.keys():
#     cont_rrn5=besthits_filtered["rrn5.8S"][0][1]
# if "rrnITS1" in besthits_filtered.keys():
#     cont_rnITS1=besthits_filtered["rrnITS1"][0][1]
# if "rrnITS2" in besthits_filtered.keys():
#     cont_rnITS2=besthits_filtered["rrnITS2"][0][1]
#
# for g in list(besthits_filtered.keys()):
#     if g=="rrn18S" or g=="rrn16S":
#         contid=besthits_filtered[g][0][1]
#         if "rrn5.8S" and "rrnITS1" in besthits_filtered.keys():
#             if contid==cont_rrn5 and contid==cont_rnITS1:
#                 besthits_filtered2[g]=besthits_filtered[g]
#             else:
#                 if cont_rrn5==cont_rnITS1:
#                     hit_found=list()
#                     for hit in stored.keys():
#                         if hit[0]==g and hit[1]==cont_rrn5:
#                             hit_found.append(hit)
#                         else:
#                             pass
#                     if len(hit_found)>0:
#                         score_max=0
#                         hit_score=[]
#                         for sub in hit_found:
#                             sub_score=int(sub[5])
#                             if sub_score>=score_max:
#                                 score_max=sub_score
#                                 hit_score=[sub]
#                         besthits_filtered2[g]=hit_score
#                     else:
#                         besthits_filtered2[g]=besthits_filtered[g]
#                 else:
#                     besthits_filtered2[g]=besthits_filtered[g]
#         else:
#             besthits_filtered2[g]=besthits_filtered[g]
#     elif g=="rrn23S" or g=="rrn26S" or g=="rrn28S":
#         contid=besthits_filtered[g][0][1]
#         if "rrn5.8S" and "rrnITS2" in besthits_filtered.keys():
#             if contid==cont_rrn5 and contid==cont_rnITS2:
#                 besthits_filtered2[g]=besthits_filtered[g]
#             else:
#                 if cont_rrn5==cont_rnITS2:
#                     hit_found=list()
#                     for hit in stored.keys():
#                         if hit[0]==g and hit[1]==cont_rrn5:
#                             hit_found.append(hit)
#                         else:
#                             pass
#                     if len(hit_found)>0:
#                         score_max=0
#                         hit_score=[]
#                         for sub in hit_found:
#                             sub_score=int(sub[5])
#                             if sub_score>=score_max:
#                                 score_max=sub_score
#                                 hit_score=[sub]
#                         besthits_filtered2[g]=hit_score
#                     else:
#                         besthits_filtered2[g]=besthits_filtered[g]
#                 else:
#                     besthits_filtered2[g]=besthits_filtered[g]
#         else:
#             besthits_filtered2[g]=besthits_filtered[g]
#     else:
#         besthits_filtered2[g]=besthits_filtered[g]


rrna_contigs=list()
for g in list(besthits_filtered2.keys()):
    if "ITS" in g:
        pass
    else:
        for hits in besthits_filtered2[g]:
            seqid=hits[1]
            rrna_contigs.append(seqid)

if len(list(besthits_filtered2.keys()))>0:
    print ("Extraction of sequences for %s" % nameofsample)
    inputs=range(len(list(besthits_filtered2.keys())))
    x=Parallel(n_jobs=num_cores)(delayed(GeneExtraction)(i) for i in inputs)

    l = [i for i in x if i is not None]
    outcontp=os.path.dirname(os.path.abspath(file))
    if len(l)==0:
        print("WARN: "+model+" - no genes following capture restrictions were extracted for %s" % nameofsample)
    else:
        pass
    if 'nucrdna' in model:
        if os.path.join(outcontp, str(nameofsample+".cont_nucrdna.log")):
            with open(os.path.join(outcontp, str(nameofsample+".cont_nucrdna.log")), 'a+') as outcontt:
                old_cont = []
                end_file=outcontt.tell()
                outcontt.seek(0)
                for line in outcontt:
                    old_cont.append(line.rstrip())
                if len(l)>0:
                    for line in set(l):
                        if not line in old_cont:
                            outcontt.seek(end_file)
                            outcontt.write(str(line)+'\n')
                        else:
                            pass
                else:
                    outcontt.seek(end_file)
                    outcontt.write(str("None")+'\n')
        else :
            with open(os.path.join(outcontp, str(nameofsample+".cont_nucrdna.log")), 'w') as outcontt:
                if len(l)>0:
                    for line in set(l):
                        outcontt.write(str(line)+'\n')
                else:
                    outcontt.write(str("None")+'\n')
    else:
        pass
else:
    print("WARN: "+model+" - no genes following capture restrictions were extracted for %s" % nameofsample)
    outcontp=os.path.dirname(os.path.abspath(file))
    if os.path.join(outcontp, str(nameofsample+".cont_nucrdna.log")):
        pass
    else:
        with open(os.path.join(outcontp, str(nameofsample+".cont_nucrdna.log")), 'w') as outcontt:
            outcontt.write("None")

print("--- %s seconds ---" % (time.time() - start_time))
