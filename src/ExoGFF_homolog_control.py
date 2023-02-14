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
import numpy as N

from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna
start_time = time.time()

parser = argparse.ArgumentParser(description='Capture regions of query reference genes from exonerate alignment. Script was writen by C. Pouchon (2019).')
parser.add_argument("-i","--infile", help="Fasta input contigs file",
                    type=str)
parser.add_argument("-m","--model", help="molecular type compartment",
                    type=str,choices=["chloroplast", "mitochondrion_CDS","mitochondrion_rRNA","mitochondrion_tRNA","nucleus_aa","nucleus_nt","chloroplast_CDS","chloroplast_rRNA","chloroplast_tRNA","busco","uce"])
parser.add_argument("-o","--outdir", help="Out directory path",
                    type=str)
parser.add_argument("-n","--namesample", help="sample name used in sequence header for output",
                    type=str)
parser.add_argument("-g","--gfftab", help="input gff table from exonerate alignment",
                    type=str)
parser.add_argument("--genetic_code", help="NCBI genetic code number used for DNA translation",
                    type=int)
parser.add_argument("-clt","--control", help="input gff table from exonerate alignment to seeds from other organelle (mitochondrial or chloroplastic)",
                    type=str)
parser.add_argument("-t","--typeofseq", help="type of sequence to extract",
                    type=str,choices=["all", "exon","intron"])
parser.add_argument("--cov_cutoff", help="coverage cut-off option for organelles according to a standard deviation of mean coverage",
                    type=str,choices=["on", "off"])
parser.add_argument("-l","--minlength", help="minimal length of alignment into reference to consider output",
                    type=int)
parser.add_argument("-rp","--refpct", help="minimal percent of reference covered allowed",
                    type=float)
parser.add_argument("--orfcov", help="minimal fraction of captured sequenced covered by ORF",
                    type=float)
parser.add_argument("-c","--mincov", help="minimal coverage of contigs allowed for genomic scan",
                    type=int)
parser.add_argument("-cl","--mincontlen", help="minimal length of contifs allowed for genomic scan",
                    type=int)
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

def mad(a, axis=None):
    """
    Compute *Median Absolute Deviation* of an array along given axis.
    """

    # Median along given axis, but *keeping* the reduced axis so that
    # result can still broadcast against a.
    med = N.median(a, axis=axis, keepdims=True)
    mad = N.median(N.absolute(a - med), axis=axis)  # MAD along given axis

    return mad

def weighted_std(values, weights):
    """
    Return the weighted average and standard deviation.
    values, weights -- Numpy ndarrays with the same shape.
    """
    average = N.average(values, weights=weights)
    # Fast and numerically precise:
    variance = N.average((values-average)**2, weights=weights)
    return N.sqrt(variance)


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

def pad_seq(sequence):
    """ Pad sequence to multiple of 3 with N """

    remainder = len(sequence) % 3

    return sequence if remainder == 0 else sequence + Seq('N' * (3 - remainder))

def GeneExtraction(genenumber):
    contpart=[]
    num_exon=[]
    num_intron=[]
    geneid=list(besthits_filtered5.keys())[genenumber]
    extract_contigs=list()

    if len(besthits_filtered5[geneid])==1:
        refcoverage_dict={}
        refcoverage_pass=[]

        seqid = besthits_filtered5[geneid][0][1]
        concat=[]

        for hits in besthits_filtered5[geneid]:
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

            if len(stored[besthits_filtered5[geneid][0]])==0:
                pass
            else:
                for parts in stored[besthits_filtered5[geneid][0]]:
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

                num_exon.append(str(stat_stored[besthits_filtered5[geneid][0]]["exon"]))
                num_intron.append(str(stat_stored[besthits_filtered5[geneid][0]]["intron"]))

    else:
        'condition pour voir quelle position il faut extraire en premier'
        lmin=list()

        refcoverage_dict={}
        refcoverage_pass=[]

        for hits in besthits_filtered5[geneid]:
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
            seqid = besthits_filtered5[geneid][idx][1]
            if seqid in refcoverage_pass:
                contpart.append(str(len(set(counthits[geneid][idx]))))
                extract_contigs.append(seqid)
                #fr = stored[besthits[geneid][idx]][0]['frame']
                dna=str(seqs[seqid][0])
                list_combo=list()
                pos_to_extract=list()
                if len(stored[besthits_filtered5[geneid][idx]])==0:
                    pass
                else:
                    for parts in stored[besthits_filtered5[geneid][idx]]:
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

                num_exon.append(str(stat_stored[besthits_filtered5[geneid][idx]]["exon"]))
                num_intron.append(str(stat_stored[besthits_filtered5[geneid][idx]]["intron"]))
            else:
                continue

    newseq="".join(concat)

    newlength = len(newseq)

    if typeseq=="exon":
        genepart=typeseq
    elif typeseq=="intron":
        genepart=typeseq
    elif typeseq=="all":
        genepart="exon+intron"

    cond_min_length=0
    if newlength>=min_length:
        cond_min_length=cond_min_length+1

    cond_toprint_var=0
    cond_orf=0

    if typeseq=="exon" and (model=="chloroplast_CDS" or model=="mitochondrion_CDS" or model=="busco" or model=="nucleus_aa"):
        #tseq=Seq(newseq).translate(table=genet_code)
        tseq=pad_seq(Seq(newseq)).translate(table=genet_code)
        if "*" in tseq:
            newseq_rev=str(Seq(newseq).reverse_complement())
            seqframes=[pad_seq(Seq(newseq)).translate(table=genet_code),pad_seq(Seq(newseq[1:len(newseq)])).translate(table=genet_code),pad_seq(Seq(newseq[2:len(newseq)])).translate(table=genet_code),pad_seq(Seq(newseq_rev)).translate(table=genet_code),pad_seq(Seq(newseq_rev[1:len(newseq_rev)])).translate(table=genet_code),pad_seq(Seq(newseq_rev[2:len(newseq_rev)])).translate(table=genet_code)]
            frame_select=0
            len_orf=0
            orf_select=""
            for p in range(len(seqframes)):
                orf_parts=[]
                start=0
                for i in range(len(str(seqframes[p]))):
                    if seqframes[p][i]!="*":
                        end=i
                    else:
                        end=i
                        orf_parts.append(str(start)+"-"+str(end))
                        start=i+1
                orf_parts.append(str(start)+"-"+str(end+1))
                for i in range(len(orf_parts)):
                    part=orf_parts[i]
                    s=int(part.split("-")[0])
                    e=int(part.split("-")[1])
                    if len(seqframes[p][s:e])>len_orf and len(seqframes[p][s:e])>=(float(orfcond)*len(seqframes[p])):
                        len_orf=len(seqframes[p][s:e])
                        orf_select=part
                        frame_select=p
            #check that orf cov refpct
            if len(orf_select)>0:
                newseq2=newseq
                cond_orf=cond_orf+1
            else:
                # taille ORF non conforme donc gene-like
                cond_toprint_var=cond_toprint_var+1
                newseq2=newseq
        else:
            newseq2=newseq
    elif typeseq=="all" and (model=="chloroplast_CDS" or model=="mitochondrion_CDS" or model=="busco" or model=="nucleus_aa"):
        if cond_min_length!=0:
            if int(max(num_exon))>1:
                # condition si intron dans la séquence, on recupere les sequences d'exon, sans cond_orf
                lmin=list()
                refcoverage_dict={}
                refcoverage_pass=[]
                for hits in besthits_filtered5[geneid]:
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
                exon_concat=[]
                for idx in orderindx:
                    seqid = besthits_filtered5[geneid][idx][1]
                    if seqid in refcoverage_pass:
                        contpart.append(str(len(set(counthits[geneid][idx]))))
                        extract_contigs.append(seqid)
                        #fr = stored[besthits[geneid][idx]][0]['frame']
                        dna=str(seqs[seqid][0])
                        list_combo=list()
                        pos_to_extract=list()
                        if len(stored[besthits_filtered5[geneid][idx]])==0:
                            pass
                        else:
                            for parts in stored[besthits_filtered5[geneid][idx]]:
                                if 'frame' in parts.keys():
                                    fr=parts['frame']
                                else:
                                    pass
                            for expart in stat_stored_position[besthits_filtered5[geneid][idx]]["exon"]:
                                start=int(expart.split("-")[0])
                                end=int(expart.split("-")[1])
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
                                exon_concat.append(extract)
                exon_seq="".join(exon_concat)
                newseq_rev=str(Seq(newseq).reverse_complement())
                tseq=pad_seq(Seq(exon_seq)).translate(table=genet_code)
                if "*" in tseq:
                    exon_seq_rev=str(Seq(exon_seq).reverse_complement())
                    seqframes=[pad_seq(Seq(exon_seq)).translate(table=genet_code),pad_seq(Seq(exon_seq[1:len(exon_seq)])).translate(table=genet_code),pad_seq(Seq(exon_seq[2:len(exon_seq)])).translate(table=genet_code),pad_seq(Seq(exon_seq_rev)).translate(table=genet_code),pad_seq(Seq(exon_seq_rev[1:len(exon_seq_rev)])).translate(table=genet_code),pad_seq(Seq(exon_seq_rev[2:len(exon_seq_rev)])).translate(table=genet_code)]
                    frame_select=0
                    len_orf=0
                    orf_select=""
                    for p in range(len(seqframes)):
                        orf_parts=[]
                        start=0
                        for i in range(len(str(seqframes[p]))):
                            if seqframes[p][i]!="*":
                                end=i
                            else:
                                end=i
                                orf_parts.append(str(start)+"-"+str(end))
                                start=i+1
                        orf_parts.append(str(start)+"-"+str(end+1))
                        for i in range(len(orf_parts)):
                            part=orf_parts[i]
                            s=int(part.split("-")[0])
                            e=int(part.split("-")[1])
                            if len(seqframes[p][s:e])>len_orf and len(seqframes[p][s:e])>=(float(orfcond)*len(seqframes[p])):
                                len_orf=len(seqframes[p][s:e])
                                orf_select=part
                                frame_select=p
                    #check that orf cov refpct
                    if len(orf_select)>0:
                        newseq2=newseq
                        cond_orf=cond_orf+1
                    else:
                        cond_toprint_var=cond_toprint_var+1
                        newseq2=newseq
                else:
                    newseq2=newseq
            else:
                tseq=pad_seq(Seq(newseq)).translate(table=genet_code)
                newseq_rev=str(Seq(newseq).reverse_complement())
                if "*" in tseq:
                    seqframes=[pad_seq(Seq(newseq)).translate(table=genet_code),pad_seq(Seq(newseq[1:len(newseq)])).translate(table=genet_code),pad_seq(Seq(newseq[2:len(newseq)])).translate(table=genet_code),pad_seq(Seq(newseq_rev)).translate(table=genet_code),pad_seq(Seq(newseq_rev[1:len(newseq_rev)])).translate(table=genet_code),pad_seq(Seq(newseq_rev[2:len(newseq_rev)])).translate(table=genet_code)]
                    frame_select=0
                    len_orf=0
                    orf_select=""
                    for p in range(len(seqframes)):
                        orf_parts=[]
                        start=0
                        for i in range(len(str(seqframes[p]))):
                            if seqframes[p][i]!="*":
                                end=i
                            else:
                                end=i
                                orf_parts.append(str(start)+"-"+str(end))
                                start=i+1
                        orf_parts.append(str(start)+"-"+str(end+1))
                        for i in range(len(orf_parts)):
                            part=orf_parts[i]
                            s=int(part.split("-")[0])
                            e=int(part.split("-")[1])
                            if len(seqframes[p][s:e])>len_orf and len(seqframes[p][s:e])>=(float(orfcond)*len(seqframes[p])):
                                len_orf=len(seqframes[p][s:e])
                                orf_select=part
                                frame_select=p
                    #check that orf cov refpct
                    if len(orf_select)>0:
                        newseq2=newseq
                        cond_orf=cond_orf+1
                    else:
                        cond_toprint_var=cond_toprint_var+1
                        newseq2=newseq
                else:
                    newseq2=newseq
        else:
            newseq2=newseq
    else:
        newseq2=newseq

    fname = geneid+".fa"

    exon_l_recovered=0
    for p in refcoverage_pass:
        exon_l_recovered=exon_l_recovered+int(refcoverage_dict[p])
    if cond_orf<1:
        gpct=round(float(exon_l_recovered)/float(reflength[geneid]),2)
    else:
        if model=="nucleus_nt":
            gpct=round(float(len(newseq2)/3)/float(reflength[geneid]/3),2)
        else:
            gpct=round(float(len(newseq2)/3)/float(reflength[geneid]),2)
    newlength=len(newseq2)

    if cond_min_length!=0:
        if cond_toprint_var<1:
            header = ">"+str(nameofsample)+"; "+"gene="+str(geneid)+"; "+"info="+str(genepart)+"; "+"type="+str(model)+"; "+"length="+str(newlength)+"; "+"match_contigs="+str("-".join(contpart))+"; "+"ref_percent="+str(gpct)+"; "+"n_exons="+str("-".join(num_exon))+"; "+"n_introns="+str("-".join(num_intron))
            if "RNA" in str(model):
                if len(contpart)>1 or len(num_exon)>1 or len(num_intron)>1:
                    pass
                else:
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
                                file.write(str(newseq2)+'\n')
                            else:
                                pass
                    else:
                        with open(os.path.join(outpath+"/"+model, fname), 'w') as out:
                            out.write(header+'\n')
                            out.write(str(newseq2)+'\n')

                    for elem in extract_contigs:
                        return elem
            else:
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
                            file.write(str(newseq2)+'\n')
                        else:
                            pass
                else:
                    with open(os.path.join(outpath+"/"+model, fname), 'w') as out:
                        out.write(header+'\n')
                        out.write(str(newseq2)+'\n')

                for elem in extract_contigs:
                    return elem
        else:
            print("WARN: min ORF length condition not filled for %s in %s. The sequence is given in %s" % (str(nameofsample),str(geneid),str(geneid)+"-like.fa"))
            fname = geneid+"-like.fa"
            header = ">"+str(nameofsample)+"; "+"gene="+str(geneid)+"-like"+"; "+"info="+str(genepart)+"; "+"type="+str(model)+"; "+"length="+str(newlength)+"; "+"match_contigs="+str("-".join(contpart))+"; "+"ref_percent="+str(gpct)+"; "+"n_exons="+str("-".join(num_exon))+"; "+"n_introns="+str("-".join(num_intron))
            if "RNA" in str(model):
                if len(contpart)>1 or len(num_exon)>1 or len(num_intron)>1:
                    pass
                else:
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
                                file.write(str(newseq2)+'\n')
                            else:
                                pass
                    else:
                        with open(os.path.join(outpath+"/"+model, fname), 'w') as out:
                            out.write(header+'\n')
                            out.write(str(newseq2)+'\n')

                    for elem in extract_contigs:
                        return elem
            else:
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
                            file.write(str(newseq2)+'\n')
                        else:
                            pass
                else:
                    with open(os.path.join(outpath+"/"+model, fname), 'w') as out:
                        out.write(header+'\n')
                        out.write(str(newseq2)+'\n')

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
control=args.control
option_cutoff=args.cov_cutoff
orfcond=args.orfcov
genet_code=args.genetic_code

stored={}
stat_stored={}
dicscore={}
besthits={}
groups={}
seqs={}
counthits={}
dict_aln={}
contigs_homolog=dict()
stat_stored_position={}

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
        seq=list()
        if genename not in refnames:
            refnames.append(genename)
        else:
            pass
        if "BUSCO" in l:
            genename = l.replace(">","").split("_")[0]+"_BUSCO"
        else:
            genename = l.replace(">","").split("_")[0]
    else:
        l=line.rstrip()
        seq.append(l)
        reflength[genename]=len("".join(seq))

'we store all sequences'
cur_genome = SeqIO.parse(input_file, "fasta")
for record in cur_genome:
    seqID=record.id
    sequence=record.seq
    seqs.setdefault(seqID, []).append(sequence)


'find each gene of ref and keep seqID,score,length,ref'
for line in tab:
    l = line.rstrip().split("\t")
    if l[0] in refnames and "Target" in l[len(l)-1]:
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

    if typeseq=="all":

        if l[2]=="gene":

            if contigcov<cov or contiglength<clen:
                continue
            else:
                pass

            score = l[5]
            frame = l[6]
            cont_min=l[3]
            cont_max=l[4]
            id=(genename,seqid,refid,refposmin,refposmax,score)

            if seqid in contigs_homolog.keys():
                toappend=dict()
                toappend["min"]=int(cont_min)
                toappend["max"]=int(cont_max)
                toappend["gene"]=genename
                toappend["score"]=score
                contigs_homolog[seqid].append(toappend)
            else:
                contigs_homolog[seqid]=[]
                toappend=dict()
                toappend["min"]=int(cont_min)
                toappend["max"]=int(cont_max)
                toappend["gene"]=genename
                toappend["score"]=score
                contigs_homolog[seqid].append(toappend)

            if id in stat_stored.keys():
                pass
            else:
                stat_stored[id]={}
                stat_stored[id]["exon"]=0
                stat_stored[id]["intron"]=0

            if id in stat_stored_position.keys():
                pass
            else:
                stat_stored_position[id]={}
                stat_stored_position[id]["exon"]=[]
                stat_stored_position[id]["intron"]=[]

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

        elif l[2]=='exon':
            if contigcov<cov or contiglength<clen:
                continue
            else:
                pass

            old_num=stat_stored[id]["exon"]
            stat_stored[id]["exon"]=old_num+1
            posexonmin=l[3]
            poseexonmax=l[4]
            stat_stored_position[id]["exon"].append(str(posexonmin)+"-"+str(poseexonmax))

        elif l[2]=="intron":
            if contigcov<cov or contiglength<clen:
                continue
            else:
                pass

            old_num=stat_stored[id]["intron"]
            stat_stored[id]["intron"]=old_num+1

    else:
        if l[2]=="gene":

            if contigcov<cov or contiglength<clen:
                continue
            else:
                pass

            score = l[5]
            frame = l[6]
            cont_min=l[3]
            cont_max=l[4]
            id=(genename,seqid,refid,refposmin,refposmax,score)

            if seqid in contigs_homolog.keys():
                toappend=dict()
                toappend["min"]=int(cont_min)
                toappend["max"]=int(cont_max)
                toappend["gene"]=genename
                toappend["score"]=score
                contigs_homolog[seqid].append(toappend)
            else:
                contigs_homolog[seqid]=[]
                toappend=dict()
                toappend["min"]=int(cont_min)
                toappend["max"]=int(cont_max)
                toappend["gene"]=genename
                toappend["score"]=score
                contigs_homolog[seqid].append(toappend)

            if id in stat_stored.keys():
                pass
            else:
                stat_stored[id]={}
                stat_stored[id]["exon"]=0
                stat_stored[id]["intron"]=0

            if id in stat_stored_position.keys():
                pass
            else:
                stat_stored_position[id]={}
                stat_stored_position[id]["exon"]=[]
                stat_stored_position[id]["intron"]=[]

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


            'we store all combination genename/seqid/refid + frame'
            #id = (genename,seqid,refid,refposmin,refposmax)
            if dict_aln[genename].count(aln_info)==1:
                dicframe=dict()
                dicframe['frame']=str(frame)
                stored.setdefault(id, []).append(dicframe)
            else:
                pass

        elif typeseq=="exon":
            if l[2]==typeseq:
                if contigcov<cov or contiglength<clen:
                    continue
                else:
                    pass

                old_num=stat_stored[id]["exon"]
                stat_stored[id]["exon"]=old_num+1
                posexonmin=l[3]
                poseexonmax=l[4]
                stat_stored_position[id]["exon"].append(str(posexonmin)+"-"+str(poseexonmax))

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
                if contigcov<cov or contiglength<clen:
                    continue
                else:
                    pass

                old_num=stat_stored[id]["intron"]
                stat_stored[id]["intron"]=old_num+1

        elif typeseq=="intron":
            if l[2]==typeseq:
                if contigcov<cov or contiglength<clen:
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
                if contigcov<cov or contiglength<clen:
                    continue
                else:
                    pass
                old_num=stat_stored[id]["exon"]
                stat_stored[id]["exon"]=old_num+1


besthits_filtered=dict()
for g in list(besthits.keys()):
    besthits_filtered[g]=[]
    bhit_partcount=0
    bhitfilt_parcount=0
    for bhits in besthits[g]:
        tmp_score=0
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
                                if tmp_score==0:
                                    tmp_score=rscore
                                    tmp_hit=hit
                                else:
                                    if rscore>tmp_score:
                                        tmp_score=rscore
                                        tmp_hit=hit
                                #besthits_filtered[g].append(hit)
                                #bhitfilt_parcount=bhitfilt_parcount+1
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
        if tmp_score>0:
            besthits_filtered[g].append(tmp_hit)
            bhitfilt_parcount=bhitfilt_parcount+1
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

if "RNA" in model:
    besthits_filtered2=dict()
    for g in list(besthits_filtered.keys()):
        subpart=[]
        for bhits in besthits_filtered[g]:
            if stat_stored[bhits]['intron']<3:
                subpart.append(bhits)
            else:
                pass
        if len(subpart)>0:
            besthits_filtered2[g]=subpart
        else:
            pass
else:
    besthits_filtered2=besthits_filtered


cstored={}
cstat_stored={}
cdicscore={}
cbesthits={}
cgroups={}
cseqs={}
ccounthits={}
cdict_aln={}

'we store the gff tab'
with open(control) as f:
    controltab = f.readlines()

'find each gene of ref and keep seqID,score,length,ref'
for line in controltab:
    l = line.rstrip().split("\t")
    if "Target" in l[len(l)-1]:
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
        cdict_aln.setdefault(genename, []).append(aln_info)

        contiglength=int(seqid.split("_")[3])
        contigcov=float(seqid.split("_")[5])

        if contigcov<cov or contiglength<clen:
            continue
        else:
            pass

        'we make groups for gene with borne values'
        if genename not in cgroups.keys():
            dic=[]
            dic.append(rmin)
            dic.append(rmax)
            cgroups.setdefault(genename, []).append(dic)
        else:
            'if existed, we scan each subgrp to check if the hits overlapp others hits on reference sequence'
            loverlap=list()
            lnoverlap=list()

            for i in range(len(cgroups[genename])):
                m = min(cgroups[genename][i])
                M = max(cgroups[genename][i])
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
                novergp.append(cgroups[genename][i])
            novergp.append(overgp)
            cgroups[genename]=novergp


        contname=[]
        contname.append(seqid)
        if genename not in ccounthits.keys():
             ccounthits.setdefault(genename, []).append(contname)
        else:
            newcont=[]

            if len(lnoverlap)>0:
                for i in lnoverlap:
                     newcont.append(ccounthits[genename][i])

            overcont=[]
            'il faut modifier pour avoir que les element est pas une liste dans une liste'
            if len(loverlap)>0:
                for i in loverlap:
                    overcont.extend(ccounthits[genename][i])
                if seqid not in overcont:
                    overcont.extend(contname)
                newcont.append(overcont)
            else:
                newcont.append(contname)

            ccounthits[genename]=newcont

    if typeseq=="all":

        if l[2]=="gene":

            if contigcov<cov or contiglength<clen:
                continue
            else:
                pass

            score = l[5]
            frame = l[6]
            id=(genename,seqid,refid,refposmin,refposmax,score)

            if id in cstat_stored.keys():
                pass
            else:
                cstat_stored[id]={}
                cstat_stored[id]["exon"]=0
                cstat_stored[id]["intron"]=0

            if genename not in cdicscore.keys():
                cdicscore.setdefault(genename, []).append(int(score))
                #besthits.setdefault(genename, []).append(dicinfo)
                cbesthits.setdefault(genename, []).append(id)
                #counthits.setdefault(genename, []).append([seqid])
            else:
                'mettre condition pour regarder dans loverlapp index'

                'condition si nouveau group'

                newscore=[]
                newhits=[]
                overscore=[]

                for i in lnoverlap:
                    'ceux pour qui pas overlapp on les gardes dans le meme ordre (new/modif in first)'
                    newscore.append(cdicscore[genename][i])
                    newhits.append(cbesthits[genename][i])

                if len(loverlap)>0:
                    for i in loverlap:
                        overscore.append(cdicscore[genename][i])
                    'extract best score on new overlapping group'
                    highscore_index = loverlap[overscore.index(max(overscore))]
                    highscore = cdicscore[genename][highscore_index]
                    'compare hit score with the bestscore for this overlapping group'
                    if int(score) > highscore:
                        highscore = int(score)
                        #highhits = dicinfo
                        highhits = id
                    else:
                        highhits = cbesthits[genename][highscore_index]

                    'maintenant faut faire le nouveau gp'
                    newscore.append(highscore)
                    newhits.append(highhits)
                else:
                    'pas overlapp donc new'
                    newscore.append(int(score))
                    #newhits.append(dicinfo)
                    newhits.append(id)


                cdicscore[genename]=newscore
                cbesthits[genename]=newhits


            'we store all combination genename/seqid/refid + frame'
            #id = (genename,seqid,refid,refposmin,refposmax)
            if cdict_aln[genename].count(aln_info)==1:
                dicframe=dict()
                dicframe['frame']=str(frame)
                cstored.setdefault(id, []).append(dicframe)
                dicinfo = dict()
                minpos = l[3]
                maxpos = l[4]
                dicinfo['min']=str(minpos)
                dicinfo['max']=str(maxpos)
                cstored.setdefault(id, []).append(dicinfo)
            else:
                pass

        elif l[2]=='exon':
            if contigcov<cov or contiglength<clen:
                continue
            else:
                pass
            old_num=cstat_stored[id]["exon"]
            cstat_stored[id]["exon"]=old_num+1

        elif l[2]=="intron":
            if contigcov<cov or contiglength<clen:
                continue
            else:
                pass
            old_num=cstat_stored[id]["intron"]
            cstat_stored[id]["intron"]=old_num+1

    else:
        if l[2]=="gene":

            if contigcov<cov or contiglength<clen:
                continue
            else:
                pass

            score = l[5]
            frame = l[6]
            id=(genename,seqid,refid,refposmin,refposmax,score)

            if id in cstat_stored.keys():
                pass
            else:
                cstat_stored[id]={}
                cstat_stored[id]["exon"]=0
                cstat_stored[id]["intron"]=0

            'dictionary to identify the best score with combination seqID/refID'
            'we have to keep index in the same order of overlapp and create new dict'
            if genename not in cdicscore.keys():
                cdicscore.setdefault(genename, []).append(int(score))
                #besthits.setdefault(genename, []).append(dicinfo)
                cbesthits.setdefault(genename, []).append(id)
                #counthits.setdefault(genename, []).append([seqid])
            else:
                'mettre condition pour regarder dans loverlapp index'

                'condition si nouveau group'

                newscore=[]
                newhits=[]
                overscore=[]

                for i in lnoverlap:
                    'ceux pour qui pas overlapp on les gardes dans le meme ordre (new/modif in first)'
                    newscore.append(cdicscore[genename][i])
                    newhits.append(cbesthits[genename][i])

                if len(loverlap)>0:
                    for i in loverlap:
                        overscore.append(cdicscore[genename][i])
                    'extract best score on new overlapping group'
                    highscore_index = loverlap[overscore.index(max(overscore))]
                    highscore = cdicscore[genename][highscore_index]
                    'compare hit score with the bestscore for this overlapping group'
                    if int(score) > highscore:
                        highscore = int(score)
                        #highhits = dicinfo
                        highhits = id
                    else:
                        highhits = cbesthits[genename][highscore_index]

                    'maintenant faut faire le nouveau gp'
                    newscore.append(highscore)
                    newhits.append(highhits)
                else:
                    'pas overlapp donc new'
                    newscore.append(int(score))
                    #newhits.append(dicinfo)
                    newhits.append(id)


                cdicscore[genename]=newscore
                cbesthits[genename]=newhits


            'we store all combination genename/seqid/refid + frame'
            #id = (genename,seqid,refid,refposmin,refposmax)
            if cdict_aln[genename].count(aln_info)==1:
                dicframe=dict()
                dicframe['frame']=str(frame)
                cstored.setdefault(id, []).append(dicframe)
            else:
                pass

        elif typeseq=="exon":

            if l[2]==typeseq:
                if contigcov<cov or contiglength<clen:
                    continue
                else:
                    pass

                old_num=cstat_stored[id]["exon"]
                cstat_stored[id]["exon"]=old_num+1

                'we keep all cds position for each combination'
                if cdict_aln[genename].count(aln_info)==1:
                    dicinfo = dict()
                    minpos = l[3]
                    maxpos = l[4]
                    dicinfo['min']=str(minpos)
                    dicinfo['max']=str(maxpos)
                    cstored.setdefault(id, []).append(dicinfo)
                else:
                    pass

            elif l[2]=="intron":
                if contigcov<cov or contiglength<clen:
                    continue
                else:
                    pass

                old_num=cstat_stored[id]["intron"]
                cstat_stored[id]["intron"]=old_num+1

        elif typeseq=="intron":
            if l[2]==typeseq:
                if contigcov<cov or contiglength<clen:
                    continue
                else:
                    pass

                old_num=cstat_stored[id]["intron"]
                cstat_stored[id]["intron"]=old_num+1

                'we keep all cds position for each combination'
                if cdict_aln[genename].count(aln_info)==1:
                    dicinfo = dict()
                    minpos = l[3]
                    maxpos = l[4]
                    dicinfo['min']=str(minpos)
                    dicinfo['max']=str(maxpos)
                    cstored.setdefault(id, []).append(dicinfo)
                else:
                    pass

            elif l[2]=="exon":
                if contigcov<cov or contiglength<clen:
                    continue
                else:
                    pass
                old_num=cstat_stored[id]["exon"]
                cstat_stored[id]["exon"]=old_num+1



#il faut faire un nouveau dict control avec les positions des contigs
control_contigs=dict()
for g in list(cbesthits.keys()):
    for hit in cbesthits[g]:
        ccontig=hit[1]
        if ccontig not in control_contigs.keys():
            control_contigs[ccontig]=[]
        else:
            pass
        cscore=int(hit[5])
        for elem in cstored[hit]:
            if 'frame' in elem:
                pass
            else:
                dicttoappend=dict()
                dicttoappend["min"]=int(elem["min"])
                dicttoappend["max"]=int(elem["max"])
                dicttoappend["gene"]=g
                dicttoappend["score"]=cscore
                control_contigs[ccontig].append(dicttoappend)


# condition control for homologs/paralog in the dataset.
#We check that the position for the best contig selected for each Gene
#don't match with overlapp position in an other gene (e.g. with two copies of the gene)
#In such case, we select only the gene copy with the best score

besthits_filtered3=dict()
for g in list(besthits_filtered2.keys()):
    subpart=[]
    for bhits in besthits_filtered2[g]:
        seqid=bhits[1]
        bscore=int(bhits[5])
        control_pass=0
        if seqid in list(contigs_homolog.keys()):
            for elem in stored[bhits]:
                if "frame" in elem:
                    pass
                else:
                    minhit=int(elem["min"])
                    maxhit=int(elem["max"])
                    rhit=set(range(minhit,maxhit))
                    size=int(maxhit-minhit+1)
                    for c in contigs_homolog[seqid]:
                        cmin=c["min"]
                        cmax=c["max"]
                        rc=set(range(cmin,cmax))
                        intersection=rhit.intersection(rc)
                        cg=c["gene"]
                        cgsize=int(cmax-cmin+1)
                        if len(intersection)/size>0.2:
                            if g == cg:
                                pass
                            else:
                                if bscore>int(c["score"]):
                                    pass
                                else:
                                    if len(intersection)/cgsize>0.2:
                                        control_pass=control_pass+1
                        else:
                            pass
        else:
            pass
        if control_pass==0:
            subpart.append(bhits)
        else:
            pass

    if len(subpart)>0:
        besthits_filtered3[g]=subpart
    else:
        pass


#control with mito/chloro seeds for the other organelle to be sure that we recover the true copy from the right plast
besthits_filtered4=dict()
for g in list(besthits_filtered3.keys()):
    subpart=[]
    for bhits in besthits_filtered3[g]:
        seqid=bhits[1]
        bscore=int(bhits[5])
        control_pass=0
        if seqid in list(control_contigs.keys()):
            # il faut les positions stored du bhits et comparer pour chacune a chaque pos du control_contigs
            for elem in stored[bhits]:
                if "frame" in elem:
                    pass
                else:
                    minhit=int(elem["min"])
                    maxhit=int(elem["max"])
                    rhit=set(range(minhit,maxhit))
                    for c in control_contigs[seqid]:
                        cmin=c["min"]
                        cmax=c["max"]
                        rc=set(range(cmin,cmax))
                        intersection=rhit.intersection(rc)
                        if len(intersection)>0:
                            if bscore>=c["score"]:
                                pass
                            else:
                                control_pass=control_pass+1
                        else:
                            pass
        else:
            pass
        if control_pass==0:
            subpart.append(bhits)
        else:
            pass

    if len(subpart)>0:
        besthits_filtered4[g]=subpart
    else:
        pass


if (option_cutoff=="on" and "chloroplast" in model) or (option_cutoff=="on" and "mitochondrion" in model):
    besthits_filtered5=dict()
    contigs_cov=[]
    contigs_name=[]
    cut_contigs=[]
    contigs_lgth=[]
    for g in list(besthits_filtered4.keys()):
        for bhits in besthits_filtered4[g]:
            cont=bhits[1]
            if cont in contigs_name:
                pass
            else:
                contigs_name.append(cont)
                cov=float(str(cont.split("cov_")[1]).split("_")[0])
                lgth=int(str(cont.split("length_")[1]).split("_")[0])
                contigs_cov.append(cov)
                contigs_lgth.append(lgth)

    if len(contigs_name)>1:
        bf=N.average(contigs_cov, weights=contigs_lgth)-3*weighted_std(contigs_cov, contigs_lgth)
        if bf<0:
            bf=N.average(contigs_cov, weights=contigs_lgth)-1*weighted_std(contigs_cov, contigs_lgth)
        else:
            pass
        if bf<0:
            bf=N.average(contigs_cov, weights=contigs_lgth)-1*mad(contigs_cov)
        else:
            pass

        for i in range(len(contigs_name)):
            if contigs_cov[i] < bf:
                cut_contigs.append(contigs_name[i])
            else:
                pass
    else:
        pass

    miss_genes=[]
    for g in list(besthits_filtered4.keys()):
        subpart=[]
        for bhits in besthits_filtered4[g]:
            seqid=bhits[1]
            control_pass=0
            if seqid in cut_contigs:
                control_pass=control_pass+1
                if g in miss_genes:
                    pass
                else:
                    miss_genes.append(g)
            else:
                pass
            if control_pass==0:
                subpart.append(bhits)
            else:
                pass
        if len(subpart)>0:
            besthits_filtered5[g]=subpart
        else:
            pass

    if len(cut_contigs)>0:
        print("WARN: %s contig(s) removed in %s according to standard deviation cut-off: %s" % (len(cut_contigs),str(nameofsample),",".join(cut_contigs)))
    else:
        pass
else:
    besthits_filtered5=besthits_filtered4

if len(list(besthits_filtered5.keys()))>0:
    print ("Extraction of sequences for %s" % nameofsample)
    inputs=range(len(list(besthits_filtered5.keys())))
    x=Parallel(n_jobs=num_cores)(delayed(GeneExtraction)(i) for i in inputs)

    l = [i for i in x if i is not None]

    outcontp=os.path.dirname(os.path.abspath(file))

    if len(l)==0:
        print("WARN: "+model+" - no genes following capture restrictions were extracted for %s" % nameofsample)
    else:
        pass

    if 'chloroplast' in model:
        if os.path.join(outcontp, str(nameofsample+".cont_cpdna.log")):
            with open(os.path.join(outcontp, str(nameofsample+".cont_cpdna.log")), 'a+') as outcontt:
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
            with open(os.path.join(outcontp, str(nameofsample+".cont_cpdna.log")), 'w') as outcontt:
                if len(l)>0:
                    for line in set(l):
                        outcontt.write(str(line)+'\n')
                else:
                    outcontt.write(str("None")+'\n')
    elif 'mitochondrion' in model:
        if os.path.join(outcontp, str(nameofsample+".cont_mtdna.log")):
            with open(os.path.join(outcontp, str(nameofsample+".cont_mtdna.log")), 'a+') as outcontt:
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
            with open(os.path.join(outcontp, str(nameofsample+".cont_mtdna.log")), 'w') as outcontt:
                if len(l)>0:
                    for line in set(l):
                        outcontt.write(str(line)+'\n')
                else:
                    outcontt.write(str("None")+'\n')
    else:
        pass
else:
    outcontp=os.path.dirname(os.path.abspath(file))
    print("WARN: "+model+" - no genes following capture restrictions were extracted for %s" % nameofsample)
    if 'chloroplast' in model:
        if os.path.join(outcontp, str(nameofsample+".cont_cpdna.log")):
            pass
        else:
            with open(os.path.join(outcontp, str(nameofsample+".cont_cpdna.log")), 'w') as outcontt:
                outcontt.write("None")
    elif 'mitochondrion' in model:
        if os.path.join(outcontp, str(nameofsample+".cont_mtdna.log")):
            pass
        else:
            with open(os.path.join(outcontp, str(nameofsample+".cont_mtdna.log")), 'w') as outcontt:
                outcontt.write("None")
    else:
        pass


print("--- %s seconds ---" % (time.time() - start_time))
