
<p align="center">
  <a>
    <img height="300" src="img/orthoskim_log.jpeg">
  </a>
</p>

OrthoSKim is a pipeline providing different tools to skim orthologous regions from whole genome low coverage sequencing data for nuclear, chloroplastic, ribosomal and mitochondrial compartments. This pipeline allow for region extracting from de novo targeted-assemblies using direct annotation when provided (in *chloroplast* and *rdnanuc* mode) or from wide-assemblies using mapping into references (in *nucleus* and *mitochondrion* mode).


### 1 . Installation
--------------------

OrthoSKim is tested on Unix environment and requires:
+ [Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)
+ [SPAdes](http://cab.spbu.ru/software/spades/)
+ [QUAST](https://github.com/ablab/quast)
+ [Diamond](https://github.com/bbuchfink/diamond)
+ [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
+ Needs Awk, Python

Some python libraries are also required:
+ Bio==0.1.0
+ joblib==0.13.2

### 2. Input files
------------------

OrthoSKim requires annotation files from sample assemblies of chloroplast and rdnanuc obtained with [ORG.Asm](https://git.metabarcoding.org/org-asm/org-asm) and [ORG.Annotate](https://git.metabarcoding.org/org-asm/org-annotate) softwares . This script has to be first run with *indexing* mode to identify samples required for skimming by searching them into **PATH_FIND_CHLORO/** directory.  


#### 2.1 - Config file (config_orthoskim.txt)


The following section describes the config file required for OrthoSKim. Users have to modify the *config_orthoskim.txt* file provided before running the pipeline. Default values are set for filtering and assembly steps.

`less config_orthoskim.txt`

```
# Global parameters :
TOOLS=/Users/pouchonc/Desktop/Scripts/OrthoSkim/tools.sh                             ## path to file with tools aliases
PATH_FIND_CHLORO=/Users/pouchonc/PhyloAlps/CDS/                                      ## path to find annotated chloroplast assemblies resulting from ORGasm
RES=/Users/pouchonc/PhyloAlps/run_orthoskim                                          ## path to directory to write output
PATHNAME_ASSEMBLY=Assembly                                                           ## name of assembly directory where contigs from SPAdes were put (before nucleus and mitochondrion modes)
EVALUE=0.001                                                                         ## evalue threshold for diamond steps
THREADS=15                                                                           ## Number of threads which will be used for diamond/SPAdes/QUAST softs
MINLENGTH=200                                                                        ## minimal length of alignment allowed mapping to reference
ANNOFMT=embl                                                                         ## format of annotated chloroplast/rdnanuc files
VERBOSE=0                                                                            ## Set verbose to TRUE (1) or FALSE (0)

# Preprocessing the data [indexing] mode
LIST_FILES=/Users/pouchonc/Desktop/Scripts/OrthoSkim/ressources/listSamples.tab      ## Parameters table for assembly path analysis according to annotated chloroplasts previously assembled. Specific format required: (1) chloroplast file name; (2) Genus of species sample\n-3: sample name\n-4: forward reads\n-5: reverse reads\n-6: Output path for assemblies samples subdirectories. This file could be obtained with the prep_assembly mode of OrthoSkim function

# Extraction steps from annotation files
# [chloroplast] mode :
CHLORO_GENES=/Users/pouchonc/Desktop/Scripts/OrthoSkim/ressources/listGenes.chloro   ## list of chloroplastic genes for extraction. Specific format of table: (1) type of gene [CDS,...]; (2) name of gene. This file could be modified by adding/removing specific lines.
# [nucrdna] mode :
NRDNA_GENES=/Users/pouchonc/Desktop/Scripts/OrthoSkim/ressources/listGenes.rdna      ## list of rdna nuclear genes for extraction. Specific format of table: (1) type of gene [CDS,...]; (2) name of gene. This file could be modified by adding/removing specific lines.

# Extraction steps from mapping assemblies into a reference
COVERAGE=5                                                                           ## Minimal contigs coverage allowed for genomic scan of mitochondrial and nuclear regions
MINCONTLENGTH=1000                                                                   ## Minimal contigs length allowed for genomic scan of mitochondrial and nuclear regions
# [SPAdes_assembly] mode:
MEMORY=30                                                                            ## Number of memory which will be used
KMER=55                                                                              ## Kmer fixed (here 55) or range values (as 21,33,55) for SPAdes assembly. Note: less than 128
# [nuclear] mode :
NUC_REF=/Users/pouchonc/Desktop/Scripts/OrthoSkim/ressources/refGenes.nu             ## list of nuclear genes of reference.  Amino acid sequence is specified. As the file contains bank of genes, gene name (header) has to be written following name_other-arguments (e.g. LFY_3702,LFY_3811 for LFY gene).
NUC_TYPE=exon                                                                        ## Type of structure extracted from the gff among "cds" (exon structure of CDS gene), "intron" (intron structure of CDS gene) and "all" (exon+intron)
# [mitochondrion] mode :
MITO_REF=/Users/pouchonc/Desktop/Scripts/OrthoSkim/ressources/refCDS.mito            ## list of mitochondrial coding genes of reference.  Amino acid sequence is specified. As the file contains bank of genes, gene name (header) has to be written following name_other-arguments (e.g. cox1_3702,cox1_3811 for cox1 gene).
MITO_TYPE=exon                                                                       ## Type of structure extracted from the gff among "cds" (exon structure of CDS gene), "intron" (intron structure of CDS gene) and "all" (exon+intron)
MITO_REF_RNA=/Users/pouchonc/Desktop/Scripts/OrthoSkim/ressources/refRNA.mito        ## list of mitochondrial non-coding regions of reference.  Nucleotidic sequence is specified. As the file contains bank of genes, gene name (header) has to be written following name_other-arguments (e.g. rrn18S_3702,rrn18S_3811 for rrn18S gene).

# Get mitochondrial sequences of reference [get_mitoRef] mode :
LENGTH_RATIO=0.75                                                                    ## length ratio of reference mitochondrial gene aligning the gene seed to consider this gene in final gene banks of reference.
SEEDS_MITO=/Users/pouchonc/Desktop/Scripts/OrthoSkim/ressources/seeds.mito           ## mitochondrial seeds of reference (from Arabidopsis_thaliana_3702_genbank) in amino acid sequence.
SEEDS_MITO_RNA=/Users/pouchonc/Desktop/Scripts/OrthoSkim/ressources/seedsRNA.mito           ## mitochondrial seeds of reference (from Arabidopsis_thaliana_3702_genbank) in amino acid sequence.
MITO_GENBANK_LOC=
```

#### 2.2 - Dependencies (tools.txt)

The access path of all dependencies required by OrthoSKim must be supplied in the *tools.sh* file, using following command:

`cat tools.sh`

```
#!/bin/bash

SPADES=/Users/pouchonc/PhyloAlps/OrthoSkim/TOOLS/SPAdes-3.13.0-Darwin/bin/spades.py
DIAMOND=/Users/pouchonc/miniconda2/bin/diamond
EXONERATE=/usr/local/bin/exonerate
QUAST=/Users/pouchonc/miniconda2/bin/quast.py
BLASTDB=/Users/pouchonc/miniconda2/bin/makeblastdb
BLASTN=/Users/pouchonc/miniconda2/bin/blastn
```

#### 2.3 - Annotation files


Annotation files must be supplied in **<LIST_FILES>** tab file.
This tab must contain for each sample:
+ the annotation file (format embl)
+ the genus named
+ the sample name
+ the forward reads
+ the reverse reads
+ pathname where additional assemblies are written (if *nucleus*/*mitochondrion* mode will be used)


`head ~/OrthoSkim/ressources/listSamples.tab`

```
/Users/pouchonc/PhyloAlps/CDS/Veronica_crassifolia:996476/Veronica_crassifolia:996476.CAR009639.BGN:NFI.chloro.embl     Veronica        Veronica_crassifolia:996476_CAR009639_BGN:NFI   /Users/pouchonc/PhyloAlps/CDS/Veronica_crassifolia:996476/BGN_NFIOSW_4_1_CA559ACXX.IND44_clean.fastq.gz /Users/pouchonc/PhyloAlps/CDS/Veronica_crassifolia:996476/BGN_NFIOSW_4_2_CA559ACXX.IND44_clean.fastq.gz /Users/pouchonc/PhyloAlps/run/Assembly/
/Users/pouchonc/PhyloAlps/CDS/Androsace_helvetica:199610/Androsace_helvetica:199610.CLA000520.BGN:ETA.chloro.embl       Androsace       Androsace_helvetica:199610_CLA000520_BGN:ETA    /Users/pouchonc/PhyloAlps/CDS/Androsace_helvetica:199610/BGN_ETAOSW_2_1_C8MR2ACXX.IND13_clean.fastq.gz  /Users/pouchonc/PhyloAlps/CDS/Androsace_helvetica:199610/BGN_ETAOSW_2_2_C8MR2ACXX.IND13_clean.fastq.gz  /Users/pouchonc/PhyloAlps/run/Assembly/
```


#### 2.4 - List of genes files

The extraction of orthologous regions for *chloroplast* and *nucrdna* modes is permitted by a given list of genes found in the annotation files. This list supplied in **<CHLORO_GENES>** and **<NRDNA_GENES>** must contain:
+ the type of gene (e.g. CDS,rRNA,tRNA)
+ the gene name

`head ~/OrthoSkim/ressources/listGenes.chloro`

```
tRNA    trnV
tRNA    trnA
tRNA    trnN
rRNA    rrn16S
rRNA    rrn23S
rRNA    rrn4.5S
rRNA    rrn5S
CDS     psbA
CDS     matK
CDS     rps16
CDS     psbK
```

By default, OrthoSkim provided tRNA, rRNA and CDS in chloroplastic list (see **<CHLORO_GENES>**). For the ribosomal complex, genes are rRNA (rrn18S, 5.8S rRNA, rrn28S) and misc_RNA (i.e. ITS1 and ITS2) (see **<NRDNA_GENES>**).

#### 2.5 - References files

The extraction of orthologous regions for *nucleus* and *mitochondrion* modes is permitted by mapping genomic assemblies into gene banks of references. For each genomic compartment, this bank is given in a multi-fasta file, containing amino acid sequences already align, with header written following name_other-arguments (e.g. cox1_3702,cox1_3811 for cox1 gene). These reference files can be modified to remove or add some genes with the specific format.

For mitochondria, the gene banks includes all CDS, and RNA (rRNA+tRNA) of available mitochondria aligned into *Arabidopsis_thaliana* seeds in **<SEEDS_MITO>** and **<SEEDS_MITO_RNA>**.

`head ~/OrthoSkim/ressources/refCDS.mito`

```
>cox2_Arabidopsis_thaliana_3702_genbank
MIVLKWLFLTISPCDAAEPWQLGSQDAATPIMQGIIDLHHDIFFFLILILVFVLWILVRALWHFHYKKNAIPQRIVHGTTIEILRTIFPSIISMFIAIPSFALLYSMDEVVVDPAITIKAIGHQWYRTYEYSDYNSSDEQSLTFDSYMIPEEDLELGQSRLLEVDNRVVVPAKTHLRIIVTSADVPHSWAVPSSGVKCDAVPGRLNQISILVQREGVYYGQCSEICGTNHAFTSIVVEAVPRKDYGSRVSNQLIPQTGEA
>cox2_Beta_macrocarpa_343494_genbank
MIVREWLFFTMAPCDAAEPWQLGFQDAATPMMQGIIDLHHDIFFFLILILVFVLWILVRALWHFHYKKNPIPQRIVHGTTIEIIWTIFPSIILMFIAIPSFALLYSMDEVVVDPAITIKAIGHQWYWTYEYSDYNSSDEQSLTFDSYMIPEDDLELGQ
LRLLEVDNRVVVPAKTHIRIIVTSADVLHSWAVPSLGVKCDAVPGRLNQTSILVQREGVYYGQCSEICGTNHAFMPIVVEAVSRKDYGSWVSNQLIPQTGEA
>cox2_Oryza_sativa_39946_genbank
XXXLECRFLTIALCDAAEPWQLGSQDAATPMMQGIIDLHHDIFFFLILILVFVSRMLVRALWHFNEQTNPIPQRIVHGTTIEIIRTIFPSVILLFIAIPSFALLYSMDGVLVDPAITIKAIGHQWYRTYEYSDYNSSDEQSLTFDSYTIPEDDPELGQ
SRLLEVDNRVVVPAKTHLRMIVTPADVPHSWAVPSSGVKCDAVPGRSNLTSISVQREGVYYGQCSEICGTNHAFTPIVVEAVTLKDYADWVSNQLILQTNXX
>cox3_Arabidopsis_thaliana_3702_genbank
MIESQRHSYHLVDPSPWPISGSLGALATTVGGVMYMHPFQGGARLLSLGLIFILYTMFVWWRDVLRESTLEGHHTKVVQLGPRYGSILFIVSEVMFFFAFFWASSHSSLAPAVEIGGIWPPKGIEVLDPWEIPFLNTPILPSSGAAVTWAHHAILAGKEKRAVYALVATVLLALVFTGFQGMEYYQAPFTISDSIYGSTFFLATGFHGFHVIIGTLFLIICGIRQYLGHLTKEHHVGFEAAAWYWHFVDVVWLFLFVSIYWWGGI
>cox3_Beta_macrocarpa_343494_genbank
MIESQRHSFHLVDPSPWPISGSLGALATTVGGVMYMHSFQGGATLLSLGLIFILYTMFVWWRDVLRESTLEGHHTKVVQLGLRYGFILFIVSEVMFFFAFFWAFFHSSLAPAVEIGGIWPPKGIWVLDPWEIPFLNTLILLSSGAAVTWAHHAILAGK
QKRAVYALVATVLLALVFTGFQGMEYYEAPFTISDSIYGSTFFLATGFHGFHVIIGTLFLIICGIRQYFGHLTKEHHVGFEAAAWYWHFVDVVWLFLFVSIYWWGGI
>cox3_Oryza_sativa_39946_genbank
MIESQRHSYHLVDPSPWPISGSLGALATTVGGVMYMHSFQGGATLLSLGLIFLLYTMFVWWRDVLRESTLEGHHTKAVQLGPRYGSILFIVSEVMFLFAFFWASSHSSLAPTVEIGGIWPPKGIGVLDPWEIPLLNTPILPSSGAAVTWAHHAILAGK
EKRAVYALVATVLLALVSTGFQGMEYYQAPSTISDSIYGSTFFLATGFHGFHVIIGTLFLIVCGIRQYLGHLTKKHHVGFEAAAWYWHFVDVVRLFPFVSIYWWGGI
>cox1_Arabidopsis_thaliana_3702_genbank
MKNLVRWLFSTNHKDIGTLYFIFGAIAGVMGTCFSVLIRMELARPGDQILGGNHQLYNVLITAHAFLMIFFMVMPAMIGGFGNWFVPILIGAPDMAFPRLNNISFWLLPPSLLLLLSSALVEVGSGTGWTVYPPLSGITSHSGGAVDLAIFSLHLSGVSSILGSINFITTIFNMRGPGMTMHRLPLFVWSVLVTAFLLLLSLPVLAGAITMLLTDRNFNTTFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHIVSTFSGKPVFGYLGMVYAMISIGVLGFLVWAHHMFTVGLDVDTRAYFTAATMIIAVPTGIKIFSWIATMWGGSIQYKTPMLFAVGFIFLFTIGGLTGIVLANSGLDIALHDTYYVVAHFHYVLSMGAVFALFAGFYYWVGKIFGRTYPETLGQIHFWITFFGVNLTFFPMHFLGLSGMPRRIPDYPDAYAGWNALSSFGSYISVVGICCFFVVVTITLSSGNNKRCAPSPWALELNSTTLEWMVQSPPAFHTFGELPAIKETKSYVK
>cox1_Beta_macrocarpa_343494_genbank
MTNLVRWLFSTNHKDIGTLYFIFGAIAGVMGTCFSVLIRMELARPGDQILGGNHQLYNVLITAHAFLMIFFMVMPAMIGGFGNWFVPILIGAPDMAFPRLNNISFWLLPPSLLLLLSSALVEVGSGTGWTVYPPLSGITSHSGGAVDLAIFSLHLSGV
SSILGSINFITTIFNMRGPGMTMHRLPLFVWSVLVTAFLLLLSLPVLAGAITMLLTDRNFNTTFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHIVSTFSGKPVFGYLGMVYAMISIGVLGFLVWAHHMFTVGLDVDTRAYFTAATMIIAV
PTGIKIFSWIATMWGGSIQYKTPMLFAVGFIFLFTVGGLTGIVLANSGLDIALHDTYYVVAHFHYVLSMGAVFALFAGFYYWVGKIFGRTYPETLGQIHFWITFFGVNLTFFPMHFLGLSGMPRRIPDYPDAYAGWNALSSFGSYISVVGICCFFVVV
TITLSSGKNKRCAPSPWAVEENSTTLEWMVQSPPAFHTFGELPAIKETKSXXX
>cox1_Oryza_sativa_39946_genbank
MTNLVRWLFSTNHKDIGTLYFIFGAIAGVMGTCFSVLIRMELARPGDQILGGNHQLYNVLITAHAFLMIFFMVMPAMIGGFGNWFVPILIGAPDMAFPRLNNISFWLLPPSLLLLLSSALVEVGSGTGWTVYPPLSGITSHSGGAVDLAIFSLHLSGV
SSILGSINFITTIFNMRGPGMTMHRLPLFVWSVLVTAFLLLLSLPVLAGAITMLLTDRNFNTTFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHIVSTFSRKPVFGYLGMVYAMISIGVLGFLVWAHHMFTVGLDVDTRAYFTAATMIIAVPTGIKIFSWIATMWGGSIQYKTPMLFAVGFIFLFTIGGLTGIVLANSGLDIALHDTYYVVAHFHYVLSMGAVFALFAGFYYWVGKIFGRTYPETLGQIHFWITFFGVNLTFFPMHFLGLSGMPRRIPDYPDAYAGWNALSSFGSYISVVGIRRFFVVVAITSSSGKNKRCAESPWAVEQNPTTLEWLVQSPPAFHTFGELPAIKETKSXXX
```

The mitochondrial list was produced by the *get_mitoRef* mode (see 3. Pipeline description below part).


The nucleus list of genes, provided with OrthoSkim, contains:
+ 1133* single copy orthologous genes ([BUSCO](https://busco.ezlab.org)) dataset for Viridiplantae (v.10) retrieved from Waterhouse et al. (2017), with 10 ancestral amino acid sequences for each gene
+ 352 ultra conserved element (UCE) designed for flowering plants and retrieved from Johnson et al. (2018).
+ additional nuclear plastid encoded and flowering-involved genes (LFY, FLC, FT, AG, ELF3,LATE)

\*Among the 1370 BUSCO of the initial dataset (which not mapped on chloroplasts and mitochondria), 237 genes mapped into the UCE dataset with diamond, and were removed to not included duplicates (here 1133 BUSCO in the file). UCE amino acid sequences were mapped into chloroplasts and mitochondria in order to keep nucleus compartment only.


`head ~/OrthoSkim/ressources/refGenes.nu`

```
>10018_0_BUSCO
IASVVSEIGLGSEPAFKVPEYDFRSPVDKLQKATGIPKAVFPVLGGLAVGLIALAYPEVLYWGFENVDILLESRPKGLSADLLLQLVAVKIVATSLCRASGLVGGYYAPSLFIGAATGMAYGKLILAEADPLLHLSILEVASPQAYGLVGMAATLAGVCQVPLTAVLLLFELTQDYRIVLPLLGAVGLSSWITSGQTKKELCKLESSLCLEDILVSEAMRTRYVTVLMSTLLVEAVSLMLAEKQSCALIVDEDNLLIGLLTLEDIQEFSKTVTPDMDLLSAEKIMGLSQLPVVVGLLDRECISL
>10018_1_BUSCO
VASVVSEIGLGSEPAFKVPEYDFRSAVDSLKKTLGLPKAVLPALGGLIVGLIALAYPEVLYWGFENVDILLESRPRGLSAELLLQLVAVKVVATSLCRASGLVGGYYAPSLFIGAATGMAYGKLIIAKADSLFDLEILEVASPQAYGLVGMAATLAGVCQVPLTAVLLLFELTQDYRIVLPLLGAVGLSSWISSKKTSKELCQLESSLCLKDVLVAEAMRTRYVTVLVTTSLAEALSLMLVEKQSLAVIVDEEDSLIGLLTLSDIQEYSKTVTPQLDLTKAEAIMELDRLAVVVGVLDRESIAL
>6500_UCE
SMQVVSALAVDHSGSRVLSGSYDYTVRMYDFQGMNSRLQSFRQLEPFEGHQVRSLSWSPTADRFLCVTGSAQAKIYDRDGLTLGEFVKGDMYIRDLKNTKGHISGLTCGEWHPKTKETILTSSEDGSLRIWDVNDFKSQKQVIKPKLARPGRVPVTTC
AWDREGKCIAGGIGDGSIQIWNLKPGWGSRPDIHVEKGHSDDITGLKFSSDGRILLSRSFDGSLKVWDLRQMKEPLKVFEDLPNHYAQTNIAFSPDEQLFLTGTSVERESTTGGLLCFYDRKLELVSRVGISPTCSVVQCAWHPKLNQIFATGDKQGG
THILYDPTLSERGALVCVARAPRKKSVDDFEAPVIHNPHALPLFRDQPSRKRQREKLKDPKSHKPELPITGPGGGRVGTKGSLLTQYLLKQGGLIKETWMEEDPREAILKYADVAAKDPKFIAPAYAQTQPEPVFAKSDSEDEE
>6506_UCE
LIKRRDVIGLGVSSLSAKGAGAALPPEKPRLCDDCEKELEKVPMVTTESGLQYKDIKVGGPSPPVGFQVAANYVAMPSGQIFDSSLEKGQFPYIFRVGSGQVIKGLDEGILSMKGGKRRLYIPGSLAFPKGLTSAPGRPRVAPNSPVIFDVSLEYIPG
LEVD
>LFY_sp|Q00958|LFY_ARATH
MDPEGFTSGLFRWNPTRALVQAPPPVPPPLQQQPVTPQTAAFGMRLGGLEGLFGPYGIRFYTAAKIAELGFTASTLVGMKDEELEEMMNSLSHIFRWELLVGERYGIKAAVRAERRRLQEEEEEESSRRRHLLLSAAGDSGTHHALDALSQEGLSEEPVQQQDQTDAAGNNGGGGSGYWDAGQGKMKKQQQQRRRKKPMLTSVETDEDVNEGEDDDGMDNGNGGSGLGTERQREHPFIVTEPGEVARGKKNGLDYLFHLYEQCREFLLQVQTIAKDRGEKCPTKVTNQVFRYAKKSGASYINKPKMRHYVHCYALHCLDEEASNALRRAFKERGENVGSWRQACYKPLVNIACRHGWDIDAVFNAHPRLSIWYVPTKLRQLCHLERNNAVAAAAALVGGISCTGSSTSGRGGCGGDDLRF
>PSBO_sp|P85194|PSBO_HELAN
MAASLQAAATFMTPTSRVQLKSSPSICKAFGIESTGSKVSCSLQADLKDFAQKCTDAAKIAGFALATSALVVSGASAEGSPKRLTYDEIQSKTYMEVKGTGTANQCPTIEGGVNGFAVKPGKYNAQKFCLEPTSFTVKAEGISKNSAPEFQKTKLMTRLTYTLDEIEGPLEVSSDGTIKFEEKDGIDYAAVTVQLPGGERVPFLFTIKELVATGKPESFGGNFLVPSYRGSSFLDPKGRGGSTGYDNAVALPAGGRGDEEELLKENIKNTAAGKGEITFSVTSSKPETGEVIGVFESIQPSDTDLGAKAPKDVKIQGVWYAQLW
>PSBO_sp|P23322|PSBO_SOLLC
QAAATLMQPTKVGVRNNLQLRSAQSVSKAFGVEQGSGRLTCSLQTEIKELAQKCTDAAKIAGFALATSALVVSGANAEGVPKRLTYDEIQSKTYMEVKGTGTANQCPTIEGGVGSFAFKPGKYTAKKFCLEPTSFTVKAEGVSKNSAPDFQKTKLMTRLTYTLDEIEGPFEVSPDGTVKFEEKDGIDYAAVTVQLPGGERVPFLFTIKQLVASGKPESFSGEFLVPSYRGSSFLDPKGRGGSTGYDNAVALPAGGRGDEEELQKENVKNTASLTGKITLSVTQSKPETGEVIGVFESIQPSDTDLGAKVPKDVKIQGIWYAQLX
```

This additional nuclear gene bank could be obtained following the code:

`nano list_add_nuclear_genes`

```
LFY_ARATH
FLC_ARATH
FT_ARATH
AG_ARATH
ELF3_ARATH
LATE_ARATH
PSBO_HELAN
PSBP_HELAN
PSBQ1_ARATH
PSBR_ARATH
PSBS_ARATH
PSBW_ARATH
PSBY_ARATH
PSAE1_ARATH
PSAF_ARATH
PSAG_ARATH
PSAH1_ARATH
PSAK_ARATH
PSAL_ARATH
PSAN_ARATH
PSAO_ARATH
ATPD_ARATH
ATPG1_ARATH
ATPX_SPIOL
```
```
$ awk '/^>/ {match($0,/RepID=([^ ]+)/,m1);match($0,/TaxID=([^ ]+)/,m2);split(m1[1],a,"=");split(m2[1],a2,"="); $0=">"a[1]" "a2[1]};{print}' ../DataBase/UniRef/UNIREF90_33090.fasta | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | perl -pe 's@>@@' | awk ' NR==FNR {a[$1]=$1;next} {if($1 in a) {split($1,spl,"_");print ">"spl[1]"_"$2"\n"$NF}}' list_add_nuclear_genes - > Add_Nuclear_genes.fa
$ diamond makedb --in Add_Nuclear_genes.fa -d db_add_genes
$ diamond blastp -d db_add_genes.dmnd -q ../DataBase/UniProt/UNIPROT_33090.fasta --outfmt 6 qseqid sseqid pident length mismatch gapopen qframe qstart qend sstart send evalue bitscore slen qseq --evalue 0.000001 -o matches_addgenes
$ awk '{print $1}' matches_addgenes | sort | uniq > hits_addgenes
$ awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' ../DataBase/UniProt/UNIPROT_33090.fasta  | perl -pe 's@>@@' | awk ' NR==FNR {a[$1]=$1;next} {if($1 in a) {print ">"$1"\n"$NF}}' hits_addgenes - > seq_hits_addgenes
$ exonerate --showcigar no --showtargetgff no --showvulgar no -q Add_Nuclear_genes.fa -t seq_hits_addgenes --ryo ">%qi gene=%qi ref=%ti length=%tl qalnlen=%qal qbal=%qab qeal=%qae qlen=%ql alnlen=%tal baln=%tab ealn=%tae score=%s\n%tas" --showalignment no | gawk '!/^Hostname:|^Command line:|^-- completed exonerate analysis/ {print $0}' > tmp.exonerate
$ python /Users/pouchonc/PhyloAlps/CDS/ExoRef.py -in tmp.exonerate -m inter -p ./ -r 0.75 -out refAddNuclearGenes
```



### 3. Pipeline description
---------------------------

OrthoSkim used different mode to index files, align and/or extract regions or check assemblies. Regardless of genomic compartment chosen mode, indexing has to be made in first. As chloroplast and nucrdna were targeted assembled, extraction can be directly processed from annotation files. For nucleus and mitochondrion extraction, a non-directional assembly has to be made with SPAdes on sample before running the pipeline in order to extract regions of interest (chloroplast regions could also be obtained from this mode by puting aligned regions into a new reference file).

**Note**: for the *nucleus*, *mitochondrion*, *nucrdna* and *chloroplast* modes, a *mode_done.log* file is created containing samples that were correctly processed, whereas unprocessed samples were added into *mode_error.log* file. This file could be used to remove processed samples from the initial **<LIST_FILES>** if the script has to be run. Command lines are also print if users want to rerun specific commands on samples.

#### 3.1 - Indexing files

This step allowed to guide the pipeline according to an indexing tab based on existing chloroplast assemblies. This tab is produced by *indexing* mode or is given in **<LIST_FILES>**. This allowed to screen each sample that will be used for extraction from the **<PATH_FIND_CHLORO>** directory. Unwanted samples must be removed from the list before processing other modes.

**Note**: If *nucleus* and/or *mitochondrion* mode will be used, *de novo* assembly has to be made on each sample using [SPAdes](http://cab.spbu.ru/software/spades/). To facilitate this, the output tab could be directly used to specify forward and reverse reads and outname of sample.


#### 3.2 - Extraction from annotation

For each sample specified in **<LIST_FILES>**, OrthoSkim will perform genes extraction directly from annotation according to a list of genes in *chloroplast* and *nucrdna* modes.

Results are output in **RES/** directory by creating subdirectories for each compartment and gene type, with a multifasta file per gene. For example, for chloroplastic CDS provided in **<CHLORO_GENES>**, OrthoSkim will output **RES/chloroplast_CDS/** subdirectory with CDS gene files.

`ls -l ~/RES/chloroplast_CDS/`

```
-rw-r--r--  1 pouchonc  staff   4874 16 avr 10:40 accD.fna
-rw-r--r--  1 pouchonc  staff   1875 16 avr 10:40 accD.pfna
-rw-r--r--  1 pouchonc  staff   4952 16 avr 10:40 atpA.fna
-rw-r--r--  1 pouchonc  staff   1901 16 avr 10:40 atpA.pfna
-rw-r--r--  1 pouchonc  staff   4853 16 avr 10:40 atpB.fna
-rw-r--r--  1 pouchonc  staff   1868 16 avr 10:40 atpB.pfna
-rw-r--r--  1 pouchonc  staff   1580 16 avr 10:40 atpE.fna
-rw-r--r--  1 pouchonc  staff    775 16 avr 10:40 atpE.pfna
-rw-r--r--  1 pouchonc  staff   2057 16 avr 10:40 atpF.fna
-rw-r--r--  1 pouchonc  staff    934 16 avr 10:40 atpF.pfna
```

**Note**: As above, if gene type is CDS, Othoskim will output both nucleotidic sequences in .fna and amino acid sequences in .pfna. For other type, only nucleotidic sequences are written.



#### 3.3 - Extraction from mapping to reference

The study of nuclear and mitochondrial compartments has to be driven on new genomic assemblies performed with [SPAdes](http://cab.spbu.ru/software/spades/) because of available assemblies for sample was conducted only on chloroplast or rDNA cluster with OrgAsm. This allowed to assembly all the genomic compartment and after extract interested regions instead of assembly specific chloroplast/rDNA.


##### 3.3.1 - SPAdes assembly run

[SPAdes](http://cab.spbu.ru/software/spades/) assembly could be run using the *SPAdes_assembly* mode from the **<LIST_FILES>** generated with *indexing* mode, by accessing to the forward and reverse reads, and by keeping the sample name providing in the file. [SPAdes](http://cab.spbu.ru/software/spades/) will be run using **THREADS** and **KMER** specified in the *config_orthoskim.txt* file.

Othoskim will output then a **samplename/** subdirectory into the **PATHNAME_ASSEMBLY** given per sample included in the **<LIST_FILES>**.  

##### 3.3.2 - Preprocessing mode

After [SPAdes](http://cab.spbu.ru/software/spades/) runs, OrthoSkim has first to preprocess SPAdes scaffolding contigs by renaming the file according to the same sample name provided in **<LIST_FILES>** and ordering them into **${RES}/${PATHNAME_ASSEMBLY}/Samples/** directory.
This is made under *SPAdes_reformate* mode.

##### 3.3.3 - Mitochondrial gene bank (for *mitochondrion* mode)

In order to extract mitochondrial CDS for our samples, a reference database must be created. To do this, annotation files for mitochondria from genbank has to be download and put into the **MITO_GENBANK_LOC/** directory.

For all these file, OrthoSkim in *get_mitoRef* mode will extract all notified CDS and align them into the CDS of *Arabidopsis_thaliana* as seeds thanks to *exonerate* (proteins *versus* proteins). Output **<MITO_REF>** file is created containing a bank of CDS genes, all well identified.
A similar approach made for non-coding regions RNAs (rRNA and tRNA with nucleotidic sequences) will come soon (a reference list of RNAs is actually provided in the **ressources/** folder).


##### 3.3.4 - Alignment extraction

###### 3.3.4.1 Coding Regions

After this, in both *mitochondrion_CDS* and *nucleus* mode, a [diamond](https://github.com/bbuchfink/diamond) database is created for each amino acid sequences provided in **<MITO_REF>** and **<NUC_REF>** (with *diamond makedb*).

Mapping is made for each sample with [diamond](https://github.com/bbuchfink/diamond) *blastx* to quickly identify contigs hits (only hits with a minimal **EVALUE** were retained).
Alignment was then conducted on these contigs from [exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate) by using the *protein2genome* mode incorporating all the appropriate gaps and frameshifts, and modelling of introns. Only contigs with a minimal coverage of **COVERAGE** and length of **MINCONTLENGTH** were considered for alignments into reference. Next, extraction was conducted through python script into [exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate) .gff output table. Type of gene structure (exon, intron or all:exon+intron) required during extraction from mapping are reported in **NUC_TYPE** or **MITO_TYPE**. This step is conducted into multiple processors using the **THREADS** specified in the the *config_orthoskim.txt* file.

###### 3.3.4.2 Non-Coding Regions

Non-coding regions (here mitochondrial RNAs) were retrieved in *mitochondrion_RNA* mode using a same approach but involving *nucleotide versus nucleotide* searches through: a [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) database compilation (*makeblastdb* program) with **<MITO_REF_RNA>** references, hits identification with *blastn* command and Alignment+Extraction from [exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate) in *genome2genome* mode.


All the output files will be in the **${RES}/mitochondrion/** and **${RES}/nucleus/** directories, as following:

`ls -l ~/RES/nucleus/`

```
-rw-r--r--  1 pouchonc  staff  1758  5 jui 11:11 10104_BUSCO.fna
-rw-r--r--  1 pouchonc  staff  1964  5 jui 11:11 10521_BUSCO.fna
-rw-r--r--  1 pouchonc  staff  5071  5 jui 11:11 10785_BUSCO.fna
-rw-r--r--  1 pouchonc  staff  1400  5 jui 11:11 11487_BUSCO.fna
-rw-r--r--  1 pouchonc  staff  2040  5 jui 11:11 11505_BUSCO.fna
-rw-r--r--  1 pouchonc  staff  1778  5 jui 11:11 1504_BUSCO.fna
```

#### 3.4 - Summary statistics of assemblies

OrthoSkim allowed to output summary statistic on contigs assemblies for both *OrgAsm* (chloroplast+rdna) and [SPAdes](http://cab.spbu.ru/software/spades/) runs thanks to [QUAST](https://github.com/ablab/quast) by specifying the *stat_[chloro,rdna,SPAdes]* modes.

The output *transposed_report.txt* tab file will be in **${RES}/report_[chloro,rdnanuc,SPAdes]_assemblies/** directories given indication on the assembly. For example, for the chloroplast, we except the assembly of a single contig sizing 150,000 Nt in average as in the following example.  

`head ~/RES/report_chloro_assemblies/transposed_report.txt`

```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                                              # contigs (>= 0 bp)  # contigs (>= 1000 bp)  # contigs (>= 5000 bp)  # contigs (>= 10000 bp)  # contigs (>= 25000 bp)  # contigs (>= 50000 bp)  Total length (>= 0 bp)  Total length (>= 1000 bp)  Total length (>= 5000 bp)  Total length (>= 10000 bp)  Total length (>= 25000 bp)  Total length (>= 50000 bp)  # contigs  Largest contig  Total length  GC (%)  N50     N75     L50  L75  # N's per 100 kbp
Veronica_crassifolia_996476.CAR009639.BGN_NFI.chloro  1                    1                       1                       1                        1                        1                        152361                  152361                     152361                     152361                      152361                      152361                      1          152361          152361        38.03   152361  152361  1    1    0.00             
Androsace_helvetica_199610.CLA000520.BGN_ETA.chloro   1                    1                       1                       1                        1                        1                        147550                  147550                     147550                     147550                      147550                      147550                      1          147550          147550        37.06   147550  147550  1    1    0.00             
Doronicum_columnae_118758.PHA003018.BGN_EEH.chloro    1                    1                       1                       1                        1                        1                        152623                  152623                     152623                     152623                      152623                      152623                      1          152623          152623        37.60   152623  152623  1    1    0.00             
```

OrthoSkim will also output the *report.pdf* file generated with [QUAST](https://github.com/ablab/quast) containing:

<figure align="center">
  <img height="240" src="img/Nx.jpeg" title="Nx values varying from 0 to 100%" />
  <figcaption>Nx values varying from 0 to 100%</figcaption>
</figure>
<br>
<figure align="center">
  <img height="240" src="img/GCall.jpeg"/>
  <figcaption>*GC content in the contigs for all samples*</figcaption>
</figure>
<figure align="center">
  <img height="240" src="img/GCsamp1.jpeg"/>
  <figcaption>*GC content in the contigs for sample1*</figcaption>
</figure>
<figure align="center">
  <img height="240" src="img/GCsamp2.jpeg"/>
  <figcaption>*GC content in the contigs for sample2*</figcaption>
</figure>
<figure align="center">
  <img height="240" src="img/GCsamp3.jpeg"/>
  <figcaption>*GC content in the contigs for sample3*</figcaption>
</figure>
<br>
<figure align="center">
  <img height="240" src="img/Covall.jpeg"/>
  <figcaption>*distribution of total contig lengths at different*</figcaption>
  <figcaption>*read coverage (only for SPAdes mode)*<figcaption>         
</figure>
<figure align="center">
  <img height="240" src="img/Covsamp1.jpeg"/>
  <figcaption>*coverage distribution for sample1*</figcaption>
</figure>
<figure align="center">
  <img height="240" src="img/Covsamp2.jpeg"/>
  <figcaption>*coverage distribution for sample2*</figcaption>
</figure>
<figure align="center">
  <img height="240" src="img/Covsamp3.jpeg"/>
  <figcaption>*coverage distribution for sample3*</figcaption>
</figure>
<br>


**Note**: see QUAST [manual](http://quast.bioinf.spbau.ru/manual.html) for more details. Outputs for [icarus](http://bioinf.spbau.ru/icarus) genome visualizer were also kept in the directory to visualize assemblies.

OrthoSkim allowed also to get statistic from gene extraction by using the *get_stat* mode for sequences (\*.fna) found in given **-p PATH**. The pipeline output a *report.tab* into this path containing:
+ gene name (gene_name)
+ samples found for this gene (#taxa)
+ the biggest length of sequence found (maxlen)
+ the number of sample covering this length at 100% (#1.0), 75% (#0.75), 50% (#0.5) and 25% (#0.25).

`head ~/PATH/report.log`

```
gene_name       #taxa   maxlen  #1.0    #0.75   #0.5    #0.25
psaB    3       2205    3       3       3       3
psbN    3       132     3       3       3       3
rps11   3       417     3       3       3       3
ndhA    3       1083    3       3       3       3
rpl14   3       369     3       3       3       3
rps7    3       468     3       3       3       3
```
<br>

### 4. Running OrthoSkim
------------------

OrthoSkim could be called step by step. Recommendations about steps are given in the previous description (section 3). After edition of the *tools.sh* and *config_orthoskim.txt* files, the first step of *indexing* should be called as following:

```
./OrtoSkim_v0.0.sh -m indexing -c config_orthoskim.txt
```

To run *chloroplast* and *nucrdna* modes :

```
./OrtoSkim_v0.0.sh -m chloroplast -c config_orthoskim.txt
```
```
./OrtoSkim_v0.0.sh -m nucrdna -c config_orthoskim.txt
```

For *mitochondrion* and *nucleus* mode, the script should be called as following:

```
./OrtoSkim_v0.0.sh -m SPAdes_assembly -c config_orthoskim.txt
./OrtoSkim_v0.0.sh -m SPAdes_reformate -c config_orthoskim.txt
```


```
./OrtoSkim_v0.0.sh -m get_mitoRef -c config_orthoskim.txt
./OrtoSkim_v0.0.sh -m mitochondrion_CDS -c config_orthoskim.txt
./OrtoSkim_v0.0.sh -m mitochondrion_RNA -c config_orthoskim.txt
```
```
./OrtoSkim_v0.0.sh -m nucleus -c config_orthoskim.txt
```
**Note**: *SPAdes_assembly* and *SPAdes_reformate* have to be run just once for both *mitochondrion* and *nucleus* mode. For mitochondrion mode, the *get_mitoRef* mode was not necessary ran if reference is already provided in the *config_orthoskim.txt* file.

For summary statistic, OrthoSkim could be run using the  `-m [stat_chloro,stat_rdna,stat_SPAdes]`, or using the `-m get_stat` with `-p PATH` option as following:

```
./OrtoSkim_v0.0.sh -m stat_SPAdes -c config_orthoskim.txt
./OrtoSkim_v0.0.sh -m get_stat -c config_orthoskim.txt -p path_to_extracted_files
```



### 5. References
--------------------
+ Bankevich, A., Nurk, S., Antipov, D., Gurevich, A.A., Dvorkin, M., Kulikov, A.S., et al. (2012). SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. J. Comput. Biol., 19, 455–477.
+ Buchfink, B., Xie, C. & Huson, D.H. (2015). Fast and sensitive protein alignment using DIAMOND. Nature Methods, 12, 59–60.
+ Gurevich, A., Saveliev, V., Vyahhi, N. & Tesler, G. (2013). QUAST: quality assessment tool for genome assemblies. Bioinformatics, 29, 1072–1075.
+ Johnson, M.G., Pokorny, L., Dodsworth, S., Botigué, L.R., Cowan, R.S., Devault, A., et al. (n.d.). A Universal Probe Set for Targeted Sequencing of 353 Nuclear Genes from Any Flowering Plant Designed Using k-Medoids Clustering. Syst Biol.
+ Simão, F.A., Waterhouse, R.M., Ioannidis, P., Kriventseva, E.V. & Zdobnov, E.M. (2015). BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. Bioinformatics, 31, 3210–3212.
+ Slater, G.S.C. & Birney, E. (2005). Automated generation of heuristics for biological sequence comparison. BMC Bioinformatics, 6, 31.
+ Waterhouse, R.M., Seppey, M., Simão, F.A., Manni, M., Ioannidis, P., Klioutchnikov, G., et al. (2018). BUSCO Applications from Quality Assessments to Gene Prediction and Phylogenomics. Molecular Biology and Evolution, 35, 543–548.
