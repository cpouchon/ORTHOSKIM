
<p align="center">
  <a>
    <img height="300" src="img/orthoskim_log.jpeg">
  </a>
</p>

OrthoSKim is a pipeline providing different tools to skim orthologous regions from whole genome low coverage sequencing data for nuclear, chloroplastic, ribosomal and mitochondrial compartments. This pipeline allow for region extracting from de novo targeted-assemblies using direct annotation when provided (in *extraction_targeted* mode) or from wide-assemblies using mapping into references (in *extraction_untargeted* mode).


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
+ pandas==0.24.2
+ ete2==2.3.10
+ numpy==1.16.2

### 2. Input files
------------------

OrthoSKim can used annotation files from sample assemblies of chloroplast and rdnanuc obtained with [ORG.Asm](https://git.metabarcoding.org/org-asm/org-asm) and [ORG.Annotate](https://git.metabarcoding.org/org-asm/org-annotate) softwares with *extraction_targeted* mode. Otherwise, OrthoSkim could be performed directly wide assemblies (in *SPAdes_assembly* mode) from the sequencing reads. In both cases, a sample files and a config files are required (see 2.3 section), including reference sequences, and list of targeted genes (see 2.4 and 2.5 sections).


#### 2.1 - Config file (config_orthoskim.txt)


The following section describes the config file required for OrthoSKim. Users have to modify the *config_orthoskim.txt* file provided before running the pipeline. Default values are set for filtering and assembly steps.

```
less config_orthoskim.txt
```

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

# [indexing] mode: preprocessing the data
LIST_FILES=/Users/pouchonc/Desktop/Scripts/OrthoSkim/ressources/listSamples.tab      ## Parameters table for assembly path analysis according to annotated chloroplasts previously assembled. Specific format required: (1) chloroplast file name; (2) Genus of species sample\n-3: sample name\n-4: forward reads\n-5: reverse reads\n-6: Output path for assemblies samples subdirectories. This file could be obtained with the prep_assembly mode of OrthoSkim function

# [extraction_targeted] mode: extraction steps from annotation files
###### [chloroplast] target:
CHLORO_GENES=/Users/pouchonc/Desktop/Scripts/OrthoSkim/ressources/listGenes.chloro   ## list of chloroplastic genes for extraction. Specific format of table: (1) type of gene [CDS,...]; (2) name of gene. This file could be modified by adding/removing specific lines.
###### [nucrdna] target:
NRDNA_GENES=/Users/pouchonc/Desktop/Scripts/OrthoSkim/ressources/listGenes.rdna      ## list of rdna nuclear genes for extraction. Specific format of table: (1) type of gene [CDS,...]; (2) name of gene. This file could be modified by adding/removing specific lines.

# [extraction_untargeted]: extraction steps from mapping assemblies into a reference
COVERAGE=5                                                                           ## Minimal contigs coverage allowed for genomic scan of mitochondrial and nuclear regions
MINCONTLENGTH=1000                                                                   ## Minimal contigs length allowed for genomic scan of mitochondrial and nuclear regions
MODE_REF=distance                                                                    ## Mode to select the closest reference from the database [taxaonomy,distance]
DISTANCE_MATRIX=                                                                     ## distance matrix at genus level
TAXONOMY_UPDATE=no                                                                   ## update of NCBI taxonomy DB [no,yes]
###### [SPAdes_assembly] mode:
MEMORY=30                                                                            ## Number of memory which will be used
KMER=55                                                                              ## Kmer fixed (here 55) or range values (as 21,33,55) for SPAdes assembly. Note: less than 128
###### [nuclear] target :
NUC_REF=~/OrthoSkim-master/ressources/refGenes.nu                                    ## list of nuclear genes of reference.  Amino acid sequence is specified. As the file contains bank of genes, gene name (header) has to be written following name_other-arguments (e.g. LFY_3702,LFY_3811 for LFY gene).
NUC_TYPE=exon                                                                        ## Type of structure extracted from the gff among "cds" (exon structure of CDS gene), "intron" (intron structure of CDS gene) and "all" (exon+intron)
NUC_DB_TYPE=orthoskim                                                                ## source of nuclear genes [orthoskim,personal]. If personal option is set, please see the documentation to see the required format of the DB.
###### [mitochondrion] target :
MITO_GENES=listGenes.mito                                                            ## list of mitochondrial genes that will be processed (CDS+rRNA+tRNA)
SEEDS_MITO_CDS=~/OrthoSkim-master/ressources/mitoCDS.seeds                           ## CDS mitochondrial seeds of reference (from Arabidopsis_thaliana_3702_genbank), proteic sequence.
SEEDS_MITO_rRNA=~/OrthoSkim-master/ressources/mitorRNA.seeds                         ## rRNA mitochondrial seeds of reference (from Arabidopsis_thaliana_3702_genbank), nucleotidic sequence.
SEEDS_MITO_tRNA=~/OrthoSkim-master/ressources/mitotRNA.seeds                         ## tRNA mitochondrial seeds of reference (from Arabidopsis_thaliana_3702_genbank), nucleotidic sequence.
MITO_REF_CDS=~/OrthoSkim-master/ressources/mit_CDS_unaligned.fa                      ## list of mitochondrial coding genes of reference.  Amino acid sequence is specified. As the file contains bank of genes, gene name (header) has to be written following name_other-arguments (e.g. cox1_3702_Genus_species,cox1_3811_Genus_species for cox1 gene).
MITO_REF_tRNA=~/OrthoSkim-master/ressources/mit_tRNA_unaligned.fa                    ## list of mitochondrial tRNA non-coding regions of reference.  Nucleotidic sequence is specified. As the file contains bank of genes, gene name (header) has to be written following name_other-arguments (e.g. trnI_3702_Genus_species,trnI_3811_Genus_species for trnI gene).
MITO_REF_rRNA=~/OrthoSkim-master/ressources/mit_rRNA_unaligned.fa                    ## list of mitochondrial rRNA non-coding regions of reference.  Nucleotidic sequence is specified. As the file contains bank of genes, gene name (header) has to be written following name_other-arguments (e.g. rrn18S_3702_Genus_species,rrn18S_3811_Genus_species for rrn18S gene).
MITO_TYPE=exon                                                                       ## Type of structure extracted from the gff among "exon" (exon structure of CDS gene), "intron" (intron structure of CDS gene) and "all" (exon+intron)
###### [chloroplast] target:
CHLORO_GENES=~/OrthoSkim-master/ressources/listGenes.chloro                          ## list of chloroplastic genes that will be processed (CDS+rRNA+tRNA)
SEEDS_CHLORO_CDS=~/OrthoSkim-master/ressources/chloroCDS.seeds                       ## chloroplastic CDS seeds of reference (from Arabidopsis_thaliana_3702_genbank), proteic sequence
SEEDS_CHLORO_rRNA=~/OrthoSkim-master/ressources/chlororRNA.seeds                     ## chloroplastic rRNA seeds of reference (from Arabidopsis_thaliana_3702_genbank), nucleotidic sequence.
SEEDS_CHLORO_tRNA=~/OrthoSkim-master/ressources/chlorotRNA.seeds                     ## chloroplastic tRNA seeds of reference (from Arabidopsis_thaliana_3702_genbank), nucleotidic sequence.
CHLORO_REF_CDS=~/OrthoSkim-master/ressources/chloro_CDS_unaligned.fa
CHLORO_REF_rRNA=~/OrthoSkim-master/ressources/chloro_rRNA_unaligned.fa
CHLORO_REF_tRNA=~/OrthoSkim-master/ressources/chloro_tRNA_unaligned.fa
CHLORO_TYPE=exon
###### [nucrdna] target:
NRDNA_REF_misc_RNA=~/OrthoSkim-master/ressources/
NRDNA_REF_rRNA=~/OrthoSkim-master/ressources/

# [DB_mitochondrion] mode: mitochondrial sequences of reference
MITO_GENBANK=~/OrthoSkim-master/ressources/                                          ## file with mitochondrial annotations from genbank

# Selection of taxa from extractions
TAXALIST=~/OrthoSkim-master/ressources/selTaxa_Primulaceae.tab                       ## list of taxa to select from extractions
EXTENSION=fa     
```

#### 2.2 - Dependencies (tools.txt)

The access path of all dependencies required by OrthoSKim must be supplied in the *tools.sh* file, using following command:

```
cat tools.sh
```
```
#!/bin/bash

SPADES=/Users/pouchonc/PhyloAlps/OrthoSkim/TOOLS/SPAdes-3.13.0-Darwin/bin/spades.py
DIAMOND=/Users/pouchonc/miniconda2/bin/diamond
EXONERATE=/usr/local/bin/exonerate
QUAST=/Users/pouchonc/miniconda2/bin/quast.py
BLASTDB=/Users/pouchonc/miniconda2/bin/makeblastdb
BLASTN=/Users/pouchonc/miniconda2/bin/blastn
```

#### 2.3 - Samples file


Samples file must be supplied in **<LIST_FILES>** tab file.
This tab must contain for each sample:
+ the pathfile to annotation file (format embl)
+ the genus named
+ the sample name with Genus_species_taxid_others order
+ the pathfile to forward reads
+ the pathfile reverse reads

**Note**: if annotations files are missing, the first column requires a "NA" value.

```
head ~/OrthoSkim/ressources/listSamples.tab
```

```
/Users/pouchonc/PhyloAlps/CDS/Veronica_crassifolia:996476/Veronica_crassifolia:996476.CAR009639.BGN:NFI.chloro.embl     Veronica        Veronica_crassifolia_996476_CAR009639_BGN_NFI   /Users/pouchonc/PhyloAlps/CDS/Veronica_crassifolia:996476/BGN_NFIOSW_4_1_CA559ACXX.IND44_clean.fastq.gz /Users/pouchonc/PhyloAlps/CDS/Veronica_crassifolia:996476/BGN_NFIOSW_4_2_CA559ACXX.IND44_clean.fastq.gz
/Users/pouchonc/PhyloAlps/CDS/Androsace_helvetica:199610/Androsace_helvetica:199610.CLA000520.BGN:ETA.chloro.embl       Androsace       Androsace_helvetica_199610_CLA000520_BGN_ETA    /Users/pouchonc/PhyloAlps/CDS/Androsace_helvetica:199610/BGN_ETAOSW_2_1_C8MR2ACXX.IND13_clean.fastq.gz  /Users/pouchonc/PhyloAlps/CDS/Androsace_helvetica:199610/BGN_ETAOSW_2_2_C8MR2ACXX.IND13_clean.fastq.gz
```


#### 2.4 - List of genes files

The extraction of orthologous regions for *chloroplast*, *mitochondrion* and *nucrdna* targets is permitted by a given list of genes. This list supplied in **<CHLORO_GENES>**, **<MITO_GENES>** and **<NRDNA_GENES>** must contain:
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

By default, OrthoSkim provided tRNA, rRNA and CDS in chloroplastic and mitochondrial list (see **<CHLORO_GENES>** and **<MITO_GENES>**). For the ribosomal complex, genes are rRNA (rrn18S, 5.8S rRNA, rrn28S) and misc_RNA (i.e. ITS1 and ITS2) (see **<NRDNA_GENES>**).

#### 2.5 - References files (database)

The extraction of orthologous regions for *nucleus*, *chloroplast* and *mitochondrion* target is permitted by mapping genomic assemblies into gene banks of references. For each genomic compartment, this bank is given in a multi-fasta file, containing amino acid sequences with header written following genename_taxid_Genus_species_other-arguments (e.g. cox1_3702,cox1_3811 for cox1 gene). Location of these files are set in the *config_orthoskim.txt* file at **<[CHLORO,MITO]_REF_[CDS,tRNA,rRNA]>**.

OrthoSkim requires also **<SEEDS_MITO_CDS>**, **<SEEDS_MITO_rRNA>**, **<SEEDS_MITO_tRNA>**, **<SEEDS_CHLORO_CDS>**, **<SEEDS_CHLORO_rRNA>** and **<SEEDS_CHLORO_tRNA>** containing the sequences of a model organism (*Arabidopsis_thaliana* seeds are supplied in OrthoSkim for the chloroplast and the mitochondrion).

For mitochondria, a reference bank of genes is supplied in OrthoSkim including all CDS, and RNA (rRNA+tRNA) of available mitochondria from the [NCBI](https://www.ncbi.nlm.nih.gov/genbank/), which are mapped into the seeds in **<SEEDS_MITO_CDS>**, **<SEEDS_MITO_rRNA>** and **<SEEDS_MITO_tRNA>** to keep a standard gene name.

```
head ~/OrthoSkim/ressources/mit_CDS_unaligned.fa
```
```
>cox2_103999_Codonopsis_lanceolata
MRELEKKNTHDFILPAPADAAEPWQLGFQDGATPIMQGIIDLHHDIFFFLIMILVLVLWILVRALWLFSSKRNPIPQRIVHGTTIEILRTIFPSIILMFIAIPSFALLYSMDEVVVDPAITIKAIGHQWYWTYEYSDYNSSDEESLTFDSYMIPEDDLELGQLRLLEVDNRVVVPANCHLRLIVTSADVPHSWAVPSLGVKCDAVPGRLNQVSISVLREGVYYGQCSEICGTNHAFMPIVVEAVSMKDYASRVSNQLIPQTGH
>cox2_104537_Roya_obtusa
MILKSLFQVVYCDAAEPWQLGFQDAATPMMQGIIDLHHDIMFFITIIITFVLWMLVRVLWHFHYKKNPIPQRFVHGTTIEIIWTIIPSIILMFIAIPSFALLYSMDEVVDPAITIKAIGHQWYWSYEYSDYSTSDEESLAFDSYMIPEDDLELGQLRLLEVDNRVVVPAKTHLRFIITSADVLHSWAVPSLGVKCDAVPGRLNQTSIFIKREGVYYGQCSEICGTNHAFMPIVVEAVSLDDYVSWVSNKME
>cox2_111617_Ulva_fasciata
MKNFSFSYCILITLFNISVISSCDAPLSATSAMLDRFGFQEPASPLMEGLIALHSDIWAIMLFVAGFVLYMMCAILYNFSASSSEISYKVHHHSLIEIVWTTIPALILCVIAIPSFTLLYSLDEVIEPSLTIKAIGRQWYWSYEYGDYEVHDGLITNGITFDSNVLQDDDLEQGQLRLLDVDNRLVLPVNRHIRLLTSGGDVIHSFAVPSLGVKLDAIPGRLNQTMVFIKRQGVFYGQCSELCGSSHGMMPIALEAVREQDYVDWVNIKLQEM
>cox1_112509_Hordeum_vulgare_subsp._vulgare
MTNLVRWLFSTNHKDIGTLYFIFGAIAGVMGTCFSVLIRMELARPGDQILGGNHQLYNVLITAHAFLMIFFMVMPAMIGGFGNWFVPILIGAPDMAFPRLNNISFWLLPPSLLLLLSSALVEVGSGTGWTVYPPLSGITSHSGGAVDLAIFSLHLSGISSILGSINFITTIFNMRGPGMTMHRLPLFVWSVLVTAFLLLLSLPVLAGAITMLLTDRNFNTTFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHIVSTFSRKPVFGYLGMVYAMISIGVLGFLVWAHHMFTVGLDVDTRAYFTAATMIIAVPTGIKIFSWIATMWGGSIQYKTPMLFAVGFIFLFTIGGLTGIVLANSGLDIALHDTYYVVAHFHYVLSMGAVFALFAGFYYWVGKIFGRTYPETLGQIHFWITFFGVNLTFFPMHFLGLSGMPRRIPDYPDAYAGWNALSSFGSYISVVGIRRFFVVVAITSSSGKNKKCAESPWAVEQNPTTLEWLVQSPPAFHTFGELPAVKETKNLS
>cox1_113636_Populus_tremula
MINLVRWLFSTNHKDIGTLYFIFGAIAGVMGTCFSVLIRMELARPGDQILGGNHQLYNVLITAHAFLMIFFMVMPAMIGGFGNWFVPILIGAPDMAFPRLNNISFWLLPPSLLLLLSSALVEVGSGTGWTVYPPLSGITSHSGGAVDLAIFSLHLSGVSSILGSINFITTIFNMRGPGMTMHRLPLFVWSVLVTAFLLLLSLPVLAGAITMLLTDRNFNTTFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHIVSTFSGKPVFGYLGMVYAMISIGVLGFLVWAHHMFTVGLDVDTRAYFTAATMIIAVPTGIKIFSWIATMWGGSIQYKTPMLFAVGFIFLFTIGGLTGIVLANSGLDIALHDTYYVVAHFHYVLSMGAVFALFAGFYYWVGKITGRTYPETLGKIHFWITFFGVNLTFFPMHFLGLSGMPRRIPDYPDAYAGWNALSSFGSYISVVGICCFFVVVTITLSSGNQNKCAPSPWALEQNSTTLEWMVQSPPAFHTFGELPAIKETKSYVK
MKNLVRWLFSTNHKDIGTLYFIFGAIAGVMGTCFSVLIRMELARPGDQILGGNHQLYNVLITAHAFLMIFFMVMPAMIGGFGNWFVPILIGAPDMAFPRLNNISFWLLPPSLLLLLSSALVEVGSGTGWTVYPPLSGITSHSGGAVDLAIFSLHLSGVSSILGSINFITTIFNMRGPGMTMHRLPLFVWSVLVTAFLLLLSLPVLAGAITMLLTDRNFNTTFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHIVSTFSGKPVFGYLGMVYAMISIGVLGFLVWAHHMFTVGLDVDTRAYFTAATMIIAVPTGIKIFSWIATMWGGSIQYKTPMLFAVGFIFLFTIGGLTGIVLANSGLDIALHDTYYVVAHFHYVLSMGAVFALFAGFYYWVGKIFGRTYPETLGQIHFWITFFGVNLTFFPMHFLGLSGMPRRIPDYPDAYAGWNALSSFGSYISVVGICCFFVVVTITLSSGNNKRCAPSPWALELNSTTLEWMVQSPPAFHTFGELPAIKETKSYVK
>cox1_Beta_macrocarpa_343494_genbank
MTNLVRWLFSTNHKDIGTLYFIFGAIAGVMGTCFSVLIRMELARPGDQILGGNHQLYNVLITAHAFLMIFFMVMPAMIGGFGNWFVPILIGAPDMAFPRLNNISFWLLPPSLLLLLSSALVEVGSGTGWTVYPPLSGITSHSGGAVDLAIFSLHLSGV
SSILGSINFITTIFNMRGPGMTMHRLPLFVWSVLVTAFLLLLSLPVLAGAITMLLTDRNFNTTFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHIVSTFSGKPVFGYLGMVYAMISIGVLGFLVWAHHMFTVGLDVDTRAYFTAATMIIAV
PTGIKIFSWIATMWGGSIQYKTPMLFAVGFIFLFTVGGLTGIVLANSGLDIALHDTYYVVAHFHYVLSMGAVFALFAGFYYWVGKIFGRTYPETLGQIHFWITFFGVNLTFFPMHFLGLSGMPRRIPDYPDAYAGWNALSSFGSYISVVGICCFFVVV
TITLSSGKNKRCAPSPWAVEENSTTLEWMVQSPPAFHTFGELPAIKETKSXXX
>nad1_119543_Anomodon_attenuatus
MRLYIIGILAKILGIIIPLLLGVAFLVLAERKIMASMQRRKGPNVVGLFGLLQPLADGLKLMIKEPILPSSANLFIFLMAPVMTFMLSLVAWAVIPFDYGMVLSDLNVGILYLFAISSLGVYGIITAGWSSNSKYAFLGALRSAAQMVSYEVSIGLIIITVLICVGSRNFSEIVIAQKQIWFAAPLFPVFIMFFISCLAETNRAPFDLPEAEAESVAGYNVEYSSMGFALFFLGEYANMILMSSLCTLLFLGGWLPILDIPIFYVIPGSIRFSIKVLFFLFVYIWVRAAFPRYRYDQLMRLGWKVFLPLSLAWVVFVSGVLVAFDWLP
```

The mitochondrial list was produced by the *DB_mitochondrion* mode (see 3. Pipeline description below part).


The nucleus list of genes, provided with OrthoSkim, contains:
+ 1133* single copy orthologous genes ([BUSCO](https://busco.ezlab.org)) dataset for Viridiplantae (v.10) retrieved from Waterhouse et al. (2017), with 10 ancestral amino acid sequences for each gene.
+ 352 ultra conserved element (UCE) designed for flowering plants and retrieved from Johnson et al. (2018).

\*Among the 1370 BUSCO of the initial dataset (which not mapped on chloroplasts and mitochondria), 237 genes mapped into the UCE dataset with diamond, and were removed to not included duplicates (here 1133 BUSCO in the file). UCE amino acid sequences were mapped into chloroplasts and mitochondria in order to keep nucleus compartment only.


```
head ~/OrthoSkim/ressources/refGenes.nu
```
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
```

**NOTE**: Personal nuclear genes list could be used in upcoming version of OrthoSkim. Users can easily change the nuclear sequence for other BUSCO at different level using the ancestral reconstructed sequences (see [BUSCO](https://busco.ezlab.org) website). For the chloroplast bank of genes, users has to put their own database in **<CHLORO_REF_CDS>**, **<CHLORO_REF_rRNA>** and **<CHLORO_REF_tRNA>** files following the restriction of header names.


### 3. Pipeline description
---------------------------

OrthoSkim uses different mode to align and extract targeted regions or to check assemblies. Regardless of genomic compartment chosen mode, a sample list files is required.
Genomic extraction can be performed directly for chloroplast and nucrdna target from annotations in `-m extraction_targeted` mode.
Otherwise, a non-directional assembly can be made with SPAdes on sample before running the pipeline in order to extract regions of interest using the `-m extraction_untargeted` mode.

**Note**: A *mode_done.log* file is created containing samples that were correctly processed, whereas unprocessed samples were added into *mode_error.log* file. This file could be used to remove processed samples from the initial **<LIST_FILES>** if the script has to be run. Command lines are also print if users want to rerun specific commands on samples.


#### 3.1 - Targeted mode: extraction from annotation

For each sample specified in **<LIST_FILES>**, OrthoSkim will perform genes extraction directly from annotation with `-m extraction_targeted` mode, according to a list of genes for `-t [chloroplast,nucrdna]` targets.

Results are output in **RES/** directory by creating subdirectories for each compartment and gene type, with a multifasta file per gene. For example, for chloroplastic CDS provided in **<CHLORO_GENES>**, OrthoSkim will output **RES/chloroplast_CDS/** subdirectory with CDS gene files.

```
ls -l ~/RES/chloroplast_CDS/
```
```
-rw-r--r--  1 pouchonc  staff   4874 16 avr 10:40 accD.fa
-rw-r--r--  1 pouchonc  staff   4952 16 avr 10:40 atpA.fa
-rw-r--r--  1 pouchonc  staff   4853 16 avr 10:40 atpB.fa
-rw-r--r--  1 pouchonc  staff   1580 16 avr 10:40 atpE.fa
-rw-r--r--  1 pouchonc  staff   2057 16 avr 10:40 atpF.fa
```


#### 3.2 - Untargeted mode: extraction from mapping to reference

The extraction of chloroplastic, nuclear and mitochondrial genes can  be driven on genomic assemblies performed with [SPAdes](http://cab.spbu.ru/software/spades/). This allowed to assembly all the genomic compartment and after to extract interested genes from these assemblies thanks to alignment into database.

##### 3.3.0 - Database (optional)

Users have to specify genes of reference into the *config_orthoskim.txt* file by following instructions for the sequence header.

OrthoSkim provides also a mode to create such database for the mitochondria in `-m DB_mitochondrion` mode. To do this, annotation files for mitochondria from genbank has to be download and put into the **<MITO_GENBANK>** single file.

For all these file, OrthoSkim will extract all notified CDS, rRNA and tRNA genes and align them into the CDS of *Arabidopsis_thaliana* as seeds thanks to *exonerate* (proteins *versus* proteins) to keep a standard gene name. Output **<MITO_REF_[CDS,rRNA,tRNA]>** files are created containing a bank of genes, all well identified.

**NOTE**: A new mode will be created soon to create a chloroplastic database from genbank annotations.

##### 3.3.1 - SPAdes assembly run

[SPAdes](http://cab.spbu.ru/software/spades/) assembly could be run in first using the `-m SPAdes_assembly` mode from the **<LIST_FILES>**, by accessing to the forward and reverse reads, and by keeping the sample name providing in the file. [SPAdes](http://cab.spbu.ru/software/spades/) will be run using **<THREADS>** and **<KMER>** specified in the *config_orthoskim.txt* file.

Othoskim will output then a **samplename/** subdirectory into the **PATHNAME_ASSEMBLY** given per sample included in the **<LIST_FILES>**.  

##### 3.3.2 - Preprocessing mode

After [SPAdes](http://cab.spbu.ru/software/spades/) runs, OrthoSkim has to preprocess SPAdes scaffolding contigs by renaming the file according to the same sample name provided in **<LIST_FILES>** and ordering them into **${RES}/${PATHNAME_ASSEMBLY}/Samples/** directory.
This is made under `-m SPAdes_reformate` mode. The extraction of genes will be made only on these renamed samples files.


##### 3.3.3 - Alignment extraction

The extraction of genomic regions of interested is made by choosing a target with `-t [nucleus,chloroplast_CDS, chloroplast_tRNA,chloroplast_rRNA,mitochondrion_CDS,mitochondrion_tRNA,mitochondrion_rRNA]` according to three steps:

1. **Reference selection**:

   For all chloroplastic and mitochondrial targets, OrthoSkim will first select the closest reference for each gene of our interested taxa from the given database of references.

   To achieve this, the selection is made according to the NCBI taxonomy (**<MODE_REF=taxonomy>**) according to taxid or by the phylogenetic distance (**<MODE_REF=distance**>) if a genus-level phylogenetic distance matrix is given in **<DISTANCE_MATRIX>**. In this mode, if the sample genus is not given into the matrix, the selection will automatically made on taxonomy. Finally, if the sample taxid does not exist in the NCBI taxonomy, OrthoSkim will use seeds as references.

   For the nucleus, no selection is made into the sequences as only BUSCO and UCE ancestral variants (already aligned) are used for the reference (a new mode using personal nuclear genes will come soon).   

   After this, if CDS are targeted, a [diamond](https://github.com/bbuchfink/diamond) database is created for each amino acid sequences provided in the retained sequences (with *diamond makedb*). Otherwise, a [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) database compilation (*makeblastdb* program) is formated.

2. **Alignment**:
   Mapping is made for each sample in two steps:

   **a.** Contigs are first aligned into the reference with [diamond](https://github.com/bbuchfink/diamond) *blastx* (for CDS targets) or with  with *blastn* command to quickly identify which of them align to the targeted genes (only hits with a minimal **<EVALUE>** were retained).

   **b.** Alignments are then conducted on these contigs from [exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate) incorporating all the appropriate gaps and frameshifts, and modelling of introns, by using the *protein2genome* mode for CDS target or the *genome2genome* mode for other targets. Only contigs with a minimal kmer coverage of **<COVERAGE>** and length of **<MINCONTLENGTH>** were considered for alignments into reference. A *gff* output table is created in **${RES}/Mapping/[nucleus,mitochondrion,chloroplast]/** folder for each sample.


3. **Genes extraction**:

   Extraction is conducted through python script from the gff  table. Type of gene structure (exon, intron or all:exon+intron) required during extraction from mapping are reported in **<NUC_TYPE>**,**<CHLORO_TYPE>** or **<MITO_TYPE>**. This step is conducted into multiple processors using the **<THREADS>** specified in the the *config_orthoskim.txt* file.

   Genes output files are created in the **${RES}/Extraction/[nucleus,mitochondrion,chloroplast]_[CDS,tRNA,rRNA]/** as following:

   ```
ls -l ~/RES/nucleus/
```

   ```
-rw-r--r--  1 pouchonc  staff  1758  5 jui 11:11 10104_BUSCO.fna
-rw-r--r--  1 pouchonc  staff  1964  5 jui 11:11 10521_BUSCO.fna
-rw-r--r--  1 pouchonc  staff  5071  5 jui 11:11 10785_BUSCO.fna
-rw-r--r--  1 pouchonc  staff  1400  5 jui 11:11 11487_BUSCO.fna
-rw-r--r--  1 pouchonc  staff  2040  5 jui 11:11 11505_BUSCO.fna
-rw-r--r--  1 pouchonc  staff  1778  5 jui 11:11 1504_BUSCO.fna
```

#### 3.4 - Summary statistics of assemblies

OrthoSkim allowed to output summary statistic on contigs assemblies for both *OrgAsm* (chloroplast+rdna) and [SPAdes](http://cab.spbu.ru/software/spades/) runs thanks to [QUAST](https://github.com/ablab/quast) by specifying the `-m stat_[chloro,rdna,SPAdes]` modes.

The output *transposed_report.txt* tab file will be in **${RES}/report_[chloro,rdnanuc,SPAdes]_assemblies/** directories given indication on the assembly. For example, for the chloroplast, we except the assembly of a single contig sizing 150,000 Nt in average as in the following example.  

```
head ~/RES/report_chloro_assemblies/transposed_report.txt
```
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

OrthoSkim allowed also to get statistic from gene extraction by using the `-m get_stat` mode for sequences (\*.fna) found in given `-p PATH`. The pipeline output a *report.tab* into this path containing:
+ gene name (gene_name)
+ samples found for this gene (#taxa)
+ the biggest length of sequence found (maxlen)
+ the number of sample covering this length at 100% (#1.0), 75% (#0.75), 50% (#0.5) and 25% (#0.25).

```
head ~/PATH/report.log
```

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

OrthoSkim is called step by step. Recommendations about steps are given in the previous description (section 3). After edition of the *tools.sh* and *config_orthoskim.txt* files (with all reuired files and formats), orthoskim is called by using the differents modes.

To run *extraction_targeted* mode if annotations files are given:

```
./orthoskim -m extraction_targeted -t chloroplast -c config_orthoskim.txt
```
   ```
./orthoskim -m extraction_targeted -t nucrdna -c config_orthoskim.txt
```

For the *extraction_untargeted* mode, the script should be called as following with first:

```
./orthoskim -m SPAdes_assembly -c config_orthoskim.txt
./orthoskim -m SPAdes_reformate -c config_orthoskim.txt
```

and then, for example if you want to extract all chloroplastic, mitochondrial and nuclear genes, users run OrthoSkim as following:

```
./orthoskim -m extraction_untargeted -t chloroplast_CDS -c config_orthoskim.txt
./orthoskim -m extraction_untargeted -t chloroplast_rRNA -c config_orthoskim.txt
./orthoskim -m extraction_untargeted -t chloroplast_tRNA -c config_orthoskim.txt
./orthoskim -m extraction_untargeted -t nucleus -c config_orthoskim.txt
./orthoskim -m extraction_untargeted -t mitochondrion_CDS -c config_orthoskim.txt
./orthoskim -m extraction_untargeted -t mitochondrion_rRNA -c config_orthoskim.txt
./orthoskim -m extraction_untargeted -t mitochondrion_tRNA -c config_orthoskim.txt
```
**Note**: *SPAdes_assembly* and *SPAdes_reformate* have to be run just once before the *extraction_untargeted* mode regardless of the targets.

For summary statistic, OrthoSkim could be run using the  `-m [stat_chloro,stat_rdna,stat_SPAdes]`, or using the `-m get_stat` with `-p PATH` option as following:

```
./OrtoSkim_v0.0.sh -m stat_SPAdes -c config_orthoskim.txt
./OrtoSkim_v0.0.sh -m get_stat -c config_orthoskim.txt -p path_to_extracted_files
```

### 5. Additional modes for PhyloDB users
#### 5.1 - Indexing files

This step allowed to guide the pipeline according to an indexing tab based on samples location into the [GriCAD](https://gricad-doc.univ-grenoble-alpes.fr/hpc/description/) infrastructures on the **/bettik/LECA/phyloskims/release/** folder. This tab is produced by the `-m indexing` mode. This allowed to screen each sample that will be used for extraction from the **<PATH_FIND_CHLORO>** directory. Unwanted samples must be removed from the list before processing other modes.


```
./orthoskim -m indexing -c config_orthoskim.txt
```

#### 5.2 - phyloDB chloroplast references

OrthoSkim provides  a mode to create a chloroplastic database from the annotated chloroplasts found with the *indexing* mode, using the `-m DB_chloroplast` mode. To do this, all genes found in these annotations files are extracted with the header restrictions. Output files are created in  

### 6. References
--------------------
+ Bankevich, A., Nurk, S., Antipov, D., Gurevich, A.A., Dvorkin, M., Kulikov, A.S., et al. (2012). SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. J. Comput. Biol., 19, 455–477.
+ Buchfink, B., Xie, C. & Huson, D.H. (2015). Fast and sensitive protein alignment using DIAMOND. Nature Methods, 12, 59–60.
+ Gurevich, A., Saveliev, V., Vyahhi, N. & Tesler, G. (2013). QUAST: quality assessment tool for genome assemblies. Bioinformatics, 29, 1072–1075.
+ Johnson, M.G., Pokorny, L., Dodsworth, S., Botigué, L.R., Cowan, R.S., Devault, A., et al. (n.d.). A Universal Probe Set for Targeted Sequencing of 353 Nuclear Genes from Any Flowering Plant Designed Using k-Medoids Clustering. Syst Biol.
+ Simão, F.A., Waterhouse, R.M., Ioannidis, P., Kriventseva, E.V. & Zdobnov, E.M. (2015). BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. Bioinformatics, 31, 3210–3212.
+ Slater, G.S.C. & Birney, E. (2005). Automated generation of heuristics for biological sequence comparison. BMC Bioinformatics, 6, 31.
+ Waterhouse, R.M., Seppey, M., Simão, F.A., Manni, M., Ioannidis, P., Klioutchnikov, G., et al. (2018). BUSCO Applications from Quality Assessments to Gene Prediction and Phylogenomics. Molecular Biology and Evolution, 35, 543–548.
