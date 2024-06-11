############################################### CUSTOM REPEAT LIBRARY CONSTRUCTION #######################
## Repeat library construction using RepeatModeler + LTR_retriever
# RM
BuildDatabase -name Neltuma_pallida -engine ncbi ppa_v1.asm.fasta
RepeatModeler -threads 40 -engine ncbi -database Neltuma_pallida 2>&1 | tee Repeatmodeler.log
# LTR
gt suffixerator -db ppa_v1.asm.fasta -indexname ppa_v1.asm.fasta -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index ppa_v1.asm.fasta -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > ppa_v1.asm.fasta.harvest.scn
perl ~/Bioprograms/LTR_FINDER_parallel/LTR_FINDER_parallel/LTR_FINDER_parallel -seq ppa_v1.asm.fasta -threads 30 -harvest_out -size 1000000 -time 300
cat ppa_v1.asm.fasta.harvest.scn ppa_v1.asm.fasta.finder.combine.scn > ppa_v1.asm.fasta.rawLTR.scn
LTR_retriever -genome ppa_v1.asm.fasta -inharvest ppa_v1.asm.fasta.rawLTR.scn -threads 30
# Joining outputs from RM and LTR
cat consensi.fa.classified ppa_v1.asm.fasta.LTRlib.fa > RepeatsTotal.fa
cat RepeatsTotal.fa | seqkit fx2tab | awk '{print "NelPal1_"$0}' | seqkit tab2fx > RepeatsTotal.prefix.fa
cat RepeatsTotal.prefix.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > RepeatsKnown.prefix.fa #1444 known repeats
cat RepeatsTotal.prefix.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > RepeatsUnknown.prefix.fa #1363 unknown repeats
## Mask the genome step by step, based on Daren Card suggestions on his blog
# First, soft masking simple repeats
RepeatMasker -pa 30 -a -e ncbi -dir 01_simple_out -noint -xsmall ppa_v1.asm.fasta 2>&1 | tee logs/01_simplemask.log
# Simple repeats 167011 elements, 7417610 bp, 1.84%
# Low complexity 40628 elements, 2100656 bp, 0.52%
# Second, mask using Repbase Viridiplantae
RepeatMasker -pa 40 -a -e ncbi -dir 02_viridiplantae_out -nolow -species viridiplantae 01_simple_out/ppa_v1.asm.simple_mask.masked.fasta 2>&1 | tee logs/02_viridiplantaemask.log
# Retroelements 53454 elements, 28879797 bp, 7.15%
# DNA transposons 31579 elements, 7102580 bp, 1.76%
# Rolling-circles 5017, 1008258 bp, 0.25%
# Unclassified 247, 35860 bp, 0.01%
# Small RNA, 745 elements, 116254 bp, 0.03%
# Satellites, 414 elements, 128591 bp, 0.03%
# Third, mask using the known repetitive elements
RepeatMasker -pa 40 -a -e ncbi -dir 03_known_out -nolow -lib RepeatsKnown.prefix.fa 02_viridiplantae_out/ppa_v1.asm.viridiplantae_mask.masked.fasta 2>&1 | tee logs/03_knownmask.log
# Retroelements, 209649 elements, 106138982 bp, 26.28%
# DNA transposons, 39789 elements, 16208558 bp, 4.01%
# Rolling-circles, 4286 elements, 1868879 bp, 0.46%
# Small RNA, 2389 elements, 273413 bp, 0.07%
# Simple repeats, 5 elements, 179 bp, 0.00%
# Fourth, mask using the unknown repetitive elements
RepeatMasker -pa 40 -a -e ncbi -dir 04_unknown_out -nolow -lib RepeatsUnknown.prefix.fa 03_known_out/ppa_v1.asm.known_mask.masked.fasta 2>&1 | tee logs/04_unknownmask.log
cat 01_simple_out/ppa_v1.asm.simple_mask.cat.gz 02_tetrapoda_out/ppa_v1.asm.viridiplantae_mask.cat.gz 03_known_out/ppa_v1.asm.known_mask.cat.gz 04_unknown_out/ppa_v1.asm.unknown_mask.cat.gz > 05_full_out/ppa_v1.asm.full_mask.cat.gz
cat 01_simple_out/ppa_v1.asm.simple_mask.out  <(cat 02_tetrapoda_out/ppa_v1.asm.viridiplantae_mask.out | tail -n +4) <(cat 03_known_out/ppa_v1.asm.known_mask.out | tail -n +4) <(cat 04_unknown_out/ppa_v1.asm.unknown_mask.out | tail -n +4) > 05_full_out/ppa_v1.asm.full_mask.out
cat 01_simple_out/ppa_v1.asm.simple_mask.out > 05_full_out/ppa_v1.asm.simple_mask.out
# Summarize the masking results
ProcessRepeats -a -species tetrapoda 05_full_out/ppa_v1.asm.full_mask.cat.gz 2>&1 | tee logs/05_fullmask.log
# Get repeats GFF3
./rmOutToGFF3custom -o 05_full_out/ppa_v1.asm.full_mask.out > 05_full_out/ppa_v1.asm.full_mask.gff3
./rmOutToGFF3custom -o 05_full_out/ppa_v1.asm.simple_mask.out > 05_full_out/ppa_v1.asm.simple_mask.gff3
./rmOutToGFF3custom -o 05_full_out/ppa_v1.asm.complex_mask.out > 05_full_out/ppa_v1.asm.complex_mask.gff3
# Create masked genome FASTA files for annotation
bedtools maskfasta -soft -fi ppa_v1.asm.fasta -bed 05_full_out/ppa_v1.asm.simple_mask.gff3 -fo 05_full_out/ppa_v1.asm.fasta.simple_mask.soft.fasta
bedtools maskfasta -fi 05_full_out/ppa_v1.asm.simple_mask.soft.fasta -bed 05_full_out/ppa_v1.asm.complex_mask.gff3 -fo 05_full_out/ppa_v1.asm.simple_mask.soft.complex_mask.hard.fasta

############################################### TRANSCRIPTS EVIDENCE: TRANSCRIPTOME #######################
guppy_basecaller --recursive --input_path /media/prosopis/D/2022_Prosopis/arn_recovery_2024/RNA_Prosopis --save_path Lib1 --config dna_r9.4.1_450bps_sup.cfg --records_per_fastq 0 --trim_strategy dna --disable_trim_barcodes --barcode_kits 'SQK-PCB109' --compress_fastq --device 'auto'
guppy_basecaller --recursive --input_path /media/prosopis/D/2022_Prosopis/arn_recovery_2024/LibTogether --save_path Lib2 --config dna_r9.4.1_450bps_sup.cfg --records_per_fastq 0 --trim_strategy dna --disable_trim_barcodes --barcode_kits 'SQK-PCB109' --compress_fastq --device 'auto'
cat */*/pass/*.fastq > cDNA_Alltissues.fastq
## Read preprocessing
~/Bioprograms/porechop/Porechop/porechop-runner.py -i cDNA_Alltissues.fastq -o cDNA_Alltissues_porechop.fastq -t 30
NanoFilt -q 9 -l 500 cDNA_Alltissues_porechop.fastq > cDNA_Alltissues_porechop_nanofilt_q9_l500.fastq
## Mapping reads to genome
minimap2 -N 1 -a -x splice -g2000 -G5k -t 35 ppa_v1.asm.fasta cDNA_Alltissues_porechop_nanofilt_q9_l500.fastq > cDNA_alignment_ref.sam
# N -> Output secondary alignments. When using long read cDNA it can output unspliced pseudogenes.
#-x splic means orientation of transcript is unknown and does two rounds of alignment to infer the orientation
# Do not take into account extremely long introns
# samtools flagstat -@ 30 cDNA_alignment_ref.sam
# 3952198 + 0 in total (QC-passed reads + QC-failed reads)
# 3547405 + 0 primary
# 257119 + 0 secondary
# 147674 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 3860426 + 0 mapped (97.68% : N/A)
# 3455633 + 0 primary mapped (97.41% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)
samtools sort -@ 30 -o cDNA_alignment_ref.sorted.bam cDNA_alignment_ref.sam
## Construct reference-guided transcriptome using StringTie
stringtie cDNA_alignment_ref.sorted.bam -m 500 -t -c 1.5 -f 0.05 -s 4.5 -g 250 -p 35 -o transcriptome.gtf -l NPTR
# Conservative parameters for assembly of transcripts based on reference genome. Also produces isoforms. Minimum 500 length, 4.5 per-base coverage for single exon transcript
# Results in exons, could be CDS or UTRs
# MAKER allows to include transcriptomes from close species. In this case, leaf transcriptomes from N. alba and P. cineraria were downloaded from NCBI.
############################################### PROTEIN EVIDENCE: PLANT SWISSPROT #######################
# For protein evidence, protein sequences from the database UniProt – SwissProt restricted to plants were downloaded in .dat format. To parse the sequences, the perl script from SwissKnife varsplic.pl was used.
perl varsplic.pl -input unprot_sprot_plants.dat -check_vsps -crosscheck -error varsplic.err -fasta Var_USP_plants.fasta -which full #Also do without isoforms for functional annotation
############################################### ALTERNATIVE TRANSCRIPTS #######################
# Download transcriptome for N. alba and P. cineraria
############################################### MAKER ANNOTATION: ROUND 1 #######################
# To use MAKER, first we need to create and edit control files. Using maker -CTL will create four files. The file maker_exe.ctl needs to have the correct paths for the programs used. The file maker_opts.ctl has all the configuration for the run.
maker -CTL
# For the first round, modify the following lines:
# genome= 	Add the path to the masked assembly fasta.
# altest=	Add the path to each close species transcriptome. Similar to est, separate with a comma and add the name of the species of origin.
# protein=	Add the path to the swissprot sequences.
# est2genome=	Turn this from 0 to 1
# protein2genome=	Turn this from 0 to 1
# min_protein=	Change from 0 to 50 to require a minimum protein length
# To start MAKER, use the command “maker” in the same folder where the control files are created. To run MAKER more efficiently, run three different tasks in the same directory. MAKER knows when a contig is being processed and will work on another one. Do not mask
maker -RM_off 2> maker1.err &
maker -RM_off 2> maker2.err &
maker -RM_off 2> maker3.err &
# To monitor how each run is going, you can print the log continuously using the following command:
tail -f maker1.err
## Train gene models for second annotation round
# After MAKER is finished, go to the output folder and extract the gff and sequences. These will be used to get the gene models for training with SNAP and Augustus.
gff3_merge -n -s -d ppa_v1.asm.fa_master_datastore_index.log > Algarrobo_rnd1.gff
fasta_merge -d ppa_v1.asm.fa_master_datastore_index.log
# To train for SNAP, gene models are first filtered to discard models with an AED higher than 0.25 and coding proteins shorter than 50 amino acids.
maker2zff -x 0.25 -l 50 -d /path/to/ppa_v1.asm.fa_master_datastore_index.log
# Then, the training sequences and annotations are gathered with additional 1 kb upstream and downstream.
fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
# Next, the training parameters are created, and the HMM file is assembled. This will be used for MAKER to predict gene models using the trained SNAP.
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
hmm-assembler.pl genome params > Algarrobo_rnd1.hmm
# To train for Augustus, a BUSCO feature is used, so in this case the 2023-01-busco environment needs to be activated. Make sure to index the assembly fasta before starting.
samtools faidx assembly.fasta
# The first step is to gather the regions that contain the mRNA annotations with up to 1 kb on each side.
awk -v OFS="\t" '{if ($3 == "mRNA") print $1, $4, $5 }' Algarrobo_rnd1.gff | while read rna; do scaffold=`echo ${rna} | awk '{ print $1 }'`; end=`cat ppa_v1.asm.fa.fai | awk -v scaffold="${scaffold}" -v OFS="\t" '{ if ($1 == scaffold) print $2 }'`; echo ${rna} | awk -v end="${end}" -v OFS="\t" '{ if ($2 < 1000 && (end - $3) < 1000) print $1, "0", end; else if ((end - $3) < 1000) print $1, "0", end; else if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }'; done | bedtools getfasta -fi ppa_v1.asm.fa -bed - -fo Algarrobo_rnd1.transcripts1000.fasta
# Then, BUSCO will do a reannotation of the extracted regions using BLAST and built-in HMMs for a set of conserved genes. Then, with the “--long" option, BUSCO will optimize the HMM search model to train Augustus and produce a trained HMM for MAKER. In this case, the set of conserved genes is from embryophyta_odb10, and the starting HMM model is from arabidopsis.
busco -i Algarrobo_rnd1.transcripts1000.fasta -o Algarrobo_rnd1_maker -l embryophyta_odb10 -m genome -c 30 --long --augustus_species arabidopsis --tar --augustus_parameters=--progress=true
# Once BUSCO is complete, we need to modify the names and the content of some files before copying it to the config path of Augustus.
mkdir /opt/augustus/config/species/BUSCO_Algarrobo_rnd1_maker
cp BUSCO_Algarrobo_rnd1_maker* /opt/augustus/config/species/BUSCO_Algarrobo_rnd1_maker/.
# ---------------------------------------------------
# |Results from dataset embryophyta_odb10            |
# ---------------------------------------------------
# |C:79.2%[S:75.4%,D:3.8%],F:3.9%,M:16.9%,n:1614     |
# |1278    Complete BUSCOs (C)                       |
# |1217    Complete and single-copy BUSCOs (S)       |
# |61    Complete and duplicated BUSCOs (D)          |
# |63    Fragmented BUSCOs (F)                       |
# |273    Missing BUSCOs (M)                         |
# |1614    Total BUSCO groups searched               |
# ---------------------------------------------------

############################################### MAKER ANNOTATION: ROUND 2 #######################
# First, in the directory ../1-Output, the gff was separated to obtain the est2genome, protein2genome and repeat annotations. These will be used for the following MAKER runs.
cat Algarrobo_rnd1.gff | grep -P “\test2genome” > Algarrobo_rnd1.est2genome.gff
cat Algarrobo_rnd1.gff | grep -P “\tcdna2genome” > Algarrobo_rnd1.cdna2genome.gff
cat Algarrobo_rnd1.gff | grep -P “\tprotein2genome” > Algarrobo_rnd1.protein2genome.gff
# Copy all control files from the first round and modify the following options:
# est_gff=	Add the path to the est2genome gff from the first round
# altest_gff=	Add the path to the cdna2genome gff from the first round
# protein_gff=	Add the path to the protein2genome gff from the first round
# snaphmm=	Add the path to the SNAP HMM file
# augustus_species=	Write “BUSCO_Algarrobo_rnd1_maker”. As the folder is in the config file, it should be recognized by Augustus.
# ---------------------------------------------------
# |Results from dataset embryophyta_odb10            |
# ---------------------------------------------------
# |C:86.5%[S:83.8%,D:2.7%],F:4.8%,M:8.7%,n:1614      |
# |1396    Complete BUSCOs (C)                       |
# |1353    Complete and single-copy BUSCOs (S)       |
# |43    Complete and duplicated BUSCOs (D)          |
# |78    Fragmented BUSCOs (F)                       |
# |140    Missing BUSCOs (M)                         |
# |1614    Total BUSCO groups searched               |
# ---------------------------------------------------

############################################### MAKER ANNOTATION: ROUND 3 #######################
# For the third annotation run, repeat the training with the new gene models. For Augustus, change the initial species to Neltuma_pallida, and the output as BUSCO_Algarrobo_rnd2_maker.
# Modify the following options in the control file to use the retrained gene models
# snaphmm=	Add the path to the new SNAP HMM file
# augustus_species=	Write “BUSCO_Algarrobo_rnd2_maker”
# Activate alt_splice=1
# ---------------------------------------------------
# |Results from dataset embryophyta_odb10            |
# ---------------------------------------------------
# |C:89.9%[S:78.4%,D:11.5%],F:3.0%,M:7.1%,n:1614     |
# |1451    Complete BUSCOs (C)                       |
# |1266    Complete and single-copy BUSCOs (S)       |
# |185    Complete and duplicated BUSCOs (D)         |
# |49    Fragmented BUSCOs (F)                       |
# |114    Missing BUSCOs (M)                         |
# |1614    Total BUSCO groups searched               |
# ---------------------------------------------------
#Third round: 20,069 genes, 22343 gene models

## Remove transcripts containing masked regions from the evidence-based gene models
seqkit grep -v --by-seq -p "X" ../../2024-05-29_Npal_rnd3.all.maker.proteins.fasta > 2024-06-02_Npal_rnd3.nomasked.proteins.fasta
seqkit fx2tab 2024-06-02_Npal_rnd3.nomasked.proteins.fasta | cut -f1 | cut -d' ' -f1 | seqkit grep -r -n -f - ../../2024-05-29_Npal_rnd3.all.maker.transcripts.fasta > 2024-06-02_Npal_rnd3.nomasked.transcripts.fasta
# Filter the original gff
seqkit grep --by-seq -p "X" ../../2024-05-29_Npal_rnd3.all.maker.proteins.fasta | seqkit fx2tab | cut -f1 | cut -d' ' -f1 > remove_names.txt
cat remove_names.txt | grep -v -f - ../../2024-05-29_Npal_rnd3.all.maker.gff > 2024-06-02_Npal_rnd3.filtered1.gff

cat remove_names.txt | cut -d'-' -f1-5 | grep -f - 2024-06-02_Npal_rnd3.filtered1.gff | grep -P "\tmRNA\t" | cut -f9 | sed 's/Parent=/\t/' | cut -f2 | cut -d';' -f1 | sort | uniq | grep -v -f - remove_names.txt | cut -d'-' -f1-5 | sed 's/$/;/' | grep -v -f - 2024-06-02_Npal_rnd3.filtered1.gff > 2024-06-02_Npal_rnd3.filtered2.gff
# From 20,069 to 18,787 genes
# From 22,343 to 20,990 gene models

### Filtered busco (only maker, no only predicted)
# ---------------------------------------------------
# ---------------------------------------------------
# |Results from dataset embryophyta_odb10            |
# ---------------------------------------------------
# |C:85.7%[S:74.3%,D:11.4%],F:3.1%,M:11.2%,n:1614    |
# |1383    Complete BUSCOs (C)                       |
# |1199    Complete and single-copy BUSCOs (S)       |
# |184    Complete and duplicated BUSCOs (D)         |
# |50    Fragmented BUSCOs (F)                       |
# |181    Missing BUSCOs (M)                         |
# |1614    Total BUSCO groups searched               |
# ---------------------------------------------------

## Get only predicted gene models without masked regions
cat ../../2024-05-29_Npal_rnd3.all.maker.non_overlapping_ab_initio.proteins.fasta | seqkit grep -v --by-seq -p "X" - > 2024-06-02_Npal_rnd3.nomasked.non_overlapping.proteins.fasta
# From 20,536 to 17,237 non-overlapping predicted proteins

## Get unmasked non overlapping proteins
/data/lastexpansion/Renato/Bioprograms/Interproscan/interproscan-5.67-99.0/interproscan.sh -appl PfamA -iprlookup -goterms -f tsv -i ../2024-06-02_Npal_rnd3.nomasked.non_overlapping.proteins.fasta
# 3682 proteins with Pfam domains
cat 2024-06-02_Npal_rnd3.nomasked.non_overlapping.proteins.fasta.tsv | cut -f1 | sort | uniq | sed 's/processed/abinit/' | grep -f - ../../../2024-05-29_Npal_rnd3.all.maker.gff > ../2024-06-02_Npal_rnd3.nomasked.predicted.gff

### Another maker round to format the predicted gff
####################################################
gff3_merge -n -s -d ppa_v1.asm_master_datastore_index.log > 2024-05-31_Npal_rnd3.predicted.gff
fasta_merge -o 2024-05-31_Npal_rnd3 -d ppa_v1.asm_master_datastore_index.log

gff3_merge -o 2024-06-02_Npal_rnd3.merged.gff 2024-06-02_Npal_rnd3.filtered2.gff ppa_v1.asm.maker.output/2024-06-02_Npal_rnd3.predicted.gff
cat 2024-06-02_Npal_rnd3.nomasked.transcripts.fasta ppa_v1.asm.maker.output/2024-06-02_Npal_rnd3.all.maker.transcripts.fasta > 2024-06-02_Npal_rnd3.merged.transcripts.fasta
cat 2024-06-02_Npal_rnd3.nomasked.proteins.fasta ppa_v1.asm.maker.output/2024-06-02_Npal_rnd3.all.maker.proteins.fasta > 2024-06-02_Npal_rnd3.merged.proteins.fasta
# Standard gene set: 22,469 genes and 24,672 gene models
# ---------------------------------------------------
# |Results from dataset embryophyta_odb10            |
# ---------------------------------------------------
# |C:89.1%[S:77.6%,D:11.5%],F:3.2%,M:7.7%,n:1614     |
# |1437    Complete BUSCOs (C)                       |
# |1252    Complete and single-copy BUSCOs (S)       |
# |185    Complete and duplicated BUSCOs (D)         |
# |52    Fragmented BUSCOs (F)                       |
# |125    Missing BUSCOs (M)                         |
# |1614    Total BUSCO groups searched               |
# ---------------------------------------------------

### Rename annotations
gff3sort.pl --precise --chr_order natural 2024-06-02_Npal_rnd3.merged.gff > 2024-06-02_Npal_rnd3.merged.sorted.gff
maker_map_ids --prefix NPAL_ --suffix _R --iterate 1 --justify 6 2024-06-02_Npal_rnd3.merged.sorted.gff > 2024-06-02_Npal_rnd3.merged.map
map_gff_ids 2024-06-02_Npal_rnd3.merged.map 2024-06-02_Npal_rnd3.merged.sorted.gff
map_fasta_ids 2024-06-02_Npal_rnd3.merged.map 2024-06-02_Npal_rnd3.merged.proteins.fasta
map_fasta_ids 2024-06-02_Npal_rnd3.merged.map 2024-06-02_Npal_rnd3.merged.transcripts.fasta

## Functional Annotation
#BLAST
makeblastdb -in 2024-03-27-uniprot_sprot_plants.fasta -input_type fasta -dbtype prot
blastp -db 2024-03-27-uniprot_sprot_plants.fasta -query ../2024-05-30_Npal_rnd3.merged.proteins.fasta -out maker2uni.blastp -evalue .000001 -outfmt 6 -num_alignments 1 -seg yes -soft_masking true -lcase_masking -max_hsps 1 -num_threads 25
maker_functional_fasta 2024-03-27-uniprot_sprot_plants.fasta maker2uni.blastp ../2024-06-02_Npal_rnd3.merged.transcripts.fasta > ../2024-06-02_Npal_rnd3.merged.transcripts.blast.fasta
cat 2024-06-02_Npal_rnd3.merged.transcripts.blast.fasta | seqkit fx2tab | sort -k1,1 | seqkit tab2fx > 2024-06-02_Npal_rnd3.merged.transcripts.sorted.blast.fasta
maker_functional_fasta 2024-03-27-uniprot_sprot_plants.fasta maker2uni.blastp ../2024-06-02_Npal_rnd3.merged.proteins.fasta > ../2024-06-02_Npal_rnd3.merged.proteins.blast.fasta
cat 2024-06-02_Npal_rnd3.merged.proteins.blast.fasta | seqkit fx2tab | sort -k1,1 | seqkit tab2fx > 2024-06-02_Npal_rnd3.merged.proteins.sorted.blast.fasta
maker_functional_gff 2024-03-27-uniprot_sprot_plants.fasta maker2uni.blastp ../2024-06-02_Npal_rnd3.merged.sorted.gff > ../2024-06-02_Npal_rnd3.merged.sorted.blast.gff

#Interpro
/data/lastexpansion/Renato/Bioprograms/Interproscan/interproscan-5.67-99.0/interproscan.sh -appl PfamA -iprlookup -goterms -f tsv -i ../2024-06-02_Npal_rnd3.merged.proteins.fasta
ipr_update_gff ../2024-05-31_Npal_rnd3.merged2.sorted.gff 2024-06-02_Npal_rnd3.nomasked.non_overlapping.proteins.tsv 2024-05-31_Npal_rnd3.merged2.sorted.interpro.gff

#Emapper
emapper.py -i ../1-Blast/2024-05-31_Npal_rnd3.merged2.proteins.sorted.blast.fasta --output_dir EM_Viridiplantae --cpu 30 --dbmem --report_orthologs --excel --decorate_gff ../1-Blast/2024-05-31_Npal_rnd3.merged2.sorted.blast.gff -o 2024-05-31_Npal_rnd3.merged2.sorted.blast --tax_scope Viridiplantae --data_dir ../../Old_Functional/2-Emapper/EM_databases/

############################################### ANNOTATION EVALUATION #######################
# For visualizaiton, separate the GFF in parts
cat Algarrobo_rnd1.gff | grep -P “\test2genome” > EST_Tissue.gff # For all tissues
cat Algarrobo_rnd1.gff | grep -P “\cdna2genome” | grep “ALTEST_Species” > ALTEST_Species.gff # For different species
cat Algarrobo_rnd1.gff | grep -P “\tprotein2genome”> PROT_Swissprot.gff
# Can use repeat gffs divided by repeat type
# The genome and all these tracks will be uploaded to IGV for visualization.
igv.sh
# It may be needed to sort and index the genome and some of the tracks. If IGV closes unexpectedly, change the configuration to allocate more memory.
gff3sort.pl --precise --chr_order natural file.gtf/gff | bgzip > file.gtf/gff.gz
tabix -p gff file.gtf/gff.gz
# A good metric for the annotation quality is the BUSCO completeness of the transcriptome. To do this, in the Busco environment use the following command. A BUSCO completeness over 90% is good.
busco -i Algarrobo_rndN.all.maker.transcripts.fasta -o Algarrobo_rndN_maker -l embryophyta_odb10 -m transcriptome -c 30 --tar

############################################### ANNOTATION STATISTICS #######################
QIfile=Npal_QImetrics.txt
# Number of gene models without UTR: 11,863
cat $QIfile | awk '$2==0 && $9==0' | wc -l
# Number of gene models with both UTRs: 11,500
cat $QIfile | awk '$2>0 && $9>0' | wc -l
# Number of gene models with 5UTR: 12,833
cat $QIfile | awk '$2>0' | wc -l
# Number of gene models with 3UTR: 13,802
cat $QIfile | awk '$9>0' | wc -l
# Mean 5UTR length: 315.013
cat $QIfile | cut -f2 | awk '$1>0' | awk '{s+=$0}END{print s/NR}'
# Mean 3UTR length: 500.779
cat $QIfile | cut -f9 | awk '$1>0' | awk '{s+=$0}END{print s/NR}'
# Number gene models with exons covered by EST: 15,306
cat $QIfile | awk '$4>0' | wc -l
# Number of gene models with exons covered by EST, or EST or Protein: 20,601
cat $QIfile | awk '$5>0' | wc -l
# Average fraction of exons covered by EST: 0.854569
cat $QIfile | awk '$4>0' | cut -f4 | awk '$1>0' | awk '{s+=$0}END{print s/NR}'
# Average fraction of exons covered by EST or Protein: 0.932506
cat $QIfile | awk '$5>0' | cut -f5 | awk '$1>0' | awk '{s+=$0}END{print s/NR}'
# Number of gene models with more than 1 exon: 23,838
cat $QIfile | awk '$8>1' | wc -l
# Number of gene models with more than 1 exon with splice sites covered by EST: 14,761
cat $QIfile | awk '$8>1' | awk '$3>0' | wc -l
# Average fraction of splice sites covered by EST: 0.874982
cat $QIfile | awk '$8>1' | awk '$3>0' | cut -f3 | awk '{s+=$0}END{print s/NR}'
## Functional annotation statistics
# Number of gene models annotated with InterPro: 23,201
cat 2024-05-31_Npal_rnd3.merged2.proteins.sorted.blast.fasta.tsv | cut -f1 | sort | uniq | wc -l
GFFfile=2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff
# Number of gene models annotated with SwissProt
cat $GFFfile | grep -P "\tmRNA" | cut -f9 | grep -v "Note=Protein of unknown function" | wc -l
# Number of gene models annotated with EggNOG
cat $GFFfile | grep -P "\tmRNA" | cut -f9 | grep "em_target" | wc -l
# Functional annotation overlaps were drawn in a Venn diagram
# Number of gene models with COG categories: 22,440
cat $GFFfile | grep -P "\tmRNA" | grep "em_COG" | sed 's/em_COG_cat=/#/' | cut -d'#' -f2 | cut -d';' -f1 | grep -v "None" | wc -l
# Get list of COG categories for all gene models:
cat $GFFfile | grep -P "\tmRNA" | grep "em_COG" | sed 's/em_COG_cat=/#/' | cut -d'#' -f2 | cut -d';' -f1 | grep -v "None" | grep -o . | sort | uniq -c | sort -rn | head
# Get matches to SwissProt proteins
cat $GFFfile | grep -P "\tmRNA\t" | grep "Note=Similar to" | cut -f9 | sed 's/Note=/#/' | cut -d'#' -f2 | cut -d';' -f1 | sed 's/(Fragment/(fragment/' | sed 's/([A-Z][a-z][a-z]/\t/' | cut -f2 | cut -d')' -f1 | sort | uniq -c | sort -rn | head
# Get KEGG pathway list for all gene models
cat $GFFfile | grep -P "\tmRNA" | grep "em_KEGG_Pathway" | sed 's/em_KEGG_Pathway=/\t/' | cut -f10 | cut -d';' -f1 | sed 's/n/\n/g' | sort | uniq -c | sort -rn | head
# Get PFAM list for all gene models
cat 2024-05-31_Npal_rnd3.merged2.proteins.sorted.blast.fasta.tsv | cut -f1,5 | sort | uniq | cut -f2 | sort | uniq -c | sort -rn | head
# Get number of protein kinase domains identified with EggNOG: 814
cat $GFFfile | grep -P "\tmRNA" | grep "em_PFAMs=" | sed 's/em_PFAMs=/\t/' | cut -f10 | cut -d';' -f1 | grep -w "Pkinase" | wc -l
# Get number of protein kinase tyr domains identified with EggNOG:568
cat $GFFfile | grep -P "\tmRNA" | grep "em_PFAMs=" | sed 's/em_PFAMs=/\t/' | cut -f10 | cut -d';' -f1 | grep -w "Pkinase_Tyr" | wc -l
# Get number of LRR domains identified with EggNOG: 463
cat $GFFfile | grep -P "\tmRNA" | grep "em_PFAMs=" | sed 's/em_PFAMs=/\t/' | cut -f10 | cut -d';' -f1 | grep "LRR" | wc -l
# Get number of PPR domains identified with EggNOG: 465
cat $GFFfile | grep -P "\tmRNA" | grep "em_PFAMs=" | sed 's/em_PFAMs=/\t/' | cut -f10 | cut -d';' -f1 | grep "PPR" | wc -l
# Number of gene models involved in carbon fixation modules: 88
cat $GFFfile | grep -P "\tmRNA\t" | grep "map00710" | wc -l
# Number of genes involved in carbon fixation modules: 80
cat $GFFfile | grep -P "\tmRNA\t" | grep "map00710" | cut -f9 | cut -d'_' -f1,2 | sort | uniq | wc -l
# Number of gene models involved in biosynthesis of secondary metabolites: 1,775
cat $GFFfile | grep -P "\tmRNA" | grep "map01110" | wc -l
# Number of genes involved in biosynthesis of secondary metabolites: 1,658
cat $GFFfile | grep -P "\tmRNA\t" | grep "map01110" | cut -f9 | cut -d'_' -f1,2 | sort | uniq | wc -l
# Number of gene models also involved in biosynthesis of secondary metabolites: 15
cat $GFFfile | grep -P "\tmRNA\t" | grep "map00940" | grep -v "map01110" | cut -f9 | cut -d'_' -f1,2 | sort | uniq | wc -l
# Number of genes involved in signal transduction pathways: 312
cat $GFFfile | grep -P "\tmRNA\t" | grep "map04075" | cut -f9 | cut -d'_' -f1,2 | sort | uniq | wc -l
# Number of genes involved in plant-pathogen interaction: 245
cat $GFFfile | grep -P "\tmRNA\t" | grep "map04626" | cut -f9 | cut -d'_' -f1,2 | sort | uniq | wc -l
# Number of genes involved in MAPK signaling pathways: 188
cat $GFFfile | grep -P "\tmRNA\t" | grep "map04016" | cut -f9 | cut -d'_' -f1,2 | sort | uniq | wc -l
