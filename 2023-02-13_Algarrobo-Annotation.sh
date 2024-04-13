############################################### CUSTOM REPEAT LIBRARY CONSTRUCTION #######################
## Repeat library construction
# The process described here is the same as in the Advanced Library Protocol in the MAKER wiki. First, MITEs were identified using MITEHunter.
perl MITE_Hunter_manager.pl -i ppa_v1.asm.fa -g Algarrobo -n 30 -S 12345678
# The output including “Step8_*.fa” and “Step8_singlet.fa” were combined and named MITE.lib.
# Second, LTR were identified using LTRharvest and filtered with LTRdigest in the software genometools and other custom programs (CRL scripts). The protocol distinguishes between relatively recent and old LTRs by searching for candidates with 99% or 85% similarity, respectively.
gt suffixerator -db ppafile -indexname ppafileindex -tis -suf -lcp -des -ssp -dna
gt ltrharvest -index ppafileindex -out ppafile.out99 -outinner ppafile.outinner99 -gff3 ppafile.gff99 -minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 -mintsd 5 -maxtsd 5 -motif tgca -similar 99 -vic 10 > ppafile.result99
gt gff3 -sort ppafile.gff99 > ppafile.gff99.sort
gt ltrdigest -trnas eukaryotic-tRNAs.fa ppafile.gff99.sort ppafileindex > ppafile.gff99.dgt
perl CRL_Step1.pl --gff ppafile.gff99.dgt
# The output is the file CRL_Step1_Passed_Elements.txt, which Is further filtered to identify three sources of false positives:  tandem local repeats such as centromeric repeats, local gene clusters derived from recent gene duplications, two other transposable elements located in adjacent regions. It may be necessary to modify the CRL_Step3.pl script if using an updated version of muscle.
perl CRL_Step2.pl --step1 CRL_Step1_Passed_Elements.txt --repeatfile ppafile.out99 --resultfile ppafile.result99 --sequencefile ppafile --removed_repeats CRL_Step2_Passed_Elements.fasta
mkdir fasta_files
mv Repeat_*.fasta fasta_files
mv CRL_Step2_Passed_Elements.fasta fasta_files
cd fasta_files
perl CRL_Step3.pl --directory . --step2 CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25
mv CRL_Step3_Passed_Elements.fasta ..
cd ..
# Next, elements with nested insertions of other LTRs or MITEs are identified and masked.
perl ltr_library.pl --resultfile ppafile.result99 --step3 CRL_Step3_Passed_Elements.fasta --sequencefile ppafile
# The output is lLTR_Only.lib which is combined with MITE.lib to form a library to mask the internal regions of putative LTR elements.
cat lLTR_Only.lib MITE.lib > repeats_to_mask_LTR99.fasta
RepeatMasker -lib repeats_to_mask_LTR99.fasta -nolow -dir . ppafile.outinner99
# Then, the elements with nested insertions are removed.
perl cleanRM.pl ppafile.outinner99.out ppafile.outinner99.masked > ppafile.outinner99.unmasked
# If the internal region of the putative candidate element is very short, it is likely to be a false positive. These elements will be excluded.
perl rmshortinner.pl ppafile.outinner99.unmasked 50 > ppafile.outinner99.clean
# Next, BLASTX is used to identify any additional outinner sequence with nested insertions of putative autonomous DNA transposons and eliminate them.
makeblastdb -in Tpases020812DNA -dbtype prot
blastx -query ppafile.outinner99.clean -db Tpases020812DNA -evalue 1e-10 -num_descriptions 10 -out ppafile.outinner99.clean_blastx.out.txt
perl outinner_blastx_parse.pl --blastx ppafile.outinner99.clean_blastx.out.txt --outinner ppafile.outinner99
# The output is the file passed_outinner_sequence.fasta, which is then processed to build exemplars to use in MAKER.
perl CRL_Step4.pl --step3 CRL_Step3_Passed_Elements.fasta --resultfile ppafile.result99 --innerfile passed_outinner_sequence.fasta --sequencefile ppafile
makeblastdb -in lLTRs_Seq_For_BLAST.fasta -dbtype nucl
blastn -query lLTRs_Seq_For_BLAST.fasta -db lLTRs_Seq_For_BLAST.fasta -evalue 1e-10 -num_descriptions 1000 -out lLTRs_Seq_For_BLAST.fasta.out
makeblastdb -in Inner_Seq_For_BLAST.fasta -dbtype nucl
blastn -query Inner_Seq_For_BLAST.fasta -db Inner_Seq_For_BLAST.fasta -evalue 1e-10 -num_descriptions 1000 -out Inner_Seq_For_BLAST.fasta.out
CRL_Step5.pl –LTR_blast lLTRs_Seq_For_BLAST.fasta.out --inner_blast Inner_Seq_For_BLAST.fasta.out --step3 CRL_Step3_Passed_Elements.fasta --final LTR99.lib --pcoverage 90 --pidentity 80
# The output is LTR99.lib, which contains the exemplar sequences of relatively recent LTRs. For old LTRs, the first step is modified to reduce the similarity to 85%.
gt ltrharvest -index ppafileindex -out ppafile.out85 -outinner ppafile.outinner85 -gff3 ppafile.gff85 -minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 -mintsd 5 -maxtsd 5 -vic 10 > ppafile.result85
# The process afterwards is the same as for the relatively recent LTRs up to getting the file LTR85.lib. Since this last set of elements contain elements collected in LTR99.lib, highly similar exemplars are masked.
RepeatMasker -lib LTR99.lib -dir . LTR85.lib
removed_masked_sequence.pl --masked_elements LTR85.lib.masked --outfile FinalLTR85.lib
# Finally, these two libraries are combined to get the final LTR library.
cat LTR99.lib FinalLTR85.lib > allLTR.lib
# The third step is to collect repetitive sequences from the remaining regions of the genome using RepeatModeler.
cat allLTR.lib MITE.lib > allMITE_LTR.lib
RepeatMasker -lib allMITE_LTR.lib -dir . ppafile
perl rmaskedpart.pl ppafile.masked 50 > umppafile
BuildDatabase -name umppafiledb -engine ncbi umppafile
nohup RepeatModeler -database umppafiledb >& umppafile.out
# The results are in a folder named “RM…”, the sequences are in consensi.fa.classified, some identified and other not. We separate them using a script.
perl repeatmodeler_parse.pl --fasta-file consensi.fa.classified --unknowns repeatmodeler_unknowns.fasta --identities repeatmodeler_identities.fasta
# The unknown repeats are searched against a transposase database for classification to the relevant superfamily.
makeblastdb -in Tpases020812 -dbtype prot
blastx -query repeatmodeler_unknowns.fasta -db Tpases020812 -evalue 1e-10 -num_descriptions 10 -out modelerunknown_blast_results.txt
perl transposon_blast_parse.pl --blastx modelerunknown_blast_results.txt --modelerunknown repeatmodeler_unknown.fasta
# These output files are renamed for easy identification:
mv unknown_elements.txt ModelerUnknown.lib
cat identified_elements.txt repeatmodeler_identities.fasta > ModelerID.lib
# Finally, all repeats are searched against a plant protein database without transposons to exclude gene fragments. For this, we use ProtExcluder. For each library, do the following.
makeblastdb -in alluniRefprexp070416 -dbtype prot
blastx -query ModelerUnknown.lib -db alluniRefprexp070416 -evalue 1e-10 -num_descriptions 10 -out ModelerUnknown.lib_blast_results.txt
perl ProtExcluder.pl -f 50 ModelerUnknown.lib_blast_results.txt ModelerUnknown.lib
# The output is ModelerUnknown.libnoProtFinal. Repeat the same for the others library (ModelerID.lib, MITE.lib, AllLTR.lib). Combine all repeats with identifications to form KnownRepeats.lib.
cat MITE.lib allLTR.lib ModelerID.libnoProtFinal > KnownRepeats.lib
# Combine this library with repeats with unknown classifications to form allRepeats.lib.
cat KnownRepeats.lib ModelerUnknown.libnoProtFinal > allRepeats.lib
# The decision to use KnownRepeats.lib or allRepeats.lib depends on whether one wants to do a relaxed or conservative repeat masking. For N. pallida, allRepeats.lib was used.

############################################### TRANSCRIPTS EVIDENCE: TRANSCRIPTOME #######################
# For transcription evidence, a long-read assembly was performed using RATTLE for each tissue separately including leaf, branch, petiole, root, buds, flowers and fruits
porechop -i input_reads.fastq -o output_reads.fastq
seqkit seq -m 150 output_reads.fastq > output_reads_150.fastq
rattle cluster -i output_reads_150.fastq -t 25 -o .
rattle cluster_summary -i output_reads_150.fastq -c clusters.out > cluster_summary.txt
rattle correct -i output_reads_150.fastq -c clusters.out -t 25 --verbose
rattle polish -i consensi.fq -t 25 --summary --verbose
# Having assembled transcripts to reduce chance of false positives. Also, MAKER uses exonerate to identify intron regions. The problem is that some transcripts may be present in low numbers and are removed as evidence
# MAKER allows to include transcriptomes from close species. In this case, leaf transcriptomes from N. alba and P. cineraria were downloaded from NCBI.

############################################### PROTEIN EVIDENCE: PLANT SWISSPROT #######################
# For protein evidence, protein sequences from the database UniProt – SwissProt restricted to plants were downloaded in .dat format. To parse the sequences, the perl script from SwissKnife varsplic.pl was used.
perl varsplic.pl -input unprot_sprot_plants.dat -check_vsps -crosscheck -error varsplic.err -fasta USP_plants.fasta -which full #Also do without isoforms for functional annotation

############################################### ALTERNATIVE TRANSCRIPTS #######################
# Download transcriptome for N. alba and P. cineraria

############################################### MAKER ANNOTATION: ROUND 1 #######################
# To use MAKER, first we need to create and edit control files. Using maker -CTL will create four files. The file maker_exe.ctl needs to have the correct paths for the programs used. The file maker_opts.ctl has all the configuration for the run.
maker -CTL
# For the first round, modify the following lines:
# genome= 	Add the path to the assembly fasta.
# est=	Add the path to each N. pallida transcriptome. Separate each one with a comma. After each path indicate a code for the tissue for easy visualization later. For example, “/path/to/leaf_transcriptome.fa:EST_Leaf”.
# altest=	Add the path to each close species transcriptome. Similar to est, separate with a comma and add the name of the species of origin.
# protein=	Add the path to the swissprot sequences.
# model_org=	Remove whatever is here since we are not using RepBase
# rmlib=	Add the path to the custom repeat library
# est2genome=	Turn this from 0 to 1
# protein2genome=	Turn this from 0 to 1
# trna=	Turn this from 0 to 1
# cpus=	For the “prosopis” workstation, set to 10
# min_protein=	Change from 0 to 50 to require a minimum protein length
# To start MAKER, use the command “maker” in the same folder where the control files are created. To run MAKER more efficiently, run three different tasks in the same directory. MAKER knows when a contig is being processed and will work on another one.
maker 2> maker1.err &
maker 2> maker2.err &
maker 2> maker3.err &
# To monitor how each run is going, you can print the log continuously using the following command:
tail -f maker1.err

## Train gene models for second annotation round
# After MAKER is finished, go to the output folder and extract the gff and sequences. These will be used to get the gene models for training with SNAP and Augustus.
gff3_merge -n -s -d ppa_v1.asm.fa_master_datastore_index.log > Algarrobo_rnd1.all.maker.noseq.gff
fasta_merge -d ppa_v1.asm.fa_master_datastore_index.log
# To train for SNAP, gene models are first filtered to discard models with an AED higher than 0.25 and coding proteins shorter than 50 amino acids.
maker2zff -x 0.25 -l 50 -d /path/to/ppa_v1.asm.fa_master_datastore_index.log
rename ‘s/genome/Algarrobo_rnd1.zff.length50_aed0.25/g’ *
# Then, the training sequences and annotations are gathered with additional 1 kb upstream and downstream.
fathom Algarrobo_rnd1.zff.length50_aed0.25.ann Algarrobo_rnd1.zff.length50_aed0.25.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
# Next, the training parameters are created, and the HMM file is assembled. This will be used for MAKER to predict gene models using the trained SNAP.
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
hmm-assembler.pl Algarrobo_rnd1.zff.length50_aed0.25 params > Algarrobo_rnd1.zff.length50_aed0.25.hmm
# To train for Augustus, a BUSCO feature is used, so in this case the 2023-01-busco environment needs to be activated. Make sure to index the assembly fasta before starting.
samtools faidx assembly.fasta
# The first step is to gather the regions that contain the mRNA annotations with up to 1 kb on each side.
awk -v OFS="\t" '{if ($3 == "mRNA") print $1, $4, $5 }' ../../../1-Output/Algarrobo_rnd1.all.maker.noseq.gff | while read rna; do scaffold=`echo ${rna} | awk '{ print $1 }'`; end=`cat ../../../Genome/ppa_v1.asm.fa.fai | awk -v scaffold="${scaffold}" -v OFS="\t" '{ if ($1 == scaffold) print $2 }'`; echo ${rna} | awk -v end="${end}" -v OFS="\t" '{ if ($2 < 1000 && (end - $3) < 1000) print $1, "0", end; else if ((end - $3) < 1000) print $1, "0", end; else if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }'; done | bedtools getfasta -fi ../../../Genome/ppa_v1.asm.fa -bed - -fo Algarrobo_rnd1.all.maker.transcripts1000.fasta
# Then, BUSCO will do a reannotation of the extracted regions using BLAST and built-in HMMs for a set of conserved genes. Then, with the “--long" option, BUSCO will optimize the HMM search model to train Augustus and produce a trained HMM for MAKER. In this case, the set of conserved genes is from embryophyta_odb10, and the starting HMM model is from arabidopsis.
busco -i Algarrobo_rnd1.all.maker.transcripts1000.fasta -o Algarrobo_rnd1_maker -l embryophyta_odb10 -m genome -c 30 --long --augustus_species arabidopsis --tar --augustus_parameters=--progress=true
# Once BUSCO is complete, we need to modify the names and the content of some files before copying it to the config path of Augustus.
rename ‘s/BUSCO_Algarrobo_rnd1_maker/Neltuma_pallida/g’ *
sed -i ‘s/BUSCO_Algarrobo_rnd1_maker/Neltuma_pallida/g’ Neltuma_pallida_parameters.cfg
sed -i ‘s/BUSCO_Algarrobo_rnd1_maker/Neltuma_pallida/g’ Neltuma_pallida_parameters.cfg.orig1
mkdir /opt/augustus/config/species/Neltuma_pallida
cp Neltuma_pallida_* /opt/augustus/config/species/Neltuma_pallida/.

############################################### MAKER ANNOTATION: ROUND 2 #######################
# First, in the directory ../1-Output, the gff was separated to obtain the est2genome, protein2genome and repeat annotations. These will be used for the following MAKER runs.
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P “\test2genome” > Algarrobo_rnd1.all.maker.est2genome.gff
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P “\tprotein2genome” > Algarrobo_rnd1.all.maker.protein2genome.gff
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P “\trepeat” > Algarrobo_rnd1.all.maker.repeats.gff
# Copy all control files from the first round and modify the following options:
# est=	Remove the paths
# altest=	Remove the paths
# est_gff=	Add the path to the est2genome gff from the first round
# altest_gff=	Add the path to the cdna2genome gff from the first round
# protein=	Remove the path
# protein_gff=	Add the path to the protein2genome gff from the first round
# rmlib=	Remove the path
# repeat_protein=	Remove the path
# rm_gff=	Add the path to the repeats gff from the first round
# snaphmm=	Add the path to the SNAP HMM file
# augustus_species=	Write “Neltuma_pallida”. As the folder is in the config file, it should be recognized by Augustus.
# max_dna_len=	Change to 300000 to increase memory usage

############################################### MAKER ANNOTATION: ROUND 3 #######################
# For the third annotation run, repeat the training with the new gene models. For Augustus, change the initial species to Neltuma_pallida, and the output as Neltuma_pallida_2.
# Modify the following options in the control file to use the retrained gene models and to keep unsupported gene models, which will later be filtered based on the presence of functional domains.
# snaphmm=	Add the path to the new SNAP HMM file
# augustus_species=	Write “Neltuma_pallida_2”

############################################### MAKER ANNOTATION: ROUND 4 #######################
# Rescue predicted gene models with functional domains identified by InterProScan.
interproscan.sh -appl PfamA -iprlookup -goterms -f stv -i Algarrobo_rnd3.all.maker.proteins.fasta
cat Algarrobo_rnd3.all.maker.proteins.fasta.tsv | cut -f1 | sort | uniq | sed ‘s/abinitio/predicted/’ > names.txt
cat Algarrobo_rnd3.all.maker.gff | grep -f names.txt > Predicted.gff
# Run an extra annotation round passing the final GFF except the predictions. Provide the predictions gff file in the Gene prediction section and keep them.
# maker_gff=	Add the path to the MAKER gff
# est_pass=	Change to 1
# altest_pass=	Change to 1
# protein_pass=	Change to 1
# rm_pass=	Change to 1
# model_pass=	Change to 1
# other_pass=	Change to 1
# est_gff=	Remove the path
# altest_gff=	Remove the path
# protein_gff=	Remove the path
# rm_gff=	Remove the path
# snaphmm=	Remove the path
# augustus_species=	Remove the name
# pred_gff=	Add the path to the predicted gene models with domains
# keep_preds=	Change to 1

############################################### MAKER ANNOTATION: ROUND 5 #######################
# Reassess and improve gene models using all RNASeq data, copy the maker control file from the third round
# est_gff= Change path to all RNASeq data (or merged with transcriptome)
# snaphmm= Remove the path
# augustus_species= Remove the name
# pred_gff= Add the path to all gene models from fourth round
# keep_preds= Change to 1
# Retrieve the GFF and sequence files as explained above and follow the functional annotation

############################################### FUNCTIONAL ANNOTATION #######################
# First get the GFF and fasta files as explained above. Change name of genes
maker_map_ids --prefix NPAL_ --justify 6 Algarrobo_rnd5.all.maker.gff > Algarrobo_rnd5.all.maker.map
map_gff_ids Algarrobo_rnd5.all.maker.map Algarrobo_rnd5.all.maker.gff
map_fasta_ids Algarrobo_rnd5.all.maker.map Algarrobo_rnd5.proteins.fasta
map_fasta_ids Algarrobo_rnd5.all.maker.map Algarrobo_rnd5.transcripts.fasta

# Run InterProScan to identify functional domains and append the matches to the last column of the GFF.
~/Bioprograms/interproscan/interproscan-5.65-97.0/interproscan.sh -appl PfamA -iprlookup -goterms -f tsv -i Algarrobo_rnd5.proteins.fasta
ipr_update_gff Algarrobo_rnd5.all.maker.gff Algarrobo_rnd5.proteins.fasta > Algarrobo_rnd5.all.maker.ipr.gff

# BLAST with SwissProt plants using BLAST and append most similar protein
~/Bioprograms/blast/ncbi-blast-2.2.28+/bin/makeblastdb -in USP_plants_parsed.fasta -input_type fasta -dbtype prot #It has to be protein file without isoform variants
~/Bioprograms/blast/ncbi-blast-2.2.28+/bin/blastp -db USP_plants_parsed.fasta -query Algarrobo_rnd5.proteins.fasta -out maker2uni.blastp -evalue .000001 -outfmt 6 -num_alignments 1 -seg yes -soft_masking true -lcase_masking -max_hsps_per_subject 1 -num_threads 25
maker_functional_gff USP_plants_parsed.fasta maker2uni.blastp Algarrobo_rnd5.all.maker.ipr.gff > Algarrobo_rnd5.all.maker.ipr.blast.gff
maker_functional_fasta USP_plants_parsed.fasta maker2uni.blastp Algarrobo_rnd5.proteins.fasta > Algarrobo_rnd5.proteins.blast.fasta
maker_functional_fasta USP_plants_parsed.fasta maker2uni.blastp Algarrobo_rnd5.transcripts.fasta > Algarrobo_rnd5.transcripts.blast.fasta

# Lastly, annotate with eggNOG mapper
emapper.py -i Algarrobo_rnd5.proteins.blast.fasta --output_dir EM_Viridiplantae --cpu 20 --dbmem --report_orthologs --excel --decorate_gff Algarrobo_rnd5.all.maker.ipr.blast.gff -o Algarrobo_rnd5.all.maker.ipr.blast --tax_scope Viridiplantae

############################################### MAKER ANNOTATION: ROUND 6 #######################
# Additional round to identify alternative splice isoforms. Copy the maker control file from round 5
# snaphmm= Add path to second trained hmm file
# augustus_species= Specify second species name
# pred_gff= Remove the path
# model_gff= Add path to gene models from previous round
# keep_preds= Change to 0
# alt_splice= Change to 1
# Retrieve the GFF and sequence files as explained above, and filter out unsupported gene models
quality_filter.pl -s Algarrobo_rnd6.all.maker.noseq.ipr.blast.emapper.decorated.gff > Algarrobo_rnd6.all.maker.noseq.ipr.blast.emapper.decorated.aed1.gff
grep -P "\tmRNA" Algarrobo_rnd6.all.maker.noseq.ipr.blast.emapper.decorated.aed1.gff | cut -f9 | cut -d ';' -f1 | seqkit grep -f - Algarrobo_rnd6.proteins.fasta > Algarrobo_rnd6.proteins.aed1.fasta
grep -P "\tmRNA" Algarrobo_rnd6.all.maker.noseq.ipr.blast.emapper.decorated.aed1.gff | cut -f9 | cut -d ';' -f1 | seqkit grep -f - Algarrobo_rnd6.transcripts.fasta > Algarrobo_rnd6.transcripts.aed1.fasta

# Then rename and redo the functional annotation as explained above. Be careful that previous named IDs may be contained in this round's "Name" key. Advise to use ID= before the code for greping gene models. Also, when creating map to rename, maker can get confused if using same "NPAL" prefix, try using "Npal" and then sed to change it to caps. Also, add sufix for different transcripts.
maker_map_ids --prefix NPAL_ --suffix _R --iterate 1 --justify 6 2024_Algarrobo_rnd6.genes.aed1.gff > 2024_Algarrobo_rnd6.genes.aed1.map
# After functional annotation, it seems that annotation from previous round persist for some genes. Not changed yet because it may be useful for later comparisons.
# For comparative genomic analyses and statistics, be careful to first remove tRNA genes in GFF. Also keeping only maker models and not repeats and evidence in the GFF to make file lighter and reading faster.
grep "trnascan" Algarrobo_rnd6.genes.aed1.ipr.blast.emapper.decorated.gff | cut -f9 | cut -d';' -f1 > trnaIDs.names
cat trnaIDs.names | grep -v -f - Algarrobo_rnd6.genes.aed1.ipr.blast.emapper.decorated.gff > Algarrobo_rnd6.PCgenes.aed1.ipr.blast.emapper.decorated.gff

############################################### ANNOTATION EVALUATION #######################
# For visualizaiton, separate the GFF in parts
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P “\test2genome” | grep “EST_Tissue” > EST_Tissue.gff # For different tissues
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P “\cdna2genome” | grep “ALTEST_Species” > ALTEST_Species.gff # For different species
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P “\tprotein2genome”> PROT_Swissprot.gff
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P “\trepeat”> REP_All.gff # Can be further separated into Simple repeats, Low complexity, MITEs, LTRs, DNA, LINE, SINE, RC, Others and Unknown. LTRs can be separated in Old and New if named differently in custom repeat library construction
cat Algarrobo_rndN.all.maker.noseq.gff | grep -P “\maker\t” | grep -v “trnascan” > MAKER_rndN.gff
cat Algarrobo_rndN.all.maker.noseq.gff | grep -P “\maker\t” | grep “trnascan” > MAKER_rndN_tRNA.gff
# The genome and all these tracks will be uploaded to IGV for visualization.
igv.sh
# It may be needed to sort and index the genome and some of the tracks. If IGV closes unexpectedly, change the configuration to allocate more memory.
gff3sort.pl --precise --chr_order natural file.gtf/gff | bgzip > file.gtf/gff.gz
tabix -p gff file.gtf/gff.gz
# To count the number of genes and the mean gene length, use the following command. Note that trnascan annotations are being excluded. Can also count the number of tRNAs identified.
cat Algarrobo_rndN.all.maker.noseq.gff | grep -v “trnascan” | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }' # Number of genes and mean gene length
# The distribution of the AED score is a good way to evaluate the annotation results. The more support a gene model has, the lower the AED score is. A good annotation typically has 90% of the gene models with an AED less than 0.5. To get the cumulative fraction of gene models under a given AED score use:
AED_cdf_generator.pl -b 0.025 Algarrobo_rndN.all.maker.noseq.gff # AED distribution
# Finally, a good metric for the annotation quality is the BUSCO completeness of the transcriptome. To do this, in the Busco environment use the following command. A BUSCO completeness over 90% is good.
busco -i Algarrobo_rndN.all.maker.transcripts.fasta -o Algarrobo_rndN_maker -l embryophyta_odb10 -m transcriptome -c 30 --tar --augustus_parameters=--progress=true

############################################### REPEAT ANALYSIS #######################
# Repeats were separated from the general gff into sub-gffs corresponding to MITE, LTR, DNA, RC, LINE, SINE and Other repeats.
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P "\trepeat" | grep -P "match\t" | sed 's/Name=/\t/' | grep -Pv "\trepeatrunner" | sed 's/|genus/\tgenus/' | sed 's/;Target/\tTarget/' | cut -f1-8,10,11 | uniq | grep "genus:Simple_repeat" | cut -f1-9 | sortBed | gff2bed | cut -f1-4,6 > Simple_repeat.bed
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P "\trepeat" | grep -P "match\t" | sed 's/Name=/\t/' | grep -Pv "\trepeatrunner" | sed 's/|genus/\tgenus/' | sed 's/;Target/\tTarget/' | cut -f1-8,10,11 | uniq | grep "genus:Low_complexity" | cut -f1-9 | sortBed | gff2bed | cut -f1-4,6 > Low_complexity.bed
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P "\trepeat" | grep -P "match\t" | sed 's/Name=/\t/' | grep -Pv "\trepeatrunner" | sed 's/|genus/\tgenus/' | sed 's/;Target/\tTarget/' | cut -f1-8,10,11 | uniq | grep "genus:RC" | cut -f1-9 | sortBed | gff2bed | cut -f1-4,6 > RC.bed
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P "\trepeat" | grep -P "match\t" | sed 's/Name=/\t/' | grep -Pv "\trepeatrunner" | sed 's/|genus/\tgenus/' | sed 's/;Target/\tTarget/' | cut -f1-8,10,11 | uniq | grep "genus:LINE" | cut -f1-9 | sortBed | gff2bed | cut -f1-4,6 > LINE.bed
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P "\trepeat" | grep -P "match\t" | sed 's/Name=/\t/' | grep -Pv "\trepeatrunner" | sed 's/|genus/\tgenus/' | sed 's/;Target/\tTarget/' | cut -f1-8,10,11 | uniq | grep "genus:LTR" | cut -f1-9 | sortBed | gff2bed | cut -f1-4,6 > LTR.bed
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P “\trepeat” | grep -P “match\t” | sed 's/Name=/\t/' | grep -Pv “\trepeatrunner” | sed 's/|genus/\tgenus/' | sed 's/;Target/\tTarget/' | cut -f1-8,10,11 | uniq | grep “genus:Unspecified” | grep “dbseq-nr” | cut -f1-9 | sortBed | gff2bed | cut -f1-4,6 >> LTR.bed
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P “\trepeat” | grep -P “match\t” | sed 's/Name=/\t/' | grep -Pv “\trepeatrunner” | sed 's/|genus/\tgenus/' | sed 's/;Target/\tTarget/' | cut -f1-8,10,11 | uniq | grep “genus:Unspecified” | grep “Algarrobo” | cut -f1-9 | sortBed | gff2bed | cut -f1-4,6 > MITE.bed
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P “\trepeat” | grep -P “match\t” | sed 's/Name=/\t/' | grep -Pv “\trepeatrunner” | sed 's/|genus/\tgenus/' | sed 's/;Target/\tTarget/' | cut -f1-8,10,11 | uniq | grep “genus:Unspecified” | grep “|SINE|” | cut -f1-9 | sortBed | gff2bed | cut -f1-4,6 > SINE.bed
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P “\trepeatrunner” | grep -P “match\t” | uniq | cut -f1-9 | sortBed | gff2bed | cut -f1-4,6 > RRunner.bed
# To count the number of repeats and the base pairs covered in each scaffold, we need to create a general bed file of the assembly. Then, we use bedmap for each repeat bed file separately.
cat Algarrobo_rnd1.all.maker.noseq.gff | grep -P "\tcontig" | sortBed | gff2bed | cut -f1,2,3,4,6 > ContigSizes.bed
bedmap --echo --ec --bases-uniq --bases-uniq-f --count ContigSizes.bed repeat.bed > repeat.info
# To count the number of repeats and bases covered in general:
cat Simple_repeats.info | sed 's/|/\t/g' | cut -f6 | awk '{s+=$1} END {print s}' #Bases covered
cat Simple_repeats.info | sed 's/|/\t/g' | cut -f8 | awk '{s+=$1} END {print s}' #Number of repeats
# To obtain these measures for all repeats combined, concatenate all bed files and then sort appropriately.
cat Simple_repeats.bed Low_complexity.bed LTR.bed LINE.bed SINE.bed DNA.bed MITE.bed RC.bed RRunner.bed > All.bed
cat All.bed | sort -k1,1 -k2,2n -k3,3n > All_sortedUnix.bed
cat ContigSizes.bed | sort -k1,1 -k2,2n -k3,3n > ContigSizes_sortedUnix.bed
bedmap --echo --ec --bases-uniq --bases-uniq-f --count ContigSizes_sortedUnix.bed All_sortedUnix.bed

############################################### TRANSCRIPTOME ANALYSIS #######################
# For each rattle output, a summary file was generated using the following:
cat transcriptome.fq | grep "@cluster" | cut -d' ' -f2 | sort | uniq -c >> Custom_Summ.txt
cat transcriptome.fq | grep "@cluster" | cut -d' ' -f3 | cut -d'=' -f2 | awk '{sum += $1} END {print sum}' >> Custom_Summ.txt
X='tail -n1 Custom_Summ.txt'
cat transcriptome.fq | grep "@cluster" | cut -d' ' -f1,3 | sed 's/total_reads=//' | sed 's/ /\t/' | awk '{$2 /= $X}1' | sed 's/ /\t/' | sort -gr -k2,2 > Proportions.txt
# Next, filter the GFF file to get the proportions of aligned transcripts and sort them in decreasing order.
cat Algarrobo_rnd1.all.maker.est2genome.gff | grep "EST_Tissue" | grep -P "\texpressed_sequence_match" | cut -f9 | sed 's/Name=/\t/' | sed 's/;Target/\t/' | cut -f2 | sort | uniq -c | sed 's/ /\t/g' | cut -f7,8 > Tissue_GFFpresence.txt
while read line; do grep -w $line Proportions_Tissue.txt; done < <(cut -f2 Tissue_GFFpresence.txt) | sort -gr -k2,2 > Tissue_GFFproportions.txt
# Additionally, modify the names of the aligned transcripts in the original GFF or the visualization tracks to reflect the tissue of origin, since different tissues could have the same cluster name for different transcripts.
sed -i 's/cluster_/Tissue_cluster_/g' Algarrobo_rnd3.all.maker.noseq.gff

############################################### PSEUDOGENE IDENTIFICATION #######################
# To update with last annotation round gff
# First, the regions covering exons in the final gene models are hard masked so the pseudogenes are searched in introns and intergenic regions.
cat Algarrobo_rnd3.filter.gff | grep -P “\texon” > Exons.gff
bedtools maskfasta -fi Algarrobo_genome.fa -bed Exons.gff -fo Masked_algarrobo.fa
# Then, a blast search using the protein sequences. Make sure to use the proper version of blast (PATH).
makeblastdb -in Masked_algarrobo.fa -dbtype nucl
tblastn -query Algarrobo_rnd3.filter.proteins.fasta -db Masked_algarrobo.fa -out Algarrobo_prot_tblastn6 -seg no -outfmt 6 -num_threads 15
# Change the test_parameter_file in the same directory to include the file names and paths.
# p_seq= 	Add the name of the protein sequences file
# g_seq=	Add the name of the masked genome file
# b_out=	Add the name of the blast results file
# p_codes=	Add the full path to the scripts directory of the Pseudogene Pipeline
# f_dir=	Add the full path to the bin of the Pseudogene conda environment
# blosum=	Add the full path to the blosum matrix sample file
# gff=	Add the full path to the full GFF annotation file
# r_dir=	Add the full path to the bin of the Pseudogene conda environment
# Once all files are in the correct directory and the parameter file is adjusted, run the pipeline.

############################################### ANNOTATION STATISTICS #######################
# Number of gene models without UTR: 7,046
cat 2024_Algarrobo_rnd6.all.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | sed 's/_QI=/\t/' | cut -f2 | cut -d';' -f1 | cut -d'|' -f1,8 | sed 's/|/\t/'| awk '$1==0 && $2==0' | wc -l
# Number of gene models with both UTRs: 24,049
cat 2024_Algarrobo_rnd6.all.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | sed 's/_QI=/\t/' | cut -f2 | cut -d';' -f1 | cut -d'|' -f1,8 | sed 's/|/\t/'| awk '$1>0 && $2>0' | wc -l
# Number of gene models with 5UTR: 27,836
cat 2024_Algarrobo_rnd6.all.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | sed 's/_QI=/\t/' | cut -f2 | cut -d';' -f1 | cut -d'|' -f1,8 | sed 's/|/\t/'| awk '$1>0' | wc -l
# Number of gene models with 3UTR: 29,358
cat 2024_Algarrobo_rnd6.all.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | sed 's/_QI=/\t/' | cut -f2 | cut -d';' -f1 | cut -d'|' -f1,8 | sed 's/|/\t/'| awk '$2>0' | wc -l
# Mean 5UTR length: 318.738
cat 2024_Algarrobo_rnd6.all.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | sed 's/_QI=/\t/' | cut -f2 | cut -d';' -f1 | cut -d'|' -f1 | awk '$1>0' | awk '{s+=$0}END{print s/NR}'
# Mean 3UTR length: 465.764
cat 2024_Algarrobo_rnd6.all.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | sed 's/_QI=/\t/' | cut -f2 | cut -d';' -f1 | cut -d'|' -f8 | awk '$1>0' | awk '{s+=$0}END{print s/NR}'
# Number gene models with exons covered by EST: 33,893
cat 2024_Algarrobo_rnd6.all.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | sed 's/_QI=/\t/' | cut -f2 | cut -d';' -f1 | cut -d'|' -f3 | awk '$1>0' | wc -l
# Number of gene models with exons covered by EST, or EST or Protein: 39,964
cat 2024_Algarrobo_rnd6.all.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | sed 's/_QI=/\t/' | cut -f2 | cut -d';' -f1 | cut -d'|' -f3,4,6 | sed 's/|/\t/g' | awk '$1>0 || $2>0' | wc -l
# Average fraction of exons covered by EST: 0.84293
cat 2024_Algarrobo_rnd6.all.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | sed 's/_QI=/\t/' | cut -f2 | cut -d';' -f1 | cut -d'|' -f3 | awk '$1>0' | awk '{s+=$0}END{print s/NR}'
# Average fraction of exons covered by EST or Protein: 0.951426
cat 2024_Algarrobo_rnd6.all.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | sed 's/_QI=/\t/' | cut -f2 | cut -d';' -f1 | cut -d'|' -f4 | awk '$1>0' | awk '{s+=$0}END{print s/NR}'
# Number of gene models with more than 1 exon: 38,274
cat 2024_Algarrobo_rnd6.all.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | sed 's/_QI=/\t/' | cut -f2 | cut -d';' -f1 | cut -d'|' -f2,7 | sed 's/|/\t/g' | awk '$2>1' | cut -f1 | wc -l
# Number of gene models with more than 1 exon with splice sites covered by EST: 32,574
cat 2024_Algarrobo_rnd6.all.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | sed 's/_QI=/\t/' | cut -f2 | cut -d';' -f1 | cut -d'|' -f2,7 | sed 's/|/\t/g' | awk '$2>1' | cut -f1 | awk '$1 > 0 '| wc -l
# Average fraction of splice sites covered by EST: 0.858279
cat 2024_Algarrobo_rnd6.all.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | sed 's/_QI=/\t/' | cut -f2 | cut -d';' -f1 | cut -d'|' -f2,7 | sed 's/|/\t/g' | awk '$2>1' | cut -f1 | awk '$1 > 0 '| awk '{s+=$0}END{print s/NR}'
## Functional annotation statistics
# Number of gene models annotated with InterPro
cat 2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | grep "Dbxref=InterPro" | wc -l
# Print names of gene models annotated with InterPro
cat 2024_Algarrobo_rnd5.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | grep-v "Dbxref=InterPro" | cut -d'-' -f1 | sed 's/ID=//' > Functional_Interpro.names
# Number of gene models annotated with SwissProt
cat 2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | grep -v "Note=Protein of unknown function" | wc -l
# Print names of gene models annotated with SwissProt
cat 2024_Algarrobo_rnd5.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | grep -v "Note=Protein of unknown function" | cut -d'-' -f1 | sed 's/ID=//' > Functional_BLAST.names
# Number of gene models annotated with EggNOG
cat 2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | grep "em_target" | wc -l
# Print names of gene models annotated with EggNOG
cat 2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | grep "em_target" | cut -d'-' -f1 | sed 's/ID=//' > Functional_EggNOG.names
# Functional annotation overlaps were drawn in a Venn diagram
# Number of gene models with COG categories: 36,031
cat 2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "em_COG" | sed 's/em_COG_cat=/#/' | cut -d'#' -f2 | cut -d';' -f1 | grep -v "None" | wc -l
# Get list of COG categories for all gene models
cat 2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "em_COG" | sed 's/em_COG_cat=/#/' | cut -d'#' -f2 | cut -d';' -f1 | grep -v "None" | grep -o . | sort | uniq -c | sort -rn | head
# Get matches to SwissProt proteins
cat 2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "Note=Similar to" | cut -f9 | sed 's/Note=/#/' | cut -d'#' -f2 | cut -d';' -f1 | sed 's/(Fragment/(fragment/' | sed 's/([A-Z][a-z][a-z]/\t/' | cut -f2 | cut -d')' -f1 | sort | uniq -c | sort -rn
# Get KEGG pathway list for all gene models
cat 2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "em_KEGG_Pathway" | sed 's/em_KEGG_Pathway=/\t/' | cut -f10 | cut -d';' -f1 | sed 's/,/\n/g' | sort | uniq -c | sort -rn
# Get PFAM list for all gene models
cat 2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "PFAM:" | sed 's/PFAM:/\t/' | cut -f10 | cut -d';' -f1 | sed 's/PFAM://g' | sed 's/,/\n/g' | sort | uniq -c | sort -rn
# Get number of protein kinase domains identified with EggNOG
cat 2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "em_PFAMs=" | sed 's/em_PFAMs=/\t/' | cut -f10 | cut -d';' -f1 | grep -w "Pkinase" | wc -l
# Get number of protein kinase tyr domains identified with EggNOG
cat 2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "em_PFAMs=" | sed 's/em_PFAMs=/\t/' | cut -f10 | cut -d';' -f1 | grep -w "Pkinase_Tyr" | wc -l
# Get number of LRR domains identified with EggNOG
cat 2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "em_PFAMs=" | sed 's/em_PFAMs=/\t/' | cut -f10 | cut -d';' -f1 | grep "LRR" | wc -l
# Get number of PPR domains identified with EggNOG
cat 2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "em_PFAMs=" | sed 's/em_PFAMs=/\t/' | cut -f10 | cut -d';' -f1 | grep "PPR" | wc -l
# Number of gene models involved in carbon fixation modules
cat 2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "map00710" | wc -l
# Number of genes involved in carbon fixation modules
cat 2024_Algarrobo_rnd6.PCgenes.aed1.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "map00710" | cut -f9 | cut -d'_' -f1,2 | sort | uniq | wc -l
# Number of gene models involved in biosynthesis of secondary metabolites: 2,786
cat 2024_Algarrobo_rnd6.genes.maker.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "map01110" | wc -l
# Number of genes involved in biosynthesis of secondary metabolites: 1,781
cat 2024_Algarrobo_rnd6.PCgenes.aed1.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "map01110" | cut -f9 | cut -d'_' -f1,2 | sort | uniq | wc -l
# Number of gene models also involved in biosynthesis of secondary metabolites: 23
cat 2024_Algarrobo_rnd6.PCgenes.aed1.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "map00940" | grep -v "map01110" | cut -f9 | cut -d'_' -f1,2 | sort | uniq | wc -l
# Number of genes involved in signal transduction pathways: 347
cat 2024_Algarrobo_rnd6.PCgenes.aed1.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "map04075" | cut -f9 | cut -d'_' -f1,2 | sort | uniq | wc -l
# Number of genes involved in plant-pathogen interaction: 317
cat 2024_Algarrobo_rnd6.PCgenes.aed1.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "map04626" | cut -f9 | cut -d'_' -f1,2 | sort | uniq | wc -l
# Number of genes involved in MAPK signaling pathways: 225
cat 2024_Algarrobo_rnd6.PCgenes.aed1.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "map04016" | cut -f9 | cut -d'_' -f1,2 | sort | uniq | wc -l
## MAKER evidence overlap
cat 2024_Algarrobo_rnd6.PCgenes.aed1.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | cut -d';' -f1,6 | sed 's/ID=//' | sed 's/-RA;_QI=/\t/' | sed 's/|/\t/g' > Genes_QImaker.tsv
# Names of gene models with exons covered by EST or Protein
cat Genes_QImaker.tsv | awk '$4 > 0' | cut -f1 > Exons_EST_Prot_cov.names
# Names of gene models with exons covered by Prediction
cat Genes_QImaker.tsv | awk '$6 > 0' | cut -f1 > Exons_Prediction_cov.names
# GFF of genes overlapped with alternative EST evidence
bedtools intersect -wa -a 2024_Algarrobo_rnd6.PCgenes.aed1.ipr.blast.emapper.decorated.gff -b Algarrobo_rnd1.cdna2genome.gff > Overlap_AltEST_bedtools.bed
# Names of gene models covered by alternative EST evidence
cat Overlap_AltEST_bedtools.bed | grep "ID=NPAL" | cut -f9 | grep "RA" | sed 's/-RA/\t/' | cut -f1 | cut -d';' -f1 | sed 's/ID=//' | sort | uniq > Exons_AltEST_cov.names
############################################### FAMILY ANALYSIS #######################
grep -P "CYP[0-9]" ../2024_Algarrobo_rnd6.genes.aed1.ipr.blast.emapper.decorated.gff | grep -v "Peptidyl-prolyl" | grep -P "\tmRNA" | cut -f9 | cut -d'_' -f1,2 | sort | uniq | grep -f - ../2024_Algarrobo_rnd6.genes.aed1.ipr.blast.emapper.decorated.gff > CYP450s.gff
