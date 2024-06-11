############################################### HYBRID ASSEMBLY #######################
# All nanopore raw fast5 files were basecalled using Guppy with SUP model
guppy_basecaller -r -I Raw_reads.fast5 -s Raw_reads.fastq –config dna_r9.4.1_450bps_sup.cfg -q 0 --trim_strategy dna --disable_trim_barcodes --compress_fastq --calib_detect --device ‘auto’ &> stdout.txt
# Reads were filtered to retain reads longer than 10 kb only. Only reads in the pass folder were used, that is with quality > 9
seqkit seq -m 10000 Raw_reads.fastq > Filtered_reads.fastq

## Draft assembly using Flye. Polishing step deactivated, using Medaka and pilon for that
flye --nano-hq Filtered_reads.fastq --out-dir HQ_10k --threads 30 -i 0

## Polishing
# First, the assembly was corrected using all nanopore reads with Medaka
medaka_consensus -d assembly.fasta -l Raw_reads.fastq -m r941_min_sup_g507 -o Med_output_1 -t 20
medaka_consensus -d assembly_Med1.fasta -l Raw_reads.fastq -m r941_min_sup_g507 -o Med_output_2 -t 20
# Second, short reads were used to polish the corrected assembly
# Mapping with subsample (~50X coverage with assembly length 456357777), otherwise it takes too long
seqtk sample -s100 Short_reads_1.fq 75555924 > Forward50X.fastq
seqtk sample -s100 Short_reads_2.fq 75555924 > Reverse50X.fastq
bwa index assembly_Med2.fasta
bwa mem -t 30 assembly_Med2.fasta Forward50X.fastq Reverse50X.fastq | samtools sort -@ 30 -o bwa_mapping_sort.bam
samtools index bwa_mapping_sort.bam
# Polish with Pilon. With "--fix" it only scans for SNPs and indels
pilon --genome assembly_Med2.fasta --frags bwa_mapping_sort.bam --changes --diploid --fix bases -Xmx30G --output assembly_Med2_Pil1.fasta --outdir Pilon_1 2>&1
pilon --genome assembly_Med2_Pil1.fasta --frags bwa_mapping_sort.bam --changes --diploid --fix bases -Xmx30G --output assembly_Med2_Pil2.fasta --outdir Pilon_2 2>&1

## BUSCO assessment
busco assembly_Med2_Pil2.fasta -db embryophyta_odb10
seqstats assembly_Med2_Pil2.fasta

## Haplotig removal. Contigs shorter than 50 kb were removed
seqkit seq -m 50000 assembly_Med2_Pil2.fasta > assembly_filt.fasta
# Haplotigs were removed using long reads with purge_dups in 1 round
pd_config.py assembly_filt.fasta Filtered_Reads.fastq
run_purge_dups.py -p bash config.json purge_dups/src/ Npallida
hist_plot.py -c cutoffs -t "Read depth histogram" PB.stat Figure.png

## Hi-C Scaffolding
# Was done in UGA with Juicebox and resulted in 14 chromosomes + 67 unanchored scaffolds
