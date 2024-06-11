############################################### KMER ANALYSIS #######################
## Create kmer histograms
kat hist -t 60 -o kat_hist_SR ../86-SA1_back_up_?.fastq.gz
kat plot spectra-hist -x 300 -o kat_hist_SR.png kat_hist_SR
kat hist -t 60 -o kat_hist_LR -H 5000000000 ../All_Nanopore_lib12.fastq
kat plot spectra-hist -x 300 -o kat_hist_LR.png kat_hist_LR
## Detect GC bias. For clean data you see a clean circle at the expected GC level
kat gcp -t 50 -o kat_gcp '../86-SA_back_up_?.fastq.gz' #Short reads
kat gcp -n -t 80 -o kat_gcp_long ../All_Nanopore_lib12.fastq -J 5000000000 &> kat_gcp_long.log #Long reads
## Compare the k-mer spectra of read 1 to read 2 files of a PE sequencing run.
# SR1 vs SR2
kat comp -t 40 -n -o kat_comp ../86-SA1_back_up_?.fastq.gz # "-n": makes a density plot
# Generate spectra plot
kat plot spectra-mx -i -o kat_plot_spectramx kat_comp-main.mx # "-i" intersection mode, plots shared and exclusive content found in the matrix
# SR vs LR
kat comp -n -t 80 -o kat_comp_IlvsNa '../86-SA1_back_up_?.fastq.gz' '../All_Nanopore_lib12.fastq' -J 5000000000 &> kat_comp_IlvsNa.log
# Generate spectra plot
kat plot spectra-mx -i -o kat_comp_IlvsNa_shared_spectramx kat_comp_IlvsNa-main.mx
## Detect contaminants in assembly
# Mark contigs with average k-mer coverage and GC%
kat sect -t 50 -o kat_sect_SR ../../Sequences/ppa_v1.asm.fasta ../86-SA1_back_up_?.fastq.gz
#Only long reads
kat sect -t 50 -o kat_sect_LR ../../Sequences/ppa_v1.asm.fasta ../All_Nanopore_lib12.fastq -H 5000000000
#Including Long reads
kat sect -t 50 -o kat_sect_SR_LR ../../Sequences/ppa_v1.asm.fasta ../86-SA1_back_up_?.fastq.gz ../All_Nanopore_lib12.fastq -H 5000000000
## Check assembly coherence
kat comp -t 50 -o SR_vs_assembly_comp '../86-SA1_back_up_?.fastq.gz' ../../Sequences/ppa_v1.asm.fasta
# We should see a single peak with kmers appearing once
kat comp -t 50 -o LR_vs_assembly_comp '../All_Nanopore_lib12.fastq' ../../Sequences/ppa_v1.asm.fasta #Using LR
# We should see a single peak with kmers appearing once
## Distribution decomposition analysis
kat_distanalysis.py --plot spectra-cn.mx
## Finding repetitive regions in assemblies
kat sect -E -F ../../Sequences/ppa_v1.asm.fasta ../../Sequences/ppa_v1.asm.fasta
############################################### PLOIDY ESTIMATION #######################
# First calculate kmer frequencies and then generate a histogram of kmers
kmc -k21 -t90 -m500 -ci1 -cs10000 @FILES kmcdb tmp
kmc_tools transform kmcdb histogram kmcdb_k21.hist -cx10000
# Choose lower and upper coverage thresholds
L=$(smudgeplot.py cutoff kmcdb_k21.hist L)
U=$(smudgeplot.py cutoff kmcdb_k21.hist U)
# Extract kmers in the coverage range using kmc_dump
kmc_dump kmcdb_k21.hist -ci"$L" -cx"$U" dump -o kmcdb_L"$L"_U"$U".dump
kmc_tools transform kmcdb -ci"$L" -cx"$U" dump -s kmcdb_L"$L"_U"$U".dump
# Compute the set of kmer pairs
smudgeplot.py hetkmers -o kmcdb_L"$L"_U"$U" < kmcdb_L"$L"_U"$U".dump
# Report haploid kmer coverage reported by GenomeScope or let it be estimated ofr generating the smudgeplot
smudgeplot.py plot kmcdb_L"$L"_U"$U"_coverages.tsv

## For GenomeScope, copy the hist file to the webpage, select kmer 21, ploidy 2 and coverage limit 10000

############################################### COMPARATIVE GENOMICS: ORTHOLOGY INFERENCE #######################
## Download GFF and protein data 2024 from NCBI (P. cin, N. alb, S. tor, B. var) and Ensembl (G. max, M. tru, P. per, V. vin)
# Extract the protein coding genes only. For NCBI data (B.var):
cat 1-Genomes/Bvar/GCA_022379115.2.gff | grep -v "#" | grep "protein_coding" | cut -f9 | sed 's/Name=/\t/' | cut -f2 | cut -d';' -f1 | grep -f - 1-Genomes/Bvar/GCA_022379115.2.gff > PCG_Bvar.gff
# Then generate a bed file for primary transcripts (first one appearing in gff)
python3 -m jcvi.formats.gff bed --type=mRNA --key=Parent --primary_only PCG_Bvar.gff -o PCG_Bvar.bed
# Catch the peptide sequences with the bed names
cat PCG_Bvar.bed | cut -f4 | grep -f - PCG_Bvar.gff | grep -P "\tCDS" | cut -f9 | sed 's/protein_id=/\t/' | cut -f2 | cut -d';' -f1 | sort | uniq | seqkit grep -f - 1-Sequences/Bauhinia_variegata.GCA_022379115.2.pep.fa > PCG_Bvar.pep
# Repeat the process for S. tor, N alb and P cin
## For Ensembl data (G. max, M. tru, P. per, V. vin).
# Get only protein coding genes
cat 1-Genomes/Gmax/Glycine_max.Glycine_max_v2.1.58.gff3 | grep -v "#" | grep "protein_coding" | cut -f9 | sed 's/Name=/\t/' | cut -f2 | cut -d';' -f1 | grep -f - 1-Genomes/Gmax/Glycine_max.Glycine_max_v2.1.58.gff3 > PCG_Gmax.gff
# Catch canonical transcripts (representatives)
cat PCG_Gmax.gff | grep -v "coding;transcript" > PCG_Gmax_canonical.gff
# Generate bed file with primary transcripts
python3 -m jcvi.formats.gff bed --type=mRNA --primary_only PCG_Gmax_canonical.gff -o PCG_Gmax_canonical.bed
# Catch peptide sequences with matching names in bed file
cat PCG_Gmax_canonical.bed | cut -f4 | sed 's/transcript://' | seqkit grep -f - 1-Sequences/Glycine_max.Glycine_max_v2.1.pep.all.fa > PCG_Gmax_canonical.pep
## Repeat the process for M. tru, P. per and V. vin
##For Neltuma pallida:
# Get longest transcript
grep ">" 2024-06-02_Npal_rnd3.merged.proteins.fasta | cut -d' ' -f1 | sed 's/>//' > Npal_transcript.names
cat Npal_transcript.names | cut -d'_' -f1,2 | sort | uniq > Npal_genes.names
grep ">" 2024-06-02_Npal_rnd3.merged.proteins.fasta | sed 's/QI:/\t/' | cut -f2 | sed 's/|/\t/g' | paste Npal_transcript.names - > Npal_QImetrics.txt
grep ">" ../1-Genomes/Npal/2024-05-30_Npal_rnd3.merged.proteins.blast.fasta | sed 's/QI:/\t/' | cut -f2 > Npal_QImetrics.txt
while read line; do grep $line Npal_QImetrics.txt | sort -k10,10 -rn | cut -f1 | head -n1 >> Npal_longesttranscript.txt; done < Npal_genes.names
cat Npal_longesttranscript.txt | seqkit grep -f - 2024-06-02_Npal_rnd3.merged.proteins.fasta > ../../NewGFF/PCG_Npal.pep
cat Npal_longesttranscript.txt | seqkit grep -f - 2024-06-02_Npal_rnd3.merged.transcripts.fasta > ../../NewGFF/PCG_Npal.dna
cat Npal_longesttranscript.txt | grep -f - 2024-06-02_Npal_rnd3.merged.sorted.gff > ../../NewGFF/PCG_Npal.gff
python3 -m jcvi.formats.gff bed --type=mRNA --primary_only PCG_Npal.gff -o PCG_Npal.bed

### Long transcripts only
# ---------------------------------------------------
# |Results from dataset embryophyta_odb10            |
# ---------------------------------------------------
# |C:88.8%[S:85.6%,D:3.2%],F:3.3%,M:7.9%,n:1614      |
# |1433    Complete BUSCOs (C)                       |
# |1381    Complete and single-copy BUSCOs (S)       |
# |52    Complete and duplicated BUSCOs (D)          |
# |54    Fragmented BUSCOs (F)                       |
# |127    Missing BUSCOs (M)                         |
# |1614    Total BUSCO groups searched               |
# ---------------------------------------------------

# [1] eggNOG-mapper v2: functional annotation, orthology assignments, and domain
#       prediction at the metagenomic scale. Carlos P. Cantalapiedra,
#       Ana Hernandez-Plaza, Ivica Letunic, Peer Bork, Jaime Huerta-Cepas. 2021.
#       Molecular Biology and Evolution, msab293, https://doi.org/10.1093/molbev/msab293
#
# [2] eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated
#       orthology resource based on 5090 organisms and 2502 viruses. Jaime
#       Huerta-Cepas, Damian Szklarczyk, Davide Heller, Ana Hernandez-Plaza,
#       Sofia K Forslund, Helen Cook, Daniel R Mende, Ivica Letunic, Thomas
#       Rattei, Lars J Jensen, Christian von Mering and Peer Bork. Nucleic Acids
#       Research, Volume 47, Issue D1, 8 January 2019, Pages D309-D314,
#       https://doi.org/10.1093/nar/gky1085
#
# [3] Sensitive protein alignments at tree-of-life scale using DIAMOND.
#        Buchfink B, Reuter K, Drost HG. 2021.
#        Nature Methods 18, 366–368 (2021). https://doi.org/10.1038/s41592-021-01101-x

## Copy all peptide sequences file to a new directory
mkdir NewPeps
cp *.pep Peps/
# Change pep files names to full species_names.fa for subsequent analyses

############################################### COMPARATIVE GENOMICS: ORTHOLOGY INFERENCE #######################
# An initial orthofinder run to get the single copy genes
orthofinder -t 60 -a 60 -f NewPeps/
for SCgene in *.fa; do mafft --quiet $SCgene > ${SCgene%.fa}.aligned.fa; done
for SCgene in *aligned.fa; do seqkit fx2tab $SCgene | cut -f2 > ${SCgene%.fa}.seq.tab ; done
seqkit fx2tab OG0015336.aligned.fa | cut -f1 > Headers.txt
paste *.seq.tab | sed 's/\t//g' > Pasted_seqs.txt
paste Headers.txt Pasted_seqs.txt | seqkit tab2fx > SCgenes.merged.fa
trimal -in SCgenes.merged.fa -out SCgenes.merged.trimed.fa -automated1
# Obtain a ML tree
raxml-ng --msa SCgenes.merged.trimed2.fa --model PROTGTR+G --prefix SCG3 --threads 38 --seed 1 --outgroup Vitis_vinifera
# Copy the tree file RAxML_bipartitions.RaxMLTree and use it for a new orthofinder run
orthofinder -t 60 -a 60 -s RAxML_bipartitions.RaxMLTree -f NewPeps/ #"-s" will take into account the tree for the orthology inference
# Apart from the statistics table (which is different from the HOG statistics used for CAFE analysis), upsetR was used to give a graphical output of shared and unique genes. Script is 3-Orthogroup_analysis.R

############################################### COMPARATIVE GENOMICS: GENE FAMILY DYNAMICS #######################
# For CAFE analysis, two inputs are needed: the rooted ultrametric tree and the filtered gene counts.
# Reformat the HOG table and remove orthogroups specific for 1 species and with >100 copies. The sript 3-Orthogroup_analysis.R was used
sed 's/\r$//' hog_gene_counts.tsv > hog_gene_counts.formated.tsv # This may be needed when copying from Windows
# For making the ultrametric tree, root node age was selected based on literature
make_ultrametric.py -r 125 RAxML_bipartitions.RaxMLTree
# Run CAFE analysis
cafe5 -c 50 -i hog_gene_counts.formated.tsv -t SCG3.raxml.bestTree.ultrametric.tre
# Summarize and plot results
cafeplotter -i results/ -o cafeplotter_results --ignore_branch_length
# Count changes significant at 0.05
printf "Sp\tSign\tSigExp\tSigCon\tExp\tCont\tNoChange\n"; while read line;
  do
  a=$(cat result_summary.tsv | awk -v sp=$line '$2 == sp' | awk '$5 < 0.05' | wc -l); \
  b=$(cat result_summary.tsv | awk -v sp=$line '$2 == sp' | awk '$5 < 0.05' | awk '$4 > 0' | wc -l); \
  c=$(cat result_summary.tsv | awk -v sp=$line '$2 == sp' | awk '$5 < 0.05' | awk '$4 < 0' | wc -l); \
  d=$(cat result_summary.tsv| awk -v sp=$line '$2 == sp' | awk '$4 > 0' | wc -l); \
  e=$(cat result_summary.tsv | awk -v sp=$line '$2 == sp' | awk '$4 < 0' | wc -l); \
  f=$(cat result_summary.tsv | awk -v sp=$line '$2 == sp' | awk '$4 == 0' | wc -l); \
  printf "$line\t$a\t$b\t$c\t$d\t$e\t$f\n"; \
done < <(cat result_summary.tsv | sed '1d' | cut -f2 | sort | uniq)
# Sp      Sign    SigExp  SigCon  Exp     Cont    NoChange
# 11      31      1       30      36      102     1653
# 12      33      15      18      255     81      1455
# 13      6       3       3       12      7       1772
# 14      193     172     21      367     183     1241
# 15      83      43      40      88      230     1473
# 16      186     157     29      334     123     1334
# 17      79      69      10      156     82      1553
# Bauhinia_variegata      201     142     59      432     501     858
# Glycine_max     528     467     61      887     244     660
# Medicago_truncatula     320     228     92      394     599     798
# Neltuma_alba    1352    1252    100     1260    111     420
# Neltuma_pallida 1122    94      1028    96      1151    544
# Prosopis_cineraria      141     118     23      320     364     1107
# Prunus_persica  72      52      20      174     1184    433
# Senna_tora      289     82      207     186     872     733
# Vitis_vinifera  71      64      7       277     870     644

printf "Sp\tOG\tCopies\n"; for i in $(seq 2 18); do a=$(cut -f $i Base_count.tab | head -n1); b=$(cut -f $i Base_count.tab | sed '1d' | awk '$1 > 0' | wc -l); c=$(cut -f $i Base_count.tab | sed '1d' | awk '{s+=$1}END{print s}'); printf "$a\t$b\t$c\n"; done
# Sp      OG      Copies
# Medicago_truncatula<1>  14865   24600
# Glycine_max<2>  15336   37632
# Neltuma_pallida<3>      13133   17995
# Neltuma_alba<4> 14776   30906
# Prosopis_cineraria<5>   14713   23307
# Senna_tora<6>   13915   19731
# Bauhinia_variegata<7>   15502   27701
# Prunus_persica<8>       15781   19927
# Vitis_vinifera<9>       17772   24715
# <10>    17772   25706
# <11>    17772   25208
# <12>    16930   24557
# <13>    16663   24297
# <14>    15864   25171
# <15>    15898   22458
# <16>    15432   22982
# <17>    15153   22935


### Significant orthogroups within Prosopis sl group including Npal not reducing significantly
cat ../Base_branch_probabilities.tab | sed '1d' | awk '$17 <= 0.05' | cut -f1 | grep -f - ../Base_change.tab | awk '$17 >= 0' | cut -f1 > 16SigExpanded.txt
cat ../Base_branch_probabilities.tab | sed '1d' | awk '$18 <= 0.05' | cut -f1 | grep -f - ../Base_change.tab | awk '$18 < 0' | cut -f1 | grep -vf - 16SigExpanded.txt > 16Sig17noSigExpanded.txt
cat ../Base_branch_probabilities.tab | sed '1d' | awk '$4 <= 0.05' | cut -f1 | grep -f - ../Base_change.tab | awk '$4 < 0' | cut -f1 | grep -vf - 16Sig17noSigExpanded.txt > 16Sig17noSigNpalnoSigExpanded.txt

cat ../Base_branch_probabilities.tab | sed '1d' | awk '$18 <= 0.05' | cut -f1 | grep -f - ../Base_change.tab | awk '$18 >= 0' | cut -f1 > 17SigExpanded.txt
cat ../Base_branch_probabilities.tab | sed '1d' | awk '$4 <= 0.05' | cut -f1 | grep -f - ../Base_change.tab | awk '$4 < 0' | cut -f1 | grep -vf - 17SigExpanded.txt > 17SigNpalnoSigExpanded.txt

cat ../Base_branch_probabilities.tab | sed '1d' | awk '$4 <= 0.05' | cut -f1 | grep -f - ../Base_change.tab | awk '$4 >= 0' | cut -f1 > NpalSigExpanded.txt

cat ../Base_branch_probabilities.tab | sed '1d' | awk '$5 <= 0.05' | cut -f1 | grep -f - ../Base_change.tab | awk '$5 >= 0' | cut -f1 > NalbSigExpanded.txt

cat ../Base_branch_probabilities.tab | sed '1d' | awk '$6 <= 0.05' | cut -f1 | grep -f - ../Base_change.tab | awk '$6 >= 0' | cut -f1 > PcinSigExpanded.txt

cat 16Sig17noSigNpalnoSigExpanded.txt 17SigNpalnoSigExpanded.txt NpalSigExpanded.txt | sort | uniq > Prosopis_Expanded.txt
# 411/821 genes for 120/194 gene families

cat Prosopis_Expanded.txt | grep -f - N0.tsv | cut -f8 | sed 's/, /\n/g' | sort | grep -f - Allgene_nameterms.txt > Prosopis_Expanded_geneterms.txt
# Same as before but only expanding in Npal branch
cat PcinSigExpanded.txt NalbSigExpanded.txt | sort | uniq | grep -v -f - NpalSigExpanded.txt > Algarrobo_Only_Expanded.txt
# 188/381 genes from 60/92 gene families

# Get groups present in Prosopis s.l group including Npal but absent in the others
cat N0.tsv | sed '1d' | cut -f1,11 | grep -Pv "\t$" > Groups/N0_Stor.txt # Do this for each species
cat N0_Bvar.txt N0_Gmax.txt N0_Mtru.txt N0_Pper.txt N0_Stor.txt N0_Vvin.txt | cut -f1 | sort | uniq | grep -v -f - N0_Npal.txt > N0_OtherAbsNpalPres.txt # 605/1309 genes from 452/878 gene families

# Significant orthogroups with N pallida significant
cat ../Base_family_results.txt | awk '$2 < 0.05' | cut -f1 | grep -f - ../Base_branch_probabilities.tab | awk '$4 < 0.05' > Sig_at_p.05_Npal.05.txt
# Changes of significantly rapidly evolving Families including Npal
cat Sig_at_p.05_Npal.05.txt | cut -f1 | grep -f - ../Base_change.tab > Change_at_p.05_Npal.05.txt
# Counts of significantly rapidly evolving Families including Npal
cat Sig_at_p.05_Npal.05.txt | cut -f1 | grep -f - ../Base_count.tab > Count_at_p.05_Npal.05.txt

# Get OG and CAFE p values to paste to GFF
cat Change_at_p.05_Npal.05.txt | awk '$4 > 0' | cut -f1 | grep -f - N0.tsv | cut -f8 | sed 's/, /\n/g' | cut -d'_' -f1,2 | grep -f - 2024-6-02_Npal_rnd3.merged.sorted.blast.interpro.emapper.decorated.gff | grep -P "\tmRNA\t" > Expanded_genes.gff
cat Expanded_genes.gff | cut -f9 | cut -d';' -f1 | sed 's/ID=//' > Expanded_gene_names.txt
while read line; do a=$(echo $line | cut -d'_' -f1,2 | grep -f - N0.tsv | cut -f1); b=$(grep $a ../3-CAFE/results/Base_family_results.txt | cut -f2); c=$(grep $a ../3-CAFE/results/Base_branch_probabilities.tab | cut -f4); d=$(grep $a OG_onlyNPal.txt); printf "$line\t$a\t$b\t$c\t$d\n"; done < Expanded_gene_names.txt > Expanded_gene_info.txt

############################################## COMPARATIVE GENOMICS: GO ENRICHMENT ANALYSIS #######################
# Append all gene models with GO terms
cat Npal_longesttranscript.txt | grep -f - 2024-6-02_Npal_rnd3.merged.sorted.blast.interpro.emapper.decorated.gff | grep -P "\tmRNA\t" | cut -f9 | grep "em_GOs=GO" | cut -d';' -f1 | sed 's/ID=//' > Allgene_names.txt
cat Npal_longesttranscript.txt | grep -f - 2024-6-02_Npal_rnd3.merged.sorted.blast.interpro.emapper.decorated.gff | grep -P "\tmRNA\t" | cut -f9 | grep "em_GOs=GO" | sed 's/em_GOs=GO/\tGO/' | cut -f2 | sed 's/,GO/, GO/g' > Allgene_terms.txt
paste Allgene_names.txt Allgene_terms.txt > Allgene_nameterms.txt # 12,860 genes with GO terms
awk '$4 > 0' Change_at_p.05_Npal.05.txt | cut -f1 | grep -f - N0.tsv | cut -f8 | sed 's/, /\n/g' | sort | grep -f - Allgene_nameterms.txt > Expanded_geneterms.txt #188/387 genes from 60/94 expanded families with GO terms

cat ../3-CAFE/results/Base_branch_probabilities.tab | awk '$4 <= 0.05' | awk '$2 > 0.05' | awk '$3 > 0.05' | awk '$5 > 0.05' | awk '$6 > 0.05' | awk '$7 > 0.05' | awk '$8 > 0.05' | awk '$9 > 0.05' | awk '$10 > 0.05' | cut -f1 | grep -f - Change_at_p.05_Npal.05.txt | awk '$4 > 0' | cut -f1 | grep -f - N0.tsv | cut -f8 | sed 's/, /\n/g' | sort | uniq | grep -f - Allgene_nameterms.txt > Expanded_NpalOnly_geneterms.txt
#45 OG significant only for Npal, 20 expanding, 25 contracting
#20 OG expanding with 111 genes, from which 73 had GO terms (from 17/20 orthogroups)

# In R, make the analysis using fisher exact test and weight01 algorithm with 4-TopGO_analysis.R

############################################### SYNTENY ANALYSIS #######################
## Pairwise comparisons using MCScanX. In each case will generate a dotplot and depth histogram pdf
python3 -m jcvi.compara.catalog ortholog --dbtype=prot --cscore=0.9 PCG_Npal PCG_Npal --no_strip_names
python3 -m jcvi.compara.synteny depth --histogram PCG_Npal.PCG_Npal.lifted.anchors &> PCG_Npal.histo
python3 -m jcvi.compara.synteny screen --simple PCG_Npal.PCG_Npal.anchors PCG_Npal.PCG_Npal.anchors.new
# Do the same for other comparisons
## Bed for circos (also Chr bed (no scaffold) and GC in 10kb windows) Copy the PCG_Npal.PCG_Npal.anchors.simple for circos link data
cat 2024-6-02_Npal_rnd3.merged.sorted.blast.interpro.emapper.decorated.gff | grep -P "\tgene\t" | cut -f1,4,5,9 | sed 's/;/\t/' | cut -f1-4 | sed 's/ID=//' | awk '{print $1 "\t" $2-1 "\t" $3 "\t" $4}' > Gene.bed
cat ../../../../Repeat_Library/Repeatmodeler/Masking/05_full_out/ppa_v1.asm.full_mask.gff3 | cut -f1,4,5,9 | sed 's/Target=//' | sed 's/ /_/g' | awk '{print $1 "\t" $2-1 "\t" $3 "\t" $4}' > Repeats.bed

# Colour the different links between chromosomes by adding at the start of the lines of anchors.simple a code based on the gene name and chromosome
# Define colors based on viridis palette
nano Colors.txt # Copy the following
# Chr1,#fde725
# Chr2,#cde11d
# Chr3,#98d83e
# Chr4,#67cc5c
# Chr5,#40bd72
# Chr6,#25ac82
# Chr7,#1f998a
# Chr8,#24878e
# Chr9,#2b748e
# Chr10,#34618d
# Chr11,#34618d
# Chr12,#453581
# Chr13,#481c6e
# Chr14,#440154
sed -i 's/,/\t/' Colors.txt

# Replace starts with colors
while read line;
do
  DUM=$(echo $line | cut -d' ' -f1 | grep -f - PCG_Npal.bed | cut -f1 | grep -w -f - Colors.txt | cut -f2); \
  echo $line | sed 's/^/'"$DUM"'*/'; \
done < PCG_Npal.PCG_Vvin_canonical.anchors.simple >> PCG_Npal.PCG_Vvin_canonical.anchors.simple.col

while read line;
do
  DUM=$(echo $line | cut -d' ' -f1 | grep -f - PCG_Npal.bed | cut -f1 | grep -w -f - Colors.txt | cut -f2); \
  echo $line | sed 's/^/'"$DUM"'*/'; \
done < PCG_Npal.PCG_Gmax_canonical.anchors.simple >> PCG_Npal.PCG_Gmax_canonical.anchors.simple.col

# Create seqids and layout files, chromosomes were rearranged based on visual inspection
nano seqids_Vvin.Npal.Gmax.col.txt
# Copy the following without comment
# 7,8,10,13,14,16,17,18,4,5,12,19,6,2,3,1,9,11,15
# Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12,Chr13,Chr14
# 3,10,19,16,4,6,12,14,2,11,13,8,5,20,17,1,15,18,9,7

nano layout_Vvin.Npal.Gmax.col.txt
# Copy the following without comment except for header
# # y, xstart, xend, rotation, color, label, va,  bed
#  .9,     .1,    .8,       0,      , Vitis vinifera, top, PCG_Vvin.bed
#  .5,     .1,    .8,       0,      , Neltuma pallida, top, PCG_Npal.bed
#  .1,     .1,    .8,       0,      , Glycine max, bottom, PCG_Gmax.bed
# # edges
# e, 0, 1, PCG_Npal.PCG_Vvin_canonical.anchors.simple.col
# e, 1, 2, PCG_Npal.PCG_Gmax_canonical.anchors.simple.col

# Draw the karyotype
python3 -m jcvi.graphics.karyotype seqids_Vvin.Npal.Gmax.col.txt layout_Vvin.Npal.Gmax.col.txt --figsize=32x9 -o Karyotype_Vvin.Npal.Gmax.col.pdf --nocircles

############################################### GENOME COUNTS #######################
## Assembly information was obtained from the NCBI/Ensemble webpage
# Gene characteristics using genestats script (similar for all species)
../../genestats.sh Algarrobo_rnd3_chr.gff > Algarrobo_rnd3_chr.stats &
../../genestats.sh Algarrobo_rnd3_scf.gff > Algarrobo_rnd3_scf.stats &
stats=Algarrobo_rnd3_chr.stats
# Average gene length
cat $stats | grep -P "\tgene\t" | awk '{print $3-$2}' | awk '{s+=$0} END {print s/NR}'
# Average transcript length
cat $stats | cut -f2 | awk '{s+=$0} END {print s/NR}'
# Average CDS length
cat $stats | cut -f8 | awk '{s+=$0} END {print s/NR}'
# Average number of exons
cat $stats | cut -f3 | awk '{s+=$0} END {print s/NR}'
# Average exon length
cat $stats | awk '$3>0' | awk '{print $4/$3}' | awk '{s+=$0} END {print s/NR}'
# Average intron length
cat $stats | awk '$5>0' |  awk '{print $6/$5}' | awk '{s+=$0} END {print s/NR}'
## BUSCO of all species
busco -i GCF_004799145.1.Nalb.fna -o Nalb_busco -m genome -c 20 --lineage embryophyta_odb10
busco -i GCA_014851425.1.Stor.fna -o Stor_busco -m genome -c 20 --lineage embryophyta_odb10
busco -i GCA_022379115.2.Bvar.fna -o Bvar_busco -m genome -c 20 --lineage embryophyta_odb10
busco -i Glycine_max.Glycine_max_v2.1.dna.toplevel.fa -o Gmax_busco -m genome -c 20 --lineage embryophyta_odb10
busco -i Medicago_truncatula.MedtrA17_4.0.dna.toplevel.fa -o Mtru_busco -m genome -c 20 --lineage embryophyta_odb10
busco -i Prunus_persica.Prunus_persica_NCBIv2.dna.toplevel.fa -o Pper_busco -m genome -c 20 --lineage embryophyta_odb10
busco -i Vitis_vinifera.PN40024.v4.dna.toplevel.fa -o Vvin_busco -m genome -c 20 --lineage embryophyta_odb10
busco -i ppa_v1.asm.fasta -o Npal_busco -m genome -c 20 --lineage embryophyta_odb10
busco -i GCA_029017545.1.Pcin.fna -o Pcin_busco -m genome -c 20 --lineage embryophyta_odb10
#############################################
