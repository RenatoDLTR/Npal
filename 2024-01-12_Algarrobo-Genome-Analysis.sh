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
# Similar for S. tor:
cat 1-Genomes/Stor/Senna_tora.gff | grep -v "#" | grep "protein_coding" | cut -f9 | sed 's/Name=/\t/' | cut -f2 | cut -d';' -f1 | grep -f - 1-Genomes/Stor/Senna_tora.gff > PCG_Stor.gff
python3 -m jcvi.formats.gff bed --type=mRNA --primary_only PCG_Stor.gff -o PCG_Stor.bed
cat PCG_Stor.bed | cut -f4 | grep -f - PCG_Stor.gff | grep -P "\tCDS" | cut -f9 | sed 's/protein_id=/\t/' | cut -f2 | cut -d';' -f1 | sort | uniq | seqkit grep -f - 1-Sequences/Senna_tora.GCA_014851425.1.pep.fa > PCG_Stor.pep
# For N. alb:
cat 1-Genomes/Nalb/genomic.gff | grep -v "#" | grep "protein_coding" | cut -f9 | sed 's/Name=/\t/' | cut -f2 | cut -d';' -f1 | grep -f - 1-Genomes/Nalb/genomic.gff > PCG_Nalb.gff
# Generate bed file
python3 -m jcvi.formats.gff bed --type=mRNA --primary_only PCG_Nalb.gff -o PCG_Nalb.bed
# Catch peptide sequences with names
cat PCG_Nalb.bed | cut -f4 | grep -f - PCG_Nalb.gff | grep -P "\tCDS" | cut -f9 | sed 's/protein_id=/\t/' | cut -f2 | cut -d';' -f1 | sort | uniq | seqkit grep -f - 1-Sequences/Prosopis_alba.GCF_004799145.1.pep.fa > PCG_Nalb.pep
# Correct the names in bed file
while read line; do grep $line PCG_Nalb.gff | grep -P "\tCDS" | cut -f9 | cut -d ';' -f1 | sed 's/ID=cds-//' | uniq; done < <(cat PCG_Nalb.bed | cut -f4) > PCG_Nalb.names
paste PCG_Nalb.names PCG_Nalb.bed | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $6 "\t" $7}' > PCG_Nalb_2.bed
mv PCG_Nalb.bed PCG_Nalbprev.bed
mv PCG_Nalb_2.bed PCG_Nalb.bed
# Similar for P cin:
cat 1-Genomes/Pcin/genomic.gff | grep -v "#" | grep "protein_coding" | cut -f9 | sed 's/Name=/\t/' | cut -f2 | cut -d';' -f1 | grep -f - 1-Genomes/Pcin/genomic.gff > PCG_Pcin.gff
python3 -m jcvi.formats.gff bed --type=mRNA --primary_only PCG_Pcin.gff -o PCG_Pcin.bed
cat PCG_Pcin.bed | cut -f4 | grep -f - PCG_Pcin.gff | grep -P "\tCDS" | cut -f9 | sed 's/protein_id=/\t/' | cut -f2 | cut -d';' -f1 | sort | uniq | seqkit grep -f - 1-Sequences/Prosopis_cineraria.GCF_029017545.1.pep.fa > PCG_Pcin.pep
while read line; do grep $line PCG_Pcin.gff | grep -P "\tCDS" | cut -f9 | cut -d ';' -f1 | sed 's/ID=cds-//' | uniq; done < <(cat PCG_Pcin.bed | cut -f4) > PCG_Pcin.names &
paste PCG_Pcin.names PCG_Pcin.bed | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $6 "\t" $7}' > PCG_Pcin_2.bed
mv PCG_Pcin.bed PCG_Pcinprev.bed
mv PCG_Pcin_2.bed PCG_Pcin.bed
# For Ensembl data (G. max, M. tru, P. per, V. vin).
# Get only protein coding genes
cat 1-Genomes/Gmax/Glycine_max.Glycine_max_v2.1.58.gff3 | grep -v "#" | grep "protein_coding" | cut -f9 | sed 's/Name=/\t/' | cut -f2 | cut -d';' -f1 | grep -f - 1-Genomes/Gmax/Glycine_max.Glycine_max_v2.1.58.gff3 > PCG_Gmax.gff
# Catch canonical transcripts (representatives)
cat PCG_Gmax.gff | grep -v "coding;transcript" > PCG_Gmax_canonical.gff
# Generate bed file with primary transcripts
python3 -m jcvi.formats.gff bed --type=mRNA --primary_only PCG_Gmax_canonical.gff -o PCG_Gmax_canonical.bed
# Catch peptide sequences with matching names in bed file
cat PCG_Gmax_canonical.bed | cut -f4 | sed 's/transcript://' | seqkit grep -f - 1-Sequences/Glycine_max.Glycine_max_v2.1.pep.all.fa > PCG_Gmax_canonical.pep
# Similar for M. tru:
cat 1-Genomes/Mtru/Medicago_truncatula.MedtrA17_4.0.58.gff3 | grep -v "#" | grep "protein_coding" | cut -f9 | sed 's/Name=/\t/' | cut -f2 | cut -d';' -f1 | grep -f - 1-Genomes/Mtru/Medicago_truncatula.MedtrA17_4.0.58.gff3 > PCG_Mtru.gff
cat PCG_Mtru.gff | grep -v "coding;transcript" > PCG_Mtru_canonical.gff
python3 -m jcvi.formats.gff bed --type=mRNA --primary_only PCG_Mtru_canonical.gff -o PCG_Mtru_canonical.bed
cat PCG_Mtru_canonical.bed | cut -f4 | sed 's/transcript://' | seqkit grep -f - 1-Sequences/Medicago_truncatula.MedtrA17_4.0.pep.all.fa > PCG_Mtru_canonical.pep
sed -i 's/transcript://' PCG_Mtru_canonical.bed
# Similar for P. per:
cat 1-Genomes/Pper/Prunus_persica.Prunus_persica_NCBIv2.58.gff3 | grep -v "#" | grep "protein_coding" | cut -f9 | sed 's/Name=/\t/' | cut -f2 | cut -d';' -f1 | grep -f - 1-Genomes/Pper/Prunus_persica.Prunus_persica_NCBIv2.58.gff3 > PCG_Pper.gff
cat PCG_Pper.gff | grep -v "coding;transcript" > PCG_Pper_canonical.gff
python3 -m jcvi.formats.gff bed --type=mRNA --primary_only PCG_Pper_canonical.gff -o PCG_Pper_canonical.bed
cat PCG_Pper_canonical.bed | cut -f4 | sed 's/transcript://' | seqkit grep -f - 1-Sequences/Prunus_persica.Prunus_persica_NCBIv2.pep.all.fa > PCG_Pper_canonical.pep
sed -i 's/transcript://' PCG_Pper_canonical.bed
# And for V. vin:
cat 1-Genomes/Vvin/Vitis_vinifera.PN40024.v4.58.gff3 | grep -v "#" | grep "protein_coding" | cut -f9 | sed 's/Name=/\t/' | cut -f2 | cut -d';' -f1 | grep -f - 1-Genomes/Vvin/Vitis_vinifera.PN40024.v4.58.gff3 > PCG_Vvin.gff
cat PCG_Vvin.gff | grep -v "coding;transcript" > PCG_Vvin_canonical.gff
python3 -m jcvi.formats.gff bed --type=mRNA --primary_only PCG_Vvin_canonical.gff -o PCG_Vvin_canonical.bed
cat PCG_Vvin_canonical.bed | cut -f4 | sed 's/transcript://' | cut -d'_' -f1 | seqkit grep -f - 1-Sequences/Vitis_vinifera.PN40024.v4.pep.all.fa > PCG_Vvin_canonical.pep
sed -i 's/transcript://' PCG_Vvin_canonical.bed
sed -i 's/_t/_P/' PCG_Vvin_canonical.bed
#For Neltuma pallida
cat ../1-Genomes/Npal/NNew/2024_Algarrobo_rnd6.PCgenes.aed1.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | cut -d'_' -f1,2 | grep -f - ../1-Genomes/Npal/NNew/2024_Algarrobo_rnd6.PCgenes.aed1.ipr.blast.emapper.decorated.gff > PCG_Npal.gff
python3 -m jcvi.formats.gff bed --type=mRNA --primary_only PCG_Npal.gff -o PCG_Npal.bed
cat PCG_Npal.bed | cut -f4 | seqkit grep -f - ../1-Genomes/Npal/NNew/2024_Algarrobo_rnd6.proteins.aed1.blast.fasta > PCG_Npal.pep

# Copy all peptide sequences file to a new directory
mkdir NewPeps
cp *.pep Peps/
# Change pep files names to full species_names.fa for subsequent analyses

############################################### COMPARATIVE GENOMICS: ORTHOLOGY INFERENCE #######################
# Copy the tree file RAxML_bipartitions.RaxMLTree. This was generated using single-copy genes in UNMSM lab with the PROTGAMMAGTR substitution model using single-copy genes from a previous orthofinder run with the same species
# Group orthologs using OrthoFinder
orthofinder -t 60 -a 60 -s RAxML_bipartitions.RaxMLTree -f NewPeps/ #"-s" will take into account the tree for the orthology inference
# Apart from the statistics table (which is different from the HOG statistics used for CAFE analysis), upsetR was used to give a graphical output of shared and unique genes. Script is 3-Orthogroup_analysis.R

############################################### COMPARATIVE GENOMICS: GENE FAMILY DYNAMICS #######################
# For CAFE analysis, two inputs are needed: the rooted ultrametric tree and the filtered gene counts. For the latter, use 3-Orthogroup_analysis.R
# For making the ultrametric tree, root node age was selected based on literature
make_ultrametric.py -r 125 RAxML_bipartitions.RaxMLTree
# Run CAFE analysis
cafe5 -c 50 -i hog_gene_counts.tsv -t RAxML_bipartitions.RaxMLTree.ultrametric.tre
# Summarize and plot results
cafeplotter -i results/ -o cafeplotter_results --ignore_branch_length
# Count families significant at 0.05
cat ../Base_family_results.txt | awk '$2 < 0.05 {print $0}' > Sig_at_p.05.txt # 1814
# Count families with Npal p value branch significant at 0.05
cut -f1 Sig_at_p.05.txt | grep -f - ../Base_branch_probabilities.tab | awk '$5 < 0.05' > Sig_at_p.05_Npal.05.txt # 937
# Obtain changes of significantly rapidly evolving Families including Npal
cat Sig_at_p.05_Npal.05.txt | cut -f1 | grep -f - ../Base_change.tab > Change_at_p.05_Npal.05.txt
# Obtain counts of significantly rapidly evolving Families including Npal
cat Sig_at_p.05_Npal.05.txt | cut -f1 | grep -f - ../Base_count.tab > Count_at_p.05_Npal.05.txt
# Count estimated number of genes/orthogroups for each species and nodes from the cafeplotter summary. Also, from CAFE tables got number of (significant) expanded/contracted families for each species and nodes.
# Get family expansion rate for all species and each species, as well as the number of expanded/contracted/non-changing families from the CAFE report.

############################################### COMPARATIVE GENOMICS: GO ENRICHMENT ANALYSIS #######################
# Append all gene models with GO terms
cat 2024_Algarrobo_rnd6.PCgenes.aed1.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | grep "em_GOs" | cut -f9 | sed 's/;Parent/\t/' | sed 's/em_GOs=/\t/' | cut -f1,3 | grep -Pv "\t$" | sed 's/;em/\t/' | cut -f1,2 | grep -Pv "\t$" > All_gene_models_GO.names
# Get names of all genes
cat 2024_Algarrobo_rnd6.PCgenes.aed1.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" | cut -f9 | cut -d'_' -f1,2 | sort | uniq | sed 's/$/_R1;/' > AllGenes.txt
# Get names of expanded orthogroups
cat Change_at_p.05_Npal.05.txt | awk '$5 > 0' | cut -f1 > Orthogroups_expanded.names
# Get names of genes in expanded orthogroups
while read line; do grep $line N0.tsv | cut -f8 | sed 's/, /\n/g'; done < Orthogroups_expanded.names | sort > Genes_expanded.names
# Get list of genes appended to expanded orthogroups
while read line; do grep $line Results_Mar22/Phylogenetic_Hierarchical_Orthogroups/N0.tsv | cut -f1,8; done < Orthogroups_expanded.names | sort > Genes_expanded_listed.names
# Get list of GO terms for each gene (represented by one gene model)
while read line; do grep $line All_gene_models_GO.names | head -n1; done < <(cat AllGenes.txt | cut -d'_' -f1,2 | sed 's/ID=//') > GOterms.names
# Remove gene model suffix to match gene names
cat GOterms.names | sed 's/ID=//' | sed 's/_R/\t/' | cut -f1,3 > GOterms_genes.names
# Modify gene model names to match gene names
cat AllGenes.txt | cut -d'_' -f1,2 | sed 's/ID=//' > AllGenes_mod.txt
# Modify expanded gene model names to match gene names
cat Genes_expanded.names | cut -d'_' -f1,2 > Genes_expanded_mod.names
# In R, make the analysis using fisher exact test and weight01 algorithm with 4-TopGO_analysis.R
# Do the same analysis for genes from orthogroups expanded only in Neltuma pallida
while read line; do grep $line Results_Mar22/Phylogenetic_Hierarchical_Orthogroups/N0.tsv | cut -f8 | sed 's/, /\n/g'; done < Orthogroup_expanded_NpalOnly.names | sort > Genes_expanded_NpalOnly.names
cat Genes_expanded_NpalOnly.names | cut -d'_' -f1,2 > Genes_expanded_NpalOnly_mod.names

############################################### EXTRA CAFE ANALYSIS #######################
# Append gene names and orthogroup names to GO terms
cat Genes_expanded_mod.names | grep -f - GOterms_genes.names > GOterms_expanded_mod.names
# All GO terms obtained from topGO full table, copied to a new file with nano.
while read line;
do
  if grep $line GOterms_expanded_mod.names >/dev/null
  then grep $line GOterms_expanded_mod.names | cut -f1 | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/, /g' | sed 's/^/'"$line"'\t/'
else echo $line
fi; done < AllGOterms.txt >> GOterms_expanded_genesinvolved.names
# Get list of genes for all expanded orthogroups
cat N0_Npal.tsv | grep -f Orthogroups_expanded.names > Expanded_N0.tsv
# Append expanded orthogroup names to each GO term
while read line;
do
  if echo $line | sed 's/, /\n/g' | grep -f - Expanded_N0.tsv >/dev/null
  then echo $line | sed 's/, /\n/g' | grep -f - Expanded_N0.tsv | cut -f1 | sort | uniq | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/, /g'
else echo $line
fi; done < <(cut -f2 GOterms_expanded_genesinvolved.names) >> GOterms_expanded_orthogroupsinvolved.names
# Append genes expanded only in Neltuma pallida to each GO term
cat Genes_expanded_NpalOnly_mod.names | grep -f - GOterms_genes.names > GOterms_expanded_NpalOnly.names
while read line;
do
  if grep $line GOterms_expanded_NpalOnly.names >/dev/null
  then grep $line GOterms_expanded_NpalOnly.names | cut -f1 | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/, /g' | sed 's/^/'"$line"'\t/'
else echo $line
fi; done < AllGOterms.txt >> GOterms_expanded_genesinvolved_NpalOnly.names
# Get orthogroup names of genes expanded only in Neltuma pallida
cat N0_Npal.tsv | grep -f Orthogroups_expanded_NpalOnly.names > Expanded_N0_NpalOnly.tsv
while read line;
do
  if echo $line | sed 's/, /\n/g' | grep -f - Expanded_N0_NpalOnly.tsv >/dev/null
  then echo $line | sed 's/, /\n/g' | grep -f - Expanded_N0_NpalOnly.tsv | cut -f1 | sort | uniq | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/, /g'
else echo $line
fi; done < <(cut -f2 GOterms_expanded_genesinvolved_NpalOnly.names) >> GOterms_expanded_orthogroupsinvolved_NpalOnly.names
# Get full annotation for expanded genes
cat Genes_expanded_mod.names | sed 's/NPAL/ID=NPAL/' | grep -f - 2024_Algarrobo_rnd6.PCgenes.aed1.ipr.blast.emapper.decorated.gff | grep -P "\tmRNA" > Gene_models_expanded.gff
# Get names of gene models in expanded orthogroups (just gene name with -R1 suffix)
while read line; do grep $line Expanded_N0.tsv | cut -f1; done < <(cat Gene_models_expanded.gff | cut -f9 | cut -d'_' -f1,2 | sed 's/ID=//') > Orthogroup_gene_models.names
# Append family p values to expanded orthogroups
while read line; do grep $line Base_family_results.txt | cut -f2; done < Orthogroup_gene_models.names > Orthogroup_familyp.names
# Append Neltuma pallida p values to each gene of each orthogroup
while read line; do grep $line Base_branch_probabilities.tab | cut -f2; done < Orthogroup_gene_models.names > Orthogroup_branch.names
# Append a "Y" to each gene of orthogroups if they are expanded only in Neltuma pallida, "N" otherwise
while read line;
do
  if grep $line Orthogroup_expanded_NpalOnly.names >/dev/null
  then echo "Y"
else echo "N"
fi; done < Orthogroup_Genes.names > Orthogroup_OnlyNpal.names

############################################### SYNTENY ANALYSIS #######################
## Pairwise comparisons using MCScanX. In each case will generate a dotplot and depth histogram pdf
# Compare Npal genome with itself.
python3 -m jcvi.compara.catalog ortholog --dbtype=prot --cscore=0.9 PCG_Npal PCG_Npal --no_strip_names
python3 -m jcvi.compara.synteny depth --histogram PCG_Npal.PCG_Npal.lifted.anchors &> PCG_Npal.histo
#Compare Npal with Gmax
python3 -m jcvi.compara.catalog ortholog --dbtype=prot --cscore=0.9 PCG_Npal PCG_Gmax_canonical --no_strip_names
python3 -m jcvi.compara.synteny depth --histogram PCG_Npal.PCG_Gmax_canonical.lifted.anchors &> PCG_Npal_PCG_Gmax_canonical.histo
#Compare Npal with Mtru
python3 -m jcvi.compara.catalog ortholog --dbtype=prot --cscore=0.9 PCG_Npal PCG_Mtru_canonical --no_strip_names
python3 -m jcvi.compara.synteny depth --histogram PCG_Npal.PCG_Mtru_canonical.lifted.anchors &> PCG_Npal_PCG_Mtru_canonical.histo
#Compare Npal with Vvin
python3 -m jcvi.compara.catalog ortholog --dbtype=prot --cscore=0.9 PCG_Npal PCG_Vvin_canonical --no_strip_names
python3 -m jcvi.compara.synteny depth --histogram PCG_Npal.PCG_Vvin_canonical.lifted.anchors &> PCG_Npal_PCG_Vvin_canonical.histo
#Compare Npal with Nalb
python3 -m jcvi.compara.catalog ortholog --dbtype=prot --cscore=0.9 PCG_Npal PCG_Nalb --no_strip_names
python3 -m jcvi.compara.synteny depth --histogram PCG_Npal.PCG_Nalb.lifted.anchors &> PCG_Npal_PCG_Nalb.histo
#Compare Npal with Pcin
python3 -m jcvi.compara.catalog ortholog --dbtype=prot --cscore=0.9 PCG_Npal PCG_Pcin --no_strip_names
python3 -m jcvi.compara.synteny depth --histogram PCG_Npal.PCG_Pcin.lifted.anchors &> PCG_Npal_PCG_Pcin.histo
# Create comparative maps for Vvin - Npal and Gmax - Npal
python3 -m jcvi.compara.synteny screen --simple PCG_Npal.PCG_Vvin_canonical.anchors PCG_Npal.PCG_Vvin_canonical.anchors.new
python3 -m jcvi.compara.synteny screen --simple PCG_Npal.PCG_Gmax_canonical.anchors PCG_Npal.PCG_Gmax_canonical.anchors.new
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
  DUM=$(echo $line | cut -d' ' -f1 | grep -f - PCG_Npal.bed | cut -f1 | grep -w -f - Colors.txt | cut -f2)
  echo $line | sed 's/^/'"$DUM"'*/';
done < PCG_Npal.PCG_Vvin_canonical.anchors.simple >> PCG_Npal.PCG_Vvin_canonical.anchors.simple.col
while read line;
do
  DUM=$(echo $line | cut -d' ' -f1 | grep -f - PCG_Npal.bed | cut -f1 | grep -w -f - Colors.txt | cut -f2)
  echo $line | sed 's/^/'"$DUM"'*/';
done < PCG_Npal.PCG_Gmax_canonical.anchors.simple >> PCG_Npal.PCG_Gmax_canonical.anchors.simple.col
# Create seqids and layout files, chromosomes were rearranged based on visual inspection
nano seqids_Vvin.Npal.Gmax.col.txt
# Copy the following without comment
## 7,8,10,13,14,16,17,18,4,5,12,19,6,2,3,1,9,11,15
## Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12,Chr13,Chr14
## 3,10,19,16,4,6,12,14,2,11,13,8,5,20,17,1,15,18,9,7
nano layout_Vvin.Npal.Gmax.col.txt
# Copy the following without comment except for header
## # y, xstart, xend, rotation, color, label, va,  bed
##  .9,     .1,    .8,       0,      , Vitis vinifera, top, PCG_Vvin.bed
##  .5,     .1,    .8,       0,      , Neltuma pallida, top, PCG_Npal.bed
##  .1,     .1,    .8,       0,      , Glycine max, bottom, PCG_Gmax.bed
## # edges
## e, 0, 1, PCG_Npal.PCG_Vvin_canonical.anchors.simple.col
## e, 1, 2, PCG_Npal.PCG_Gmax_canonical.anchors.simple.col
# Draw the karyotype
python3 -m jcvi.graphics.karyotype seqids_Vvin.Npal.Gmax.col.txt layout_Vvin.Npal.Gmax.col.txt --figsize=32x9 -o Karyotype_Vvin.Npal.Gmax.col.pdf --nocircles

############################################### GENOME ANALYSIS #######################
## Assembly information was obtained from the NCBI/Ensemble webpage
# Gene characteristics using genestats script (similar for all species)
genestats.sh PCG_Npal.gff > PCG_Npal.stats &
# Average gene length
cat PCG_Npal.bed | awk '{print $3-$2}' | awk '{s+=$0} END {print s/NR}'
# Average transcript length
cat PCG_Npal.stats | cut -f2 | awk '{s+=$0} END {print s/NR}'
# Average CDS length
cat PCG_Npal.stats | cut -f8 | awk '{s+=$0} END {print s/NR}'
# Average number of exons
cat PCG_Npal.stats | cut -f3 | awk '{s+=$0} END {print s/NR}'
# Average exon length
cat PCG_Npal.stats | awk '$3>0' | awk '{print $4/$3}' | awk '{s+=$0} END {print s/NR}'
# Average intron length
cat PCG_Npal.stats | awk '$5>0' |  awk '{print $6/$5}' | awk '{s+=$0} END {print s/NR}'
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
