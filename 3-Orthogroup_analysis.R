### Packages
library(UpSetR)
library(stringr)
library(data.table)
library(tidyverse)
library(reshape2)
### Working dir
setwd("R:/2-NTNU-Holomuseomics/SharedAstbury/2023-12_NuclearGenome/2024_Repeats/")

orthogroups <- read.table("N0.tsv",
                          header = TRUE, sep = "\t")
orthogroups_total <- orthogroups
orthogroups <- orthogroups[,-c(2, 3)]
for(i in 1:nrow(orthogroups)){
  for(j in 2:ncol(orthogroups)){
    orthogroups[i,j] <- ifelse(orthogroups[i,j] == "", 0, str_count(orthogroups[i,j], ",") + 1)
  }
}
plot(sort(as.numeric(orthogroups$Neltuma_pallida)))
orthogroups[which(as.numeric(orthogroups$Neltuma_pallida) > 20),]
for(i in 1:nrow(orthogroups)){
  for(j in 2:(ncol(orthogroups))){
    orthogroups[i,j] <- ifelse(orthogroups[i,j] > 0, 1, 0)
  }
}
orthogroups$Bauhinia_variegata <- as.numeric(orthogroups$Bauhinia_variegata)
orthogroups$Glycine_max <- as.numeric(orthogroups$Glycine_max)
orthogroups$Medicago_truncatula <- as.numeric(orthogroups$Medicago_truncatula)
orthogroups$Neltuma_alba <- as.numeric(orthogroups$Neltuma_alba)
orthogroups$Neltuma_pallida <- as.numeric(orthogroups$Neltuma_pallida)
orthogroups$Prosopis_cineraria <- as.numeric(orthogroups$Prosopis_cineraria)
orthogroups$Prunus_persica <- as.numeric(orthogroups$Prunus_persica)
orthogroups$Senna_tora <- as.numeric(orthogroups$Senna_tora)
orthogroups$Vitis_vinifera <- as.numeric(orthogroups$Vitis_vinifera)
upset(orthogroups,keep.order=T,text.scale = 1.2,
      intersections = list(list(c("Neltuma_pallida")),
                           list(c("Neltuma_alba")),
                           list(c("Prosopis_cineraria")),
                           list(c("Senna_tora")),
                           list(c("Bauhinia_variegata")),
                           list(c("Glycine_max")),
                           list(c("Medicago_truncatula")),
                           list(c("Prunus_persica")),
                           list(c("Vitis_vinifera")),
                           list(c("Neltuma_pallida", "Neltuma_alba")),
                           list(c("Neltuma_pallida", "Prosopis_cineraria")),
                           list(c("Neltuma_alba", "Prosopis_cineraria")),
                           list(c("Neltuma_pallida", "Neltuma_alba", "Prosopis_cineraria")),
                           list(c("Prunus_persica", "Medicago_truncatula", "Neltuma_pallida", "Neltuma_alba", "Prosopis_cineraria", "Senna_tora", "Bauhinia_variegata", "Vitis_vinifera", "Glycine_max")),
                           list(c("Prunus_persica", "Medicago_truncatula", "Neltuma_alba", "Prosopis_cineraria", "Senna_tora", "Bauhinia_variegata", "Vitis_vinifera", "Glycine_max")),
                           list(c("Prunus_persica", "Medicago_truncatula", "Prosopis_cineraria", "Senna_tora", "Bauhinia_variegata", "Vitis_vinifera", "Glycine_max")),
                           list(c("Prunus_persica", "Medicago_truncatula", "Senna_tora", "Bauhinia_variegata", "Vitis_vinifera", "Glycine_max"))),
      main.bar.color = c(rep("black", 9), rep("darkgreen", 4), rep("darkblue", 3), "black")
                           )
## Filter for CAFE analysis
hog<- read.table("N0.tsv", header = TRUE, sep = "\t") #30,700
remove <- c()
for(i in 1:nrow(hog)){
  ifelse(length(which(hog[i,]==""))==8, remove <- c(remove, hog[i,1]), print(i))
}
#Exclude HOGs present in only 1 species
hog <- hog[-which((hog$HOG %in% remove)==TRUE),] #22,765
hog <- hog[,-c(2,3)]
rownames(hog) <- NULL
hog <- melt(hog, id.vars = 'HOG', variable.name = 'species', value.name = 'pid')
hog <- hog[-which(hog$pid==""),]
hog$n <- sapply(hog$pid, function(x) length(strsplit(x, ', ')[[1]]))
#Exclude HOGs with lots of genes in one or more species
#hog <- hog[-which(hog$n>100),]
remove2 <- hog[which(hog$n>100),1]
hog <- hog[-which((hog$HOG %in% remove2)==TRUE),] #22747 left
rownames(hog) <- NULL
counts <- dcast(hog, HOG ~ species, value.var = 'n', fill = 0)
counts$Desc <- 'n/a'
setcolorder(counts, 'Desc')
fwrite(counts, 'hog_gene_counts.tsv', sep = '\t')
