### Packages
library(topGO)
### Working dir
setwd("R:/2-NTNU-Holomuseomics/SharedAstbury/2023-12_NuclearGenome/")

geneID2GO <- readMappings(file="Allgene_nameterms.txt")
Expanded <- read.table(file="Expanded_geneterms.txt", sep = "\t")
Expanded <- as.character(Expanded[,1])
GO2geneID <- inverseList(geneID2GO)
geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% Expanded))
names(geneList) <- geneNames

# BP terms
BPGOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
resultBPFisher <- runTest(BPGOdata, algorithm = "weight01", statistic = "fisher")
resultBPFisher
# Description:
# Ontology: BP
# 'weight01' algorithm with the 'fisher' test
# 2758 GO terms scored: 45 terms with p < 0.01
# Annotation data:
#     Annotated genes: 11133
#     Significant genes: 164
#     Min. no. of genes annotated to a GO: 10
#     Nontrivial nodes: 602

# Put the number of GO terms in topNodes
allBPRes <- GenTable(BPGOdata, classicFisher = resultBPFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 2758)
write.table(allBPRes, file = "BP_Table_Expanded.txt", sep = "\t", quote = FALSE)
showSigOfNodes(BPGOdata, score(resultBPFisher), firstSigNodes = 10, useInfo = "all")
printGraph(BPGOdata, resultBPFisher, firstSigNodes = 10, fn.prefix = "BP_Graph_Expanded", useInfo = "all", pdfSW = TRUE)
### MF terms
MFGOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
resultMFFisher <- runTest(MFGOdata, algorithm = "weight01", statistic = "fisher")
resultMFFisher
# Description:
# Ontology: MF
# 'weight01' algorithm with the 'fisher' test
# 724 GO terms scored: 25 terms with p < 0.01
# Annotation data:
#     Annotated genes: 9599
#     Significant genes: 150
#     Min. no. of genes annotated to a GO: 10
#     Nontrivial nodes: 192

# Put the number of GO terms in topNodes
allMFRes <- GenTable(MFGOdata, classicFisher = resultMFFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 724)
write.table(allMFRes, file = "MF_Table_Expanded.txt", sep = "\t", quote = FALSE)
showSigOfNodes(MFGOdata, score(resultMFFisher), firstSigNodes = 10, useInfo = "all")
printGraph(MFGOdata, resultMFFisher, firstSigNodes = 10, fn.prefix = "MF_Graph_Expanded", useInfo = "all", pdfSW = TRUE)

### CC terms
CCGOdata <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
resultCCFisher <- runTest(CCGOdata, algorithm = "weight01", statistic = "fisher")
resultCCFisher
# Description:
# Ontology: CC
# 'weight01' algorithm with the 'fisher' test
# 401 GO terms scored: 9 terms with p < 0.01
# Annotation data:
#     Annotated genes: 10953
#     Significant genes: 131
#     Min. no. of genes annotated to a GO: 10
#     Nontrivial nodes: 103
# Put the number of GO terms in topNodes
allCCRes <- GenTable(CCGOdata, classicFisher = resultCCFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 401)
write.table(allCCRes, file = "CC_Table_Expanded.txt", sep = "\t", quote = FALSE)
showSigOfNodes(CCGOdata, score(resultCCFisher), firstSigNodes = 10, useInfo = "all")
printGraph(CCGOdata, resultCCFisher, firstSigNodes = 10, fn.prefix = "CC_Graph_Expanded", useInfo = "all", pdfSW = TRUE)

#### Expanded only in Neltuma pallida
ExpandedNO <- read.table(file="Expanded_NpalOnly_geneterms.txt", sep = "\t")
ExpandedNO <- as.character(ExpandedNO[,1])
geneID2GONO <- readMappings(file = "Allgene_nameterms.txt")
GO2geneIDNO <- inverseList(geneID2GONO)
geneNamesNO <- names(geneID2GONO)
geneListNO <- factor(as.integer(geneNamesNO %in% ExpandedNO))
names(geneListNO) <- geneNamesNO
# BP terms
BPGOdataNO <- new("topGOdata", ontology = "BP", allGenes = geneListNO, annot = annFUN.gene2GO, gene2GO = geneID2GONO, nodeSize = 10)
resultBPFisherNO <- runTest(BPGOdataNO, algorithm = "weight01", statistic = "fisher")
resultBPFisherNO
# Description:
# Ontology: BP
# 'weight01' algorithm with the 'fisher' test
# 2758 GO terms scored: 38 terms with p < 0.01
# Annotation data:
#     Annotated genes: 11133
#     Significant genes: 65
#     Min. no. of genes annotated to a GO: 10
#     Nontrivial nodes: 304
# Put the number of GO terms in topNodes
allBPResNO <- GenTable(BPGOdataNO, classicFisher = resultBPFisherNO, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 2758)
write.table(allBPResNO, file = "BP_Table_Expanded_Npal_Only.txt", sep = "\t", quote = FALSE)
showSigOfNodes(BPGOdataNO, score(resultBPFisherNO), firstSigNodes = 10, useInfo = "all")
printGraph(BPGOdataNO, resultBPFisherNO, firstSigNodes = 10, fn.prefix = "BP_Graph_Expanded_Npal_Only", useInfo = "all", pdfSW = TRUE)
### MF terms
MFGOdataNO <- new("topGOdata", ontology = "MF", allGenes = geneListNO, annot = annFUN.gene2GO, gene2GO = geneID2GONO, nodeSize = 10)
resultMFFisherNO <- runTest(MFGOdataNO, algorithm = "weight01", statistic = "fisher")
resultMFFisherNO
# Description:
# Ontology: MF
# 'weight01' algorithm with the 'fisher' test
# 724 GO terms scored: 16 terms with p < 0.01
# Annotation data:
#     Annotated genes: 9599
#     Significant genes: 54
#     Min. no. of genes annotated to a GO: 10
#     Nontrivial nodes: 80
# Put the number of GO terms in topNodes
allMFResNO <- GenTable(MFGOdataNO, classicFisher = resultMFFisherNO, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 724)
write.table(allMFResNO, file = "MF_Table_Expanded_Npal_Only.txt", sep = "\t", quote = FALSE)
showSigOfNodes(MFGOdataNO, score(resultMFFisherNO), firstSigNodes = 10, useInfo = "all")
printGraph(MFGOdataNO, resultMFFisherNO, firstSigNodes = 10, fn.prefix = "MF_Graph_Expanded_Npal_Only", useInfo = "all", pdfSW = TRUE)
### CC terms
CCGOdataNO <- new("topGOdata", ontology = "CC", allGenes = geneListNO, annot = annFUN.gene2GO, gene2GO = geneID2GONO, nodeSize = 10)
resultCCFisherNO <- runTest(CCGOdataNO, algorithm = "weight01", statistic = "fisher")
resultCCFisherNO
# Description:
# Ontology: CC
# 'weight01' algorithm with the 'fisher' test
# 401 GO terms scored: 9 terms with p < 0.01
# Annotation data:
#     Annotated genes: 10953
#     Significant genes: 54
#     Min. no. of genes annotated to a GO: 10
#     Nontrivial nodes: 66
# Put the number of GO terms in topNodes
allCCResNO <- GenTable(CCGOdataNO, classicFisher = resultCCFisherNO, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 401)
write.table(allCCResNO, file = "CC_Table_Expanded_Npal_Only.txt", sep = "\t", quote = FALSE)
showSigOfNodes(CCGOdataNO, score(resultCCFisherNO), firstSigNodes = 10, useInfo = "all")
printGraph(CCGOdataNO, resultCCFisherNO, firstSigNodes = 10, fn.prefix = "CC_Graph_Expanded_Npal_Only", useInfo = "all", pdfSW = TRUE)
