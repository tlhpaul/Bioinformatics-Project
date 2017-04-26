#Most of codes are from CIS/GCB 535 Bioinformatics in-class assignemnt.

library(GEOquery)
library(oligo)
library(limma)
phenoData <- read.AnnotatedDataFrame("phenotype_taiwan.csv", header=TRUE, sep=",")
phenoData
summary(phenoData)
pData <- read.table(file="phenotype_taiwan.csv", sep=",",header=T)
pData
celFileslist  <- list.celfiles("GSE19804",full.name=T)
celFiles <- celFileslist[]
print(celFiles)
affyRaw <- read.celfiles(celFiles)
genenorm <- rma(affyRaw)
genenorm
boxplot(genenorm)        
group <- factor(t(pData(phenoData)[,1]))
group
design <- model.matrix(~group,pData(phenoData))
design
fit1 <- lmFit(genenorm, design)
efit1 <- eBayes(fit1)
tt1 <- topTable(efit1, sort="P", n=500)
tt1[1:5,] 
allTT1 <- tt1[order(tt1$P),]
allTT1[1:100,]
with(allTT1, sum(adj.P.Val < 0.05))
library(biomaRt)
mymart=useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
mymart
head(tt1)
pidsTophits1 <- row.names(tt1)
pidsTophits1[1:5]
attributes <- listAttributes(mymart)
myannot1 <- getBM(attributes=c('affy_hg_u133_plus_2','chromosome_name','start_position','end_position', 'hgnc_symbol', 'external_gene_name','description'), filters = 'affy_hg_u133_plus_2', values = pidsTophits1, mart = mymart)
myannot1[1:50,]
write.csv(myannot1, file = "tt-taiwan.csv", quote = FALSE)
tt_ids1 <- cbind(tt1,pidsTophits1)
colnames(tt_ids1)[7] <- c("affyprobeid")
colnames(myannot1)[1] <- c("affyprobeid")
myannot1[1:50,]
tt_plus_annot1 <- merge(myannot1, tt_ids1, by="affyprobeid", all=TRUE)
tt_plus_annot1[1:50,]
tt_annot_s1 <- tt_plus_annot1[order(tt_plus_annot1$adj.P.Val),]
tt_annot_s1[1:50,]
write.csv(tt_annot_s1, file= "merged_taiwan.csv", quote = F)
tt_annot_s1[,5]
unique1(tt_annot_s1[,5])
uniqgenelist1 <- unique(tt_annot_s1[,5])
uniqgenelist1
length(uniqgenelist1)
write.table(uniqgenelist1, "genelist_taiwan.txt", quote = F, row.name = F, col.name = F)
