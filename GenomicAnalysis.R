install.packages("BioManager")
install.packages("NMF")
install.packages("pheatmap")
BiocManager::install("GenomicDataCommons")
BiocManager::install("maftools")
BiocManager::install("TCGAbiolinks")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("mclust")

library(GenomicDataCommons)
library(maftools)
library(TCGAbiolinks)
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
library('NMF')
library('pheatmap')
library(mclust)



#Recover LUAD maf file
luad.maf<-GDCquery_Maf(tumor="LUAD", save.csv = FALSE, directory = "LUAD", pipelines = "mutect2")
sort(colnames(maf))
luad.maf <- read.maf(maf = luad.maf)

#Recover LUSC maf file 
lusc.maf<-GDCquery_Maf(tumor="LUSC", save.csv = FALSE, directory = "LUSC", pipelines = "mutect2")
sort(colnames(lusc.maf))
lusc.maf <- read.maf(maf = lusc.maf)

#Shows sample summary.
getSampleSummary(luad.maf)
#Shows gene summary.
getGeneSummary(luad.maf)
#shows clinical data associated with samples
getClinicalData(luad.maf)
#Shows all fields in MAF
getFields(luad.maf)
#Writes maf summary to an output file with basename luad.
write.mafSummary(maf = luad.maf, basename = 'luad')


#LUAD
# Supplementary Figure 1
plotmafSummary(maf = luad.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

# Gene selection
#oncoplot for frequently mutated genes.(Looking at those who were frequently mutated in more than 15% of the cohort)
oncoplot(maf = luad.maf, top = 35)

# Supplementary Figure 3
OncogenicPathways(maf = luad.maf)

#LUSC
#Shows sample summary.
getSampleSummary(lusc.maf)
#Shows gene summary.
getGeneSummary(lusc.maf)
#shows clinical data associated with samples
getClinicalData(lusc.maf)
#Shows all fields in MAF
getFields(lusc.maf)
#Writes maf summary to an output file with basename lusc.
write.mafSummary(maf = lusc.maf, basename = 'lusc')

# Supplementary Figure 2
plotmafSummary(maf = lusc.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

# Gene selection
#oncoplot for frequently mutated genes.(Looking at those who were frequently mutated in more than 15% of the cohort)
oncoplot(maf = lusc.maf, top = 35)

# Supplementary Figure 4
OncogenicPathways(maf = lusc.maf)

# Supplementary Figure 13 (to cohort comparison)
#Considering only genes which are mutated in at-least in 100 samples in one of the cohort to avoid bias due to genes mutated in single sample.
luad.vs.lusc <- mafCompare(m1 = luad.maf, m2 = lusc.maf, m1Name = 'LUAD', m2Name = 'LUSC', minMut = 100)
print(luad.vs.lusc)

forestPlot(mafCompareRes = luad.vs.lusc, pVal = 0.1, geneFontSize = 0.8)


