

######################################
## 1. FOLDER AND LIBRARY DEFINITION ##
######################################

rm(list = ls())
WK_DIR = "C:/Users/franc/Desktop/TESI/EDA2R/Analisi Revisione/Revisore 2/Rat/"
setwd(WK_DIR)
Sys.setenv(LANG = "en")

neededlibraries <- c("edgeR", "RColorBrewer", "forcats", "parallel", "doParallel", "BiocParallel", "DESeq2", "ggplot2", "ggrepel", "org.Hs.eg.db", "org.Mm.eg.db", "GSVA", "writexl",  
                     "qusage", "clipr", "readtext", "scales", "sva", "amap", "ape", "msigdbr", "biomaRt", "Rfast", "caret", "recipes", "ipred", "timeDate", "monocle", "HH", "DT",
                     "prodlim", "lava", "SQUAREM", "gower", "pROC", "ModelMetrics", "limma", "VennDiagram", "psych", "corrplot", "cowplot", "gridExtra", "grid", "data.table", "DOSE",
                     "affycoretools", "plyr", "FactoMineR", "reshape2", "EGSEA", "GeneOverlap", "AnnotationDbi", "pheatmap", "rgl", "plotly", "pca3d", "R.utils", "purrr", "openxlsx",
                     "cluster", "graphics", "rafalib", "factoextra", "hrbrthemes", "reshape", "dplyr", "devtools", "CMScaller", "TCGAbiolinks", "CePa", "shiny", "regioneR", "normr",
                     "karyoploteR", "TxDb.Hsapiens.UCSC.hg38.knownGene", "sf", "rgdal", "ggVennDiagram", "ChIPseeker", "ReactomePA", "pander", "clusterProfiler", "ggupset", "Gviz",
                     "enrichplot", "ggnewscale", "pathview", "rtracklayer", "GenomicFormatExamples", "GenomicFeatures", "BSgenome.Hsapiens.UCSC.hg38", "rGADEM", "memes", "Biostrings",
                     "seqLogo", "BCRANK", "Rsamtools", "circlize", "readxl", "DiffBind", "ORFik", "TxDb.Hsapiens.UCSC.hg38.refGene", "clipr", "methylGSA", "seqTools", "stringr",
                     "methylationArrayAnalysis", "IlluminaHumanMethylation450kmanifest", "AnnotationHub", "ChIPpeakAnno", "satin", "Bolstad2", "MESS", "randomForest", "pROC", "ROCR",
                     "PharmacoGx", "S4Vectors", "SummarizedExperiment", "Biobase", "drda", "TeachingDemos", "areaplot", "scater", "easyGgplot2", "rstanemax", "broom", "DoseFinding",
                     "glmnet", "BBmisc", "gcrma", "httr", "jsonlite", "xml2", "GEOquery", "umap", "GenomicRanges", "Organism.dplyr", "BiocVersion", "a4Core", "Matrix.utils",
                     "org.Rn.eg.db")
lapply(neededlibraries, require, character.only = TRUE)



#################################
## 2. IMPORTAZIONE DATI GREZZI ##
#################################

SPECIES = "MOUSE"
if(SPECIES=="MOUSE"){initial = "ENSM"}
if(SPECIES=="HUMAN"){initial = "ENSG"}
targets = dir(pattern = "*.tab")

temp = read.table(targets[1], header = F, row.names=1, skip =4, stringsAsFactors = F)
nGenes = nrow(temp) # 32623
nSamples = length(targets) # 320
rat.rawcounts = matrix(data = NA, ncol = nSamples, nrow = nGenes)
pb <- txtProgressBar(min = 0, max = length(targets), initial = 0, char = "=", width = NA, style = 3)
for (i in 1:length(targets)) {
  setTxtProgressBar(pb, i)
  rat.rawcounts[,i] = read.table(targets[i], header = F, row.names=1, skip =4, stringsAsFactors = F)[,1]
}
rownames(rat.rawcounts) = rownames(temp)
colnames(rat.rawcounts) = targets
colnames(rat.rawcounts) <- gsub(colnames(rat.rawcounts), pattern = ".ReadsPerGene.out.tab", replacement = "")
dim(rat.rawcounts) # 32623   320


# SALVATAGGIO RDATA #
save(rat.rawcounts, file = "rat.rawcounts.RData")

# LOADING RDATA #
load("rat.rawcounts.RData")




####################################################
## 3. IMPORTAZIONI DELLE ANNOTAZIONI DEI CAMPIONI ##
####################################################


rat.annot = read.table("rat_annot.txt", header = T, row.names = 1, stringsAsFactors = F)
dim(rat.annot) # 320   4 

rat.rawcounts = rat.rawcounts[,rownames(rat.annot)]

keep = intersect(rownames(rat.annot), colnames(rat.rawcounts))
length(keep) # 320
rat.rawcounts = rat.rawcounts[,keep]



####################################
## 4. CREAZIONI DELLE ANNOTAZIONI ##
####################################

# Annotate Gene Symbols --------------------------

symbols <- mapIds(org.Rn.eg.db,
                  keys = rownames(rat.rawcounts),
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first"
)
entrez <- mapIds(org.Rn.eg.db,
                 keys = rownames(rat.rawcounts),
                 column = "ENTREZID",
                 keytype = "ENSEMBL",
                 multiVals = "first")
# columns(org.Rn.eg.db)
gene_annot = data.frame(SYMBOL=symbols, ENTREZ=entrez)
dim(gene_annot) # 32623     2



#####################################
## 5. NORMALIZZAZIONE CONTE GREZZE ##
#####################################

rat.dds = DESeqDataSetFromMatrix(countData = rat.rawcounts, colData = rat.annot, design = ~ 1)
rat.vst = vst(rat.dds)




###################################
## 6. CALCOLO DELLE CORRELAZIONI ##
###################################

tissues = c("adrenal", "brain", "heart", "kidney", "liver", "lung", "muscle", "spleen", "thymus", "testes", "uterus")

RAT = list() 
for(i in 1:length(tissues)){
  word = tissues[i]
  print(word)
  current.samples = rownames(rat.annot[grep(word, rat.annot$TISSUE),])
  print(length(current.samples))
  current.annot = rat.annot[current.samples,]
  current.vst = rat.vst[,current.samples]
  current.adj.vst = current.vst
  if(length(table(current.annot$SEX))>1) {
    print("both")
    batchcov = model.matrix(~ 0 + AGE, data = current.annot)
    assay(current.adj.vst) <- removeBatchEffect(x = assay(current.vst), batch = current.annot[colnames(current.adj.vst),]$SEX, design = batchcov)
  } 
  if(length(table(current.annot$SEX))==1) {
    print("single")
    current.adj.vst = current.vst
  }
  TRAIN = assay(current.adj.vst)
  train_size=dim(TRAIN)[2]
  lho_size=round(train_size/1, digits = 0)
  count=dim(TRAIN)[1]
  t = current.annot$AGE
  names(t) = rownames(current.annot)
  tx = current.annot$AGE
  TRAIN.PTIME = data.frame(PTIME = as.numeric(t), PTIMEX = tx, stringsAsFactors = F, row.names = colnames(TRAIN))
  corrtest = apply(TRAIN, 1, cor.test, TRAIN.PTIME$PTIME, method="pearson")
  coefs <- as.double(unlist(corrtest)[grep("estimate.cor",names(unlist(corrtest)))])
  names(coefs)<-gsub(pattern=".estimate.cor",x=rownames(TRAIN), replacement="")
  pvals <- as.double(unlist(corrtest)[grep("p.value",names(unlist(corrtest)))])
  names(pvals)<-gsub(pattern=".p.value",x=rownames(TRAIN), replacement="")
  corr_table <- data.frame(row.names = rownames(TRAIN),  pvalue = pvals, coef = coefs, stringsAsFactors = F)
  RAT[[i]]=corr_table
}
names(RAT) = tissues


# SALVATAGGIO RDATA #
save(RAT, file = "RAT.RData")

# LOADING RDATA #
load("RAT.RData")


# Annotazione Geni della tabella RAT
for(i in 1:length(RAT)){
  RAT[[i]]$SYMBOL = gene_annot[rownames(RAT[[i]]),]$SYMBOL
  RAT[[i]]$ENTREZ = gene_annot[rownames(RAT[[i]]),]$ENTREZ
}



tRAT = RAT

# Estrapolazione dei coefficiente di correlazione di ogni organo
DF = data.frame(row.names = rownames(tRAT[[1]]), matrix(ncol = length(tRAT), nrow = nrow(tRAT[[1]]), data = NA))
for(j in 1:length(tRAT)){
  DF[,j] = tRAT[[j]]$coef
  print(j)
}
names(DF) = names(tRAT)
dim(DF) # 32623    11


# Rename columns
names(DF)[1] <- "Adrenal"
names(DF)[2] <- "Brain"
names(DF)[3] <- "Heart"
names(DF)[4] <- "Kidney"
names(DF)[5] <- "Liver"
names(DF)[6] <- "Lung"
names(DF)[7] <- "Muscle"
names(DF)[8] <- "Spleen"
names(DF)[9] <- "Thymus"
names(DF)[10] <- "Testes"
names(DF)[11] <- "Uterus"

# Mediana dei coefficienti di correlazione su ogni tessuto
Median_Corr = rowMedians(data.matrix(DF), na.rm = T)
names(Median_Corr) = rownames(DF)
# Median_Corr = sort(Median_Corr, decreasing = T)
# plot(sort(Median_Corr))
# which(names(Median_Corr)==EDA2R)
Median_Corr = as.data.frame(Median_Corr)
dim(Median_Corr) # 32623     1

# Dataset Totale
DF_Total = data.frame(Median_DFR = Median_Corr, DF)
DF_Total$SYMBOL = gene_annot[rownames(DF_Total),]$SYMBOL
DF_Total = DF_Total[order(DF_Total$Median_Corr, decreasing = T),]
DF_Total$Rank_Median = 1:nrow(DF_Total)
DF_Total = as.data.frame(DF_Total)


# Add a ranking for the tisssue and Order the columns
DF_Total = DF_Total[order(DF_Total$Adrenal, decreasing = T),]
DF_Total$Rank_Adrenal = 1:nrow(DF_Total)
DF_Total = as.data.frame(DF_Total)

DF_Total = DF_Total[order(DF_Total$Brain, decreasing = T),]
DF_Total$Rank_Brain = 1:nrow(DF_Total)
DF_Total = as.data.frame(DF_Total)

DF_Total = DF_Total[order(DF_Total$Heart, decreasing = T),]
DF_Total$Rank_Heart = 1:nrow(DF_Total)
DF_Total = as.data.frame(DF_Total)

DF_Total = DF_Total[order(DF_Total$Kidney, decreasing = T),]
DF_Total$Rank_Kidney = 1:nrow(DF_Total)
DF_Total = as.data.frame(DF_Total)

DF_Total = DF_Total[order(DF_Total$Liver, decreasing = T),]
DF_Total$Rank_Liver = 1:nrow(DF_Total)
DF_Total = as.data.frame(DF_Total)

DF_Total = DF_Total[order(DF_Total$Lung, decreasing = T),]
DF_Total$Rank_Lung = 1:nrow(DF_Total)
DF_Total = as.data.frame(DF_Total)

DF_Total = DF_Total[order(DF_Total$Muscle, decreasing = T),]
DF_Total$Rank_Muscle = 1:nrow(DF_Total)
DF_Total = as.data.frame(DF_Total)

DF_Total = DF_Total[order(DF_Total$Spleen, decreasing = T),]
DF_Total$Rank_Spleen = 1:nrow(DF_Total)
DF_Total = as.data.frame(DF_Total)

DF_Total = DF_Total[order(DF_Total$Thymus, decreasing = T),]
DF_Total$Rank_Thymus = 1:nrow(DF_Total)
DF_Total = as.data.frame(DF_Total)

DF_Total = DF_Total[order(DF_Total$Testes, decreasing = T),]
DF_Total$Rank_Pancreas = 1:nrow(DF_Total)
DF_Total = as.data.frame(DF_Total)

DF_Total = DF_Total[order(DF_Total$Uterus, decreasing = T),]
DF_Total$Rank_Uterus = 1:nrow(DF_Total)
DF_Total = as.data.frame(DF_Total)

DF_Total = DF_Total[order(DF_Total$Median_Corr, decreasing = T),]
DF_Total = as.data.frame(DF_Total)

DF_Total <- DF_Total[, c(13, 1, 14, 2, 15, 3, 16, 4, 17, 5, 18, 6, 19, 7, 20, 8,
                         21, 9, 22, 10, 23, 11, 24, 12, 25)]
DF_Total = as.data.frame(DF_Total)


# SALVATAGGIO RDATA #
save(DF_Total, file = "DF_Total_RAT.RData")

# LOADING RDATA #
load("DF_Total_RAT.RData")


EDA2R = rownames(gene_annot)[which(gene_annot$SYMBOL=="Eda2r")] # "ENSRNOG00000013018"
EDA2R_Corr = DF_Total[which(rownames(DF_Total) == "ENSRNOG00000013018"),] # RANK 73




###################
## TISSUE SUBSET ##
###################

tissues_subset = c("adrenal", "brain", "heart", "kidney", "liver", "lung", "muscle", "spleen", "testes", "uterus")

DF_subset = DF[,c(1:8,10:11)]
dim(DF_subset) # 32623    10


# Mediana dei coefficienti di correlazione su ogni tessuto
Median_Corr = rowMedians(data.matrix(DF_subset), na.rm = T)
names(Median_Corr) = rownames(DF_subset)
# Median_Corr = sort(Median_Corr, decreasing = T)
# plot(sort(Median_Corr))
# which(names(Median_Corr)==EDA2R)
Median_Corr = as.data.frame(Median_Corr)
dim(Median_Corr) # 55359     1

# Dataset Totale
DF_Total_Sub = data.frame(Median_Corr = Median_Corr, DF_subset)
DF_Total_Sub$SYMBOL = gene_annot[rownames(DF_Total_Sub),]$SYMBOL
DF_Total_Sub = DF_Total_Sub[order(DF_Total_Sub$Median_Corr, decreasing = T),]
DF_Total_Sub$Rank_Median = 1:nrow(DF_Total_Sub)
DF_Total_Sub = as.data.frame(DF_Total_Sub)


# Order the columns
DF_Total_Sub = DF_Total_Sub[order(DF_Total_Sub$Adrenal, decreasing = T),]
DF_Total_Sub$Rank_Adrenal = 1:nrow(DF_Total_Sub)
DF_Total_Sub = as.data.frame(DF_Total_Sub)

DF_Total_Sub = DF_Total_Sub[order(DF_Total_Sub$Brain, decreasing = T),]
DF_Total_Sub$Rank_Brain = 1:nrow(DF_Total_Sub)
DF_Total_Sub = as.data.frame(DF_Total_Sub)

DF_Total_Sub = DF_Total_Sub[order(DF_Total_Sub$Heart, decreasing = T),]
DF_Total_Sub$Rank_Heart = 1:nrow(DF_Total_Sub)
DF_Total_Sub = as.data.frame(DF_Total_Sub)

DF_Total_Sub = DF_Total_Sub[order(DF_Total_Sub$Kidney, decreasing = T),]
DF_Total_Sub$Rank_Kidney = 1:nrow(DF_Total_Sub)
DF_Total_Sub = as.data.frame(DF_Total_Sub)

DF_Total_Sub = DF_Total_Sub[order(DF_Total_Sub$Liver, decreasing = T),]
DF_Total_Sub$Rank_Liver = 1:nrow(DF_Total_Sub)
DF_Total_Sub = as.data.frame(DF_Total_Sub)

DF_Total_Sub = DF_Total_Sub[order(DF_Total_Sub$Lung, decreasing = T),]
DF_Total_Sub$Rank_Lung = 1:nrow(DF_Total_Sub)
DF_Total_Sub = as.data.frame(DF_Total_Sub)

DF_Total_Sub = DF_Total_Sub[order(DF_Total_Sub$Muscle, decreasing = T),]
DF_Total_Sub$Rank_Muscle = 1:nrow(DF_Total_Sub)
DF_Total_Sub = as.data.frame(DF_Total_Sub)

DF_Total_Sub = DF_Total_Sub[order(DF_Total_Sub$Spleen, decreasing = T),]
DF_Total_Sub$Rank_Spleen = 1:nrow(DF_Total_Sub)
DF_Total_Sub = as.data.frame(DF_Total_Sub)

DF_Total_Sub = DF_Total_Sub[order(DF_Total_Sub$Testes, decreasing = T),]
DF_Total_Sub$Rank_Pancreas = 1:nrow(DF_Total_Sub)
DF_Total_Sub = as.data.frame(DF_Total_Sub)

DF_Total_Sub = DF_Total_Sub[order(DF_Total_Sub$Uterus, decreasing = T),]
DF_Total_Sub$Rank_Uterus = 1:nrow(DF_Total_Sub)
DF_Total_Sub = as.data.frame(DF_Total_Sub)

DF_Total_Sub = DF_Total_Sub[order(DF_Total_Sub$Median_Corr, decreasing = T),]
DF_Total_Sub = as.data.frame(DF_Total_Sub)

DF_Total_Sub <- DF_Total_Sub[, c(12, 1, 13, 2, 14, 3, 15, 4, 16, 5, 17, 6, 18, 7, 19, 8,
                         20, 9, 21, 10, 22, 11, 23)]
DF_Total_Sub = as.data.frame(DF_Total_Sub)



# SALVATAGGIO RDATA #
save(DF_Total_Sub, file = "DF_Total_Sub_RAT.RData")

# LOADING RDATA #
load("DF_Total_Sub_RAT.RData")













