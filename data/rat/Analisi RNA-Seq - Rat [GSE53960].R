

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
                  multiVals = "first")
entrez <- mapIds(org.Rn.eg.db,
                 keys = rownames(rat.rawcounts),
                 column = "ENTREZID",
                 keytype = "ENSEMBL",
                 multiVals = "first")
ensembl <- mapIds(org.Rn.eg.db,
                  keys = row.names(rat.rawcounts),
                  column = "ENSEMBL",
                  keytype = "ENSEMBL",
                  multiVals = "first")
chromosomes <- mapIds(org.Rn.eg.db,
                      keys = row.names(rat.rawcounts),
                      column = "CHR",
                      keytype = "ENSEMBL",
                      multiVals = "first")

gene_annot = data.frame(SYMBOL = symbols, ENTREZ = entrez, ENSEMBL = ensembl, CHROMOSOME = chromosomes, row.names = row.names(rat.rawcounts))
dim(gene_annot) # 32623     4

# # columns(org.Rn.eg.db)
# gene_annot = data.frame(SYMBOL=symbols, ENTREZ=entrez)
# dim(gene_annot) # 32623     2





#####################################
## 5. NORMALIZZAZIONE CONTE GREZZE ##
#####################################

rat.dds = DESeqDataSetFromMatrix(countData = rat.rawcounts, colData = rat.annot, design = ~ 0 + GROUP)
# rat.dds <- rat.dds[rowSums(counts(rat.dds)) > 0, ]
# rat.dds.cpm <- cpm(rat.dds, normalized.lib.sizes = TRUE)
# dim(rat.dds.cpm) # 52293   198
# 
# ConteperMilione = dds.cpm
# ConteperMilione_forgene <- data.frame(SYMBOL = gene_annot[rownames(ConteperMilione),]$SYMBOL, 
#                                       ENTREZ = gene_annot[rownames(ConteperMilione),]$ENTREZ,
#                                       ENSEMBL = rownames(ConteperMilione),
#                                       CHROMOSOME = gene_annot[rownames(ConteperMilione),]$CHROMOSOME, 
#                                       ConteperMilione)
# ConteperMilione_forgene = ConteperMilione_forgene[which(!is.na(ConteperMilione_forgene$SYMBOL)),]
# ConteperMilione_forgene <- ConteperMilione_forgene[!duplicated(ConteperMilione_forgene$SYMBOL),]
# rownames(ConteperMilione_forgene) = ConteperMilione_forgene$SYMBOL
# ConteperMilione_forgene = as.data.frame(ConteperMilione_forgene)
# ConteperMilione_forgene = ConteperMilione_forgene[order(ConteperMilione_forgene$CHROMOSOME),] 
# ConteperMilione_forgene = ConteperMilione_forgene[order(ConteperMilione_forgene$SYMBOL),] 
# dim(ConteperMilione_forgene) # 32700   202
# 
# # Per non avere le colonne di: SYMBOL - ENTREZ - ENSEMBL - CHROMOSOME:
# ConteperMilione_forgene_selection = data.frame(ConteperMilione_forgene[,c(5:dim(ConteperMilione_forgene)[2])])
# ConteperMilione_forgene_selection = as.matrix(ConteperMilione_forgene_selection)
# dim(ConteperMilione_forgene_selection) # 32700   198
# 
# 
# log.ConteperMilione <- log2(1 + dds.cpm)
# log.ConteperMilione_forgene <- data.frame(SYMBOL = gene_annot[rownames(log.ConteperMilione),]$SYMBOL, 
#                                           ENTREZ = gene_annot[rownames(log.ConteperMilione),]$ENTREZ,
#                                           ENSEMBL = rownames(log.ConteperMilione),
#                                           CHROMOSOME = gene_annot[rownames(ConteperMilione),]$CHROMOSOME,
#                                           log.ConteperMilione)
# log.ConteperMilione_forgene = log.ConteperMilione_forgene[which(!is.na(log.ConteperMilione_forgene$SYMBOL)),]
# log.ConteperMilione_forgene <- log.ConteperMilione_forgene[!duplicated(log.ConteperMilione_forgene$SYMBOL),]
# rownames(log.ConteperMilione_forgene) = log.ConteperMilione_forgene$SYMBOL
# log.ConteperMilione_forgene = as.data.frame(log.ConteperMilione_forgene)
# log.ConteperMilione_forgene = log.ConteperMilione_forgene[order(log.ConteperMilione_forgene$CHROMOSOME),]
# log.ConteperMilione_forgene = log.ConteperMilione_forgene[order(log.ConteperMilione_forgene$SYMBOL),] 
# dim(log.ConteperMilione_forgene) # 32700   202
# 
# # Per non avere le colonne di: SYMBOL - ENTREZ - ENSEMBL - CHROMOSOME:
# log.ConteperMilione_forgene_selection = data.frame(log.ConteperMilione_forgene[,c(5:dim(log.ConteperMilione_forgene)[2])])
# log.ConteperMilione_forgene_selection = as.matrix(log.ConteperMilione_forgene_selection)
# dim(log.ConteperMilione_forgene_selection) # 32700   198




rat.vst = vst(rat.dds)

vsd_assay.rat <- assay(rat.vst)
vsd_assay.rat.df <- data.frame(vsd_assay.rat)


EDA2R.rat = vsd_assay.rat.df["ENSRNOG00000013018",]
# EDA2R = 2^vsd_assay.rat.df["ENSRNOG00000013018",]

EDA2R_df.rat <- data.frame(EDA2R.rat)
EDA2R_df.rat <- t(EDA2R_df.rat)

age <- rat.annot["GROUP"]

EDA2R_df.tot <- merge(EDA2R_df.rat, age, by="row.names")
rownames(EDA2R_df.tot) <- EDA2R_df.tot$Row.names
EDA2R_df.tot <- EDA2R_df.tot[,2:3]
colnames(EDA2R_df.tot)[1] <- "EXP"



#########################
## 6. BOX PLOT - EDA2R ##
#########################

tissues = c("adrenal", "brain", "heart", "kidney", "liver", "lung", "muscle", "spleen", "thymus", "testes", "uterus")
library(ggpubr)

PLOT_list = list() 
# for(i in 1:length(tissues)){
  i=9
  word = tissues[i]
  print(word)
  current.samples = rownames(rat.annot[grep(word, rat.annot$TISSUE),])
  print(length(current.samples))
  current.annot = rat.annot[current.samples,]
  current.vst = rat.vst[,current.samples]
  current_assay <- assay(current.vst)
  current_assay.df <- data.frame(current_assay)
  Eda2r = current_assay.df["ENSRNOG00000013018",]
  Eda2r = t(Eda2r)
  
  age <- rat.annot["GROUP"]
  
  
  Eda2r_df.tot <- merge(Eda2r, age, by="row.names")
  rownames(Eda2r_df.tot) <- Eda2r_df.tot$Row.names
  Eda2r_df.tot <- Eda2r_df.tot[,2:3]
  colnames(Eda2r_df.tot)[1] <- "EXP"

  Eda2r_df.tot$GROUP = factor(Eda2r_df.tot$GROUP, levels = c("juvenile","adolescent","adult","elderly"))
  
  my_comparisons <- list(c("juvenile", "adolescent"),c("juvenile", "adult"),c("juvenile", "elderly"))
  
  
  p_gene = ggplot(Eda2r_df.tot, aes(x=GROUP, y=EXP, fill=GROUP)) +
    geom_boxplot(alpha = 0.4, show.legend = FALSE, notch = FALSE) +
    geom_jitter( aes( fill=GROUP), shape = 21,  size=2, alpha=1,  position = position_jitter(width = .1)) +
    stat_compare_means(comparisons = my_comparisons,  color = "black", bracket.size = 0.5, ) +
    theme_bw() + theme(axis.text.x = element_text(angle=0,vjust = 0)) +
    theme(axis.text.x=element_text(angle=0,vjust = 0))
  p_gene

  
# }  
  
  
# 2 cresce con l'eta
# 7 è crescente, juv e ado non sono significativi
# 9 è crescente, juv e ado non sono significativi

  

## EDA2R ##
  
EDA2R_df.tot$GROUP = factor(EDA2R_df.tot$GROUP, levels = c("juvenile","adolescent","adult","elderly"))


my_comparisons <- list(c("juvenile", "adolescent"),c("juvenile", "adult"),c("juvenile", "elderly"))

p_gene = ggplot(EDA2R_df.tot, aes(x=GROUP, y=EXP, fill=GROUP)) +
  #stat_compare_means(comparisons = my_comparisons) +
  #scale_shape_manual(name = "df", values = c(21)) +
  #scale_y_continuous(limits = c(40,180)) + 
  geom_boxplot(alpha = 0.4, show.legend = FALSE, notch = FALSE) +
  #geom_boxplot(aes(fill=AGE), width = 0.20, outlier.size =0) +
  #geom_point(data = DFX, size=4, aes(fill=I("yellow")), colour="black",pch=21, stroke = 0) +
  geom_jitter( aes( fill=GROUP), shape = 21,  size=2, alpha=1,  position = position_jitter(width = .1)) +
  #stat_compare_means(comparisons = my_comparisons,  color = "black", bracket.size = 0.5, ) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=-90,vjust = 0))
p_gene



p_gene + stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test') + 
  stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test') +
  stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test') 
  
p_gene +stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test", ref.group = "elderly")


ggsave('plot.svg', p_gene, width = 6, height = 8)
readtext("plot.svg")-> temp
write_clip(temp$text, object_type = "character")



################################################################
## 8. CALCOLO DELLO SCORE DEI GENI DIFFERENZIALMENTE ESPRESSI ##
################################################################

rat.dds <- estimateSizeFactors(rat.dds)
rat.dds <- DESeq(rat.dds, parallel = T, BPPARAM = SnowParam(4))


# SALVATAGGIO RDATA #
save(rat.dds, file = "rat.dds.RData")

# LOADING RDATA #
# load("rat.dds.RData")


rat.dds.design = model.matrix(~ 0 + GROUP, data = rat.annot)
colnames(rat.dds.design) = gsub(x = colnames(rat.dds.design), pattern = "GROUP", replacement = "")

contlist = list()
contlist[[1]] =  makeContrasts(adolescent - juvenile, levels = rat.dds.design)
contlist[[2]] =  makeContrasts(adult - juvenile, levels = rat.dds.design)
contlist[[3]] =  makeContrasts(elderly - juvenile, levels = rat.dds.design)



comparison_names = c("Juvenile vs Adolescent", "Juvenile vs Adult", "Juvenile vs Elderly")

names(contlist) = comparison_names



##############################################################
## 9. ANALISI DEI GENI DIFFERENZIALMENTE ESPRESSI CON DESeq ##
##############################################################

## Calcolo "res.treatment" ##

res.treatment = list()
for(i in 1:length(contlist)) {
  print(i)
  results.treatment = results(object = rat.dds, contrast = contlist[[i]], independentFiltering = TRUE, alpha = 0.1)
  results.treatment = data.frame(SYMBOL = gene_annot[rownames(results.treatment), 1], 
                                 ENTREZ = gene_annot[rownames(results.treatment), 2],
                                 ENSEMBL = gene_annot[rownames(results.treatment), 3],
                                 CHROMOSOME = gene_annot[rownames(results.treatment), 4],
                                 results.treatment)
  results.treatment = results.treatment[which(!is.na(results.treatment$SYMBOL)),]
  results.treatment <- results.treatment[!duplicated(results.treatment$SYMBOL),]
  rownames(results.treatment) = results.treatment$SYMBOL
  results.treatment = results.treatment[order(results.treatment$padj),] 
  results.treatment = results.treatment[order(results.treatment$SYMBOL),] 
  results.treatment <- results.treatment[,c(1,2,3,4,5,7,8,6,9,10)]   
  res.treatment[[i]] <- results.treatment
  
}

names(res.treatment) = comparison_names

# SALVATAGGIO RDATA #
save(res.treatment, file = "res.treatment_rat.RData")

# LOADING RDATA #
# load("res.treatment_rat.RData")


## Trasformazione p-value nulli ##

res.treat = list()
for(j in 1:length(contlist)) {
  print(j)
  r.treatment = res.treatment[[j]][which(!is.na(res.treatment[[j]]$padj)),]  
  # I p-value dei geni che hanno un valore pari 0, vengono sostituiti con il valore del p-value del gene antecedente non nullo, in modo da avere un valore finito
  length = dim(r.treatment)[1]
  for(i in 1:length) { 
    if (r.treatment$pvalue[i] == 0) {r.treatment$pvalue[i] = min(r.treatment$pvalue[r.treatment$pvalue != min(r.treatment$pvalue)])}
  }
  length = dim(r.treatment)[1]
  for(i in 1:length) { 
    if (r.treatment$padj[i] == 0) {r.treatment$padj[i] = min(r.treatment$padj[r.treatment$padj != min(r.treatment$padj)])}
  }
  dim(r.treatment)  
  res.treat[[j]] <- r.treatment
}

names(res.treat) = comparison_names

# SALVATAGGIO RDATA #
save(res.treat, file = "res.treat_rat.RData")

# LOADING RDATA #
# load("res.treat_rat.RData")


## Calcolo "results.sig" ##

res.sig = list()
for(j in 1:length(contlist)) {
  print(j)
  results.sig = res.treat[[j]][which(res.treat[[j]]$padj <= 0.1),]
  results.sig = results.sig[,-5]  
  results.sig = results.sig[,-5]
  results.sig = results.sig[order(abs(results.sig$log2FoldChange), decreasing = TRUE), ] 
  colnames(results.sig) = c("SYMBOL","ENTREZ","ENSEMBLE","CHROMOSOME","STAT","FOLD-CHANGE","P-VALUE","P-ADJ")
  dim(results.sig)
  res.sig[[j]] = results.sig
}

names(res.sig) = comparison_names


# SALVATAGGIO RDATA #
save(res.sig, file = "RNA-Seq Data (CdLS) - RAT.RData")

# LOADING RDATA #
# load("RNA-Seq Data (CdLS) - RAT.RData")



#######################
## 12. ARRICCHIMENTO ##
#######################

# Riduzione Dataframe "dds":

dim(rat.dds)  # 32623   320
rat.dds.sub = rat.dds[which(!is.na(gene_annot[rownames(rat.dds),]$SYMBOL)),]
dim(rat.dds.sub)  # 20365   320
dds.vst <- vst(rat.dds, blind = TRUE)
temp = rowMeans(assay(dds.vst)[rownames(rat.dds.sub),])
temp = temp[order(temp, decreasing = T)]
keep = names(temp[which(!duplicated(gene_annot[names(temp),]$SYMBOL))])
length(keep) # 20346
rat.dds.sub = rat.dds.sub[keep,]
dim(rat.dds.sub) # 20346   320


# Caricamento Confronti:

comlist = list()
for(i in 1:length(contlist)){
  print(i)
  temp = results(rat.dds.sub, contrast =  contlist[[i]], independentFiltering = TRUE, alpha = 0.1, parallel = T, BPPARAM = SnowParam(4))
  temp = data.frame(gene_annot[rownames(temp),],temp, stringsAsFactors = F)
  temp = temp[order(temp$pvalue),]
  comlist[[i]] = temp
}
names(comlist) = comparison_names

# comlist = res.treat



# SALVATAGGIO RDATA #
save(comlist, file = "comlist_rat.RData")

# LOADING RDATA #
# load("comlist_rat.RData")


# COMPARISON #

comparison = list()
for(i in 1:length(contlist)) {
  print(i)
  p = as.data.frame(res.treatment[[i]][,c(8,9,10)])
  colnames(p) = c(paste("log2FC_", as.character(i), sep = ""), paste("PValue_", as.character(i), sep = ""), paste("FDR_", as.character(i), sep = ""))
  comparison[[i]] = p
}
names(comparison) = comparison_names

# # comp = do.call(cbind, comparison)
# comp <- data.frame(matrix(ncol = 0, nrow = nrow(ConteperMilione_forgene)))
# for (i in 1:length(contlist)) {
#   print(i)
#   p = comparison[[i]]
#   comp <- data.frame(comp,p)
# }

Data_Expression = data.frame(comparison)
dim(Data_Expression) # 20346    12

Eda2r <- Data_Expression["Eda2r",]


df <- data.frame(group=c("2vs6", "2vs21", "2vs104"),
                 len=c(-0.01324535, -0.2698127, 29.5))
head(df)

library(ggplot2)
# Basic barplot
p<-ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity")
p

p <- p+scale_fill_brewer(palette="Blues")
p + theme(legend.position="none")




















## RANKING ANALISI DIFFERENZIALE (ALL GENES) COMPARISON 3 ##

Data_Expression$SIG = -10*log10(Data_Expression$Juvenile.vs.Elderly.FDR_3)
Data_Expression$SIGX = Data_Expression$SIG*sign(Data_Expression$Juvenile.vs.Elderly.log2FC_3)
Data_Expression = Data_Expression[order(Data_Expression$SIGX, decreasing = T),]
#gse = gse[order(gse$t, decreasing = T),]
#gse = gse[order(gse$logFC, decreasing = T),]

Data_Expression$GROUP = "GROUP"

p_gene = ggplot(Data_Expression, aes(x=GROUP, y=SIGX, fill=GROUP)) +
  #stat_compare_means(comparisons = my_comparisons) +
  #scale_shape_manual(name = "df", values = c(21)) +
  #scale_y_continuous(limits = c(40,180)) + 
  geom_violin(aes(fill=GROUP), width = 0.4, trim = F) +
  geom_boxplot(aes(fill=GROUP), width = 0.08, outlier.size =1) +
  #geom_point(data = DFX, size=4, aes(fill=I("yellow")), colour="black",pch=21, stroke = 0) +
  #geom_jitter( aes( fill=GROUP), shape = 21,  size=1, alpha=1,  position = position_jitter(width = .1)) +
  #stat_compare_means(comparisons = my_comparisons,  color = "black", bracket.size = 0.5, ) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=-90,vjust = 0))
p_gene

ggsave('plot.svg', p_gene, width = 4, height = 4)
readtext("plot.svg")-> temp
write_clip(temp$text, object_type = "character")





#############################
## 6. VOLCANO PLOT - EDA2R ##
#############################

res = comlist$`Juvenile vs Elderly`
res = res[which(abs(res$log2FoldChange)<4),]
res = res[which(res$baseMean>100),]

threshold_FC = 1
threshold_PVALUE = 0.05
Genes_Symbol = res$SYMBOL
Genes_Ensembl = rownames(res)
log2FoldChange = res$log2FoldChange
FDR = -10*log10(res$padj)

my.list_UP = as.data.frame(Genes_Symbol)[which(res$log2FoldChange > 0),]
my.list_DOWN = as.data.frame(Genes_Symbol)[which(res$log2FoldChange < 0),]

gene.selected.UP = intersect(intersect(which(res$SYMBOL %in% my.list_UP), which(res$log2FoldChange >= 1)), which(res$padj <= 0.05))
gene.selected.DOWN = intersect(intersect(which(res$SYMBOL %in% my.list_DOWN), which(res$log2FoldChange <= -1)), which(res$padj <= 0.05))

cex = rep(1.2, length(Genes_Symbol))
cex[gene.selected.UP] = 1
cex[gene.selected.DOWN] = 1

inner = rep("grey60", length(Genes_Symbol))
inner[gene.selected.UP] = "red4"
inner[gene.selected.DOWN] = "navy"

boundary = rep("grey", length(Genes_Symbol)) 
boundary[gene.selected.UP] = "red"
boundary[gene.selected.DOWN] = "blue"

Matrix = data.frame(SYMBOL = Genes_Symbol, ENSEMBL = Genes_Ensembl, LOG2FOLDCHANGE = log2FoldChange, FDR = FDR, CEX = cex, INNER = inner, BOUNDARY = boundary)

p = plot(Matrix$LOG2FOLDCHANGE, Matrix$FDR, panel.first = grid(), main = "Volcano Plot", xlab = "log2(Fold-Change)", ylab = "-log10(FDR)", cex = Matrix$CEX,
         pch = 21, col = Matrix$BOUNDARY, bg = Matrix$INNER)
abline(v = 0)
abline(v = c(-threshold_FC, threshold_FC), col = "brown")
abline(h = -10*log10(threshold_PVALUE), col = "brown")

print_text = FALSE
UP = as.data.frame(Genes_Symbol)[which(res$log2FoldChange > 1 & res$padj < threshold_PVALUE),]
DOWN = as.data.frame(Genes_Symbol)[which(res$log2FoldChange < -1 & res$padj < threshold_PVALUE),]
gene.selected.text = c(UP, DOWN)

if (print_text == TRUE) {
  gene.selected.text = which(res$SYMBOL %in% gene.selected.text)
  text(log2FoldChange[gene.selected.text], FDR[gene.selected.text], lab = Genes_Symbol[gene.selected.text], cex = 0.7, pos = 3, font = 2)
  print(p)
}  # Salvataggio pdf: 5 y - 7 x ; 900 - 400


prova <- Matrix[order(Matrix$LOG2FOLDCHANGE, decreasing = T),]




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


EDA2R = rownames(gene_annot)[which(gene_annot$SYMBOL=="Eda2r")] # "ENSRNOG00000013018"
EDA2R_Corr = DF_Total_Sub[which(rownames(DF_Total_Sub) == "ENSRNOG00000013018"),] # RANK 118




####################################
## 5. ANALISI DIFFERENZIALE [GEO] ##
####################################

# LINK ANALISI DIFFERENZIALE: https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE47792&platform=GPL14844


# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Data plots for selected GEO samples
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE47792", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL14844", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]


# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("11113333222200001111333322220000111133332222000011",
               "11333322220000111133332222000011113333222200001111",
               "33332222000011113333222200001111333322220000111133",
               "33222200001111333322220000111133332222000011113333",
               "22220000111133332222000011113333222200001111333322",
               "22000011113333222200001111333322220000111133332222",
               "00001111333322220000XXXXXXXXXXXXXXX")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]


# # log2 transformation
# ex <- exprs(gset)
# qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# LogC <- (qx[5] > 100) ||
#   (qx[6]-qx[1] > 50 && qx[2] > 0)
# if (LogC) { ex[which(ex <= 0)] <- NaN
# exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("elderly","juvenile","adult","adolescent"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""), paste(groups[2],"-",groups[3],sep=""), paste(groups[2],"-",groups[4],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)








ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
par(mar=c(7,4,2,1))
title <- paste ("GSE47792", "/", annotation(gset), sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
dev.off()

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE47792", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

# mean-variance trend
ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE47792")

# UMAP plot (multi-dimensional scaling)
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)







