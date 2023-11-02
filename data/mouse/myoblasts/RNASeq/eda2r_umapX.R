setwd("~/Desktop/REVX/RNASeq/TAB Files/")

library(BiocParallel)
rawcounts = list()
targets = dir(pattern = "*.tab", full.names = F)

genes_to_import = rownames(read.table(targets[1], header = F, row.names=1, stringsAsFactors = F, skip = 4)[,1:2])
nSamples = length(targets)
rawcounts = matrix(data = NA, ncol = nSamples, nrow = length(genes_to_import))
colnames(rawcounts)=targets
rownames(rawcounts)=genes_to_import
pb <- txtProgressBar(min = 0, max = length(targets), initial = 0, char = "=",width = NA, style = 3)
for(i in 1:length(targets)){
  setTxtProgressBar(pb, i)
  rawcounts[,i] = read.table(targets[i], header = F, row.names=1, stringsAsFactors = F, skip = 4)[,1:3][genes_to_import,3]
}
colnames(rawcounts) = targets
colnames(rawcounts) = gsub(colnames(rawcounts), pattern = ".tab", replacement = "")
keep=grep(cbind(unlist(strsplit(rownames(rawcounts),'.',fixed=TRUE))), pattern = "ENSM")
keep_names = cbind(unlist(strsplit(rownames(rawcounts),'.',fixed=TRUE)))[keep]
rownames(rawcounts) = keep_names


# PCA FARE
#DE HEATMAP


library(Seurat)
library(ggplot2)
#library(clustree)

## Create list of Seurat Object ------------------------------------------------------------------------

datasets_rawcounts = list()
datasets_rawcounts$EDA2R = rawcounts

rawsamples = datasets_rawcounts
for(i in 1:length(rawsamples)){
  print(i)
  currName = names(rawsamples)[i]
  if(i==1){
    samples = list()
    samples[[i]] = CreateSeuratObject(rawsamples[[i]], project = currName)
  }
  if(i>1){samples[[i]] = CreateSeuratObject(rawsamples[[i]], project = currName)}
}
names(samples) = names(rawsamples)

# retrieve gene annotations ------------------------------
SPECIES = "MOUSE"
if(SPECIES=="MOUSE"){
  library(org.Mm.eg.db)
  Db = org.Mm.eg.db
}
if(SPECIES=="HUMAN"){
  library(org.Hs.eg.db)
  Db = org.Hs.eg.db
}
if(SPECIES=="RAT"){
  library(org.Rn.eg.db)
  Db = org.Rn.eg.db
}
symbols <- mapIds(Db,
                  keys = row.names(rawsamples[[1]]),
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first"
)
entrez <- mapIds(Db,
                 keys = row.names(rawsamples[[1]]),
                 column = "ENTREZID",
                 keytype = "ENSEMBL",
                 multiVals = "first")
gene_annot = data.frame(SYMBOL=symbols, ENTREZ=entrez)
rownames(gene_annot) <- rownames(rawsamples[[1]])


rawcounts = rawsamples$EDA2R

library(DESeq2)
library(edgeR)
library(ggrepel)
library(readtext)
library(svglite)
samples = readRDS(file = "../eda2r_annot.RDS")
dds <- DESeqDataSetFromMatrix(countData = rawcounts, colData = samples, design = ~ 0 + GROUP + BIOREP)
# Generate alternative quantifications
dds.cpm <- cpm(dds, normalized.lib.sizes = TRUE)
dds.log2cpm <- log2(1+dds.cpm)
dds.vst <- vst(dds, blind = TRUE)
dds.assay <- assay(dds.vst)

#PCA ANALYSIS ------------------------------
#plotPCA(dds_vst, intgroup = "GROUP", ntop = 2000)
# Perform PCA analysis using top 2000 genes
ntop = 2000
rv = rowVars(dds.assay) # calcola varianza di ogni gene
rv_ordered = rownames(dds.assay)[order(rv, decreasing = TRUE)] # ordina i gene symbol in base alla varianza
selected = rv_ordered[1:ntop]
pca = prcomp( t(dds.assay[selected,]) )
percentVar <- round(100* pca$sdev^2/sum(pca$sdev^2) )
df.pca = data.frame(pca$x[,1:3], samples[rownames(pca$x),], LABEL = rownames(pca$x))
percentVar <- round(100* pca$sdev^2/sum(pca$sdev^2) )
p1 <- ggplot(df.pca, aes(x = PC1, y = PC2, color=GROUP)) +
  geom_point(size=3.5, aes(fill=GROUP ), colour="black",pch=21, stroke = 0.8) +
  coord_fixed(ratio = 1) +
  xlab(paste0("PC1"," : ",percentVar[1]," % variance")) +
  ylab(paste0("PC2"," : ",percentVar[2]," % variance")) +
  coord_fixed(ratio = 6/5)  +
  geom_text_repel(data = df.pca, aes(label = LABEL), nudge_y = 5, nudge_x = 5, segment.size = 0.4, hjust = 1, point.padding = 0.5, arrow = arrow(length = unit(0.08, "inches"), angle = 15, type = "closed", ends = "first")) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA, colour = "black", size = 1), axis.text=element_text(size=10, face="plain"), axis.title=element_text(size=12,face="plain",family = "Arial"), axis.ticks.x = element_line(colour = "black"))
p1
ggsave('plot.svg', p1, width = 7, height = 7)
readtext("plot.svg")-> temp
write_clip(temp$text, object_type = "character")

# Hierarchical Clustering ------------------------------
library(amap)
library(ape)
data <- scale(t(dds.assay))
amap_dist = amap::Dist(data,method="euclidean", nbproc = 4)
amap_hc <- hcluster(amap_dist, link = "ave", nbproc = 4)
hc = amap_hc
hc
d <- as.dendrogram(hc)
colTips= c( rep("black", ncol(dds.assay)  ) )
names(colTips) = colnames(dds.assay)
colTips[rownames(samples)[which(samples$GROUP=="LMNA_Y")]] = "blue"
colTips[rownames(samples)[which(samples$GROUP=="WT_Y")]] = "darkgreen"
colTips[rownames(samples)[which(samples$GROUP=="WT_O")]] = "red"
plot(as.phylo(hc), type = "phylogram", cex = 0.8, edge.color = "black", edge.width = 2, edge.lty = 1, tip.color = colTips)

# Define design ------------------------------
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, parallel = T, BPPARAM = MulticoreParam(4))
dds.design = model.matrix(~ 0 + GROUP + BIOREP, data=samples)
colnames(dds.design) = gsub(x = colnames(dds.design), pattern = "GROUP", replacement = "")

convert_human_to_mouse <- function(gene_list) {
  output = c()
  mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  for(gene in gene_list) {
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "human"))[['DB.Class.Key']]
    if( !identical(class_key, integer(0)) ) {
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(human_gene in human_genes) {
        output = rbind(c(gene, human_gene), output)
      }
    }
  }
  return (output)
}

convert_mouse_to_human <- function(gene_list) {
  output = c()
  mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  for(gene in gene_list) {
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "mouse, laboratory"))[['DB.Class.Key']]
    if( !identical(class_key, integer(0)) ) {
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes) {
      output = rbind(c(gene, human_gene), output)
      }
    }
  }
return (output)
}

library(dplyr)
keep_annot = gene_annot[which(!is.na(gene_annot$SYMBOL)),]
conv_annot = convert_mouse_to_human(keep_annot$SYMBOL)
keep2_annot = keep_annot[which(keep_annot$SYMBOL %in% conv_annot[,1]),]
keep2_annot$ENSEMBL = rownames(keep2_annot)
rownames(conv_annot) = conv_annot[,1]
keep2_annot$ORTHO = conv_annot[keep2_annot$SYMBOL,2]

out <- EDASeq::getGeneLengthAndGCContent(
  id = rownames(keep2_annot),
  org = 'mm',
  mode = c('org.db'))

keep2_annot$length = out[rownames(keep2_annot),1]
keepcounts = rawcounts[rownames(keep2_annot),]


samples = readRDS(file = "../eda2r_annot.RDS")

mergecounts = data.frame(
  row.names = rownames(keepcounts),
  EDA2R_BRep1 = rowMeans(keepcounts[,1:3]),
  EDA2R_BRep2 = rowMeans(keepcounts[,4:6]),
  EDA2R_BRep3 = rowMeans(keepcounts[,7:9]),
  GFP_BRep1 = rowMeans(keepcounts[,10:12]),
  GFP_BRep2 = rowMeans(keepcounts[,13:15]),
  GFP_BRep3 = rowMeans(keepcounts[,16:18])
)

library(RNAAgeCalc)
mergeFPKM = count2FPKM(mergecounts, genelength = keep2_annot$length, idtype = "ENSEMBL")

unique_annot = keep2_annot[which(!duplicated(keep2_annot$ORTHO)),]
keptFPKM = mergeFPKM[rownames(unique_annot),]
rownames(keptFPKM) = unique_annot$ORTHO

RES <- predict_age(
  exprdata = keptFPKM,
  tissue = "muscle",
  exprtype = "FPKM",
  idtype = "SYMBOL"
)
plot(RES$RNAAge)

expdesign = data.frame(
  row.names = colnames(mergeFPKM),
  BATCH = c("EXP1", "EXP2", "EXP3", "EXP1", "EXP2", "EXP3"), 
  GROUP = c("EDA2R", "EDA2R", "EDA2R", "GFP", "GFP", "GFP")
)

AREP1 = mean(c(RES$RNAAge[1], RES$RNAAge[4]))
AREP2 = mean(c(RES$RNAAge[2], RES$RNAAge[5]))
AREP3 = mean(c(RES$RNAAge[3], RES$RNAAge[6]))

RES$modAGE = NA
RES$modAGE[1] = RES$RNAAge[1]-AREP1
RES$modAGE[2] = RES$RNAAge[2]-AREP2
RES$modAGE[3] = RES$RNAAge[3]-AREP3
RES$modAGE[4] = RES$RNAAge[4]-AREP1
RES$modAGE[5] = RES$RNAAge[5]-AREP2
RES$modAGE[6] = RES$RNAAge[6]-AREP3

head(RES) 

NF = mean(RES$modAGE[4:6])
RES$fcAGE = RES$modAGE-NF

df = cbind(RES, expdesign)
library(ggplot2)
library(gghalves)
p1 <- ggplot(df, aes(y = fcAGE, x = GROUP) ) +
#geom_point(data = correlations_nosig, size=3.5, colour="black", pch=20, stroke = 0.1) +
#geom_point(data = correlations_sig, size=3.5, colour="red", pch=20, stroke = 0.1) +
#geom_point(data = correlations_cmn, size=3.5, colour="yellow", pch=20, stroke = 0.1) +
#geom_violin(side = "r", draw_quantiles = F, trim = F, scale = "width", aes(fill = GROUP)) +
geom_boxplot( aes(fill = GROUP))  +
geom_point(data = df, size=8, colour="blue", pch=20, stroke = 0.4) + 
#stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test', method.args = list(alternative = "greater")) +
#stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test', symnum.args = symnum.args) +
#geom_boxplot() +
#geom_half_violin() +
#geom_jitter(data = correlations_nosig, pch = 20,  size=0.2, alpha=1,  position = position_jitter(width = 0.1)) +
theme_bw() +
#coord_flip() +
theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA, colour = "black", size = 1), axis.text=element_text(size=10, face="plain"), axis.title=element_text(size=12,face="plain",family = "Arial"), axis.ticks.x = element_line(colour = "black"))
p1









adjFPKM = removeBatchEffect(mergeFPKM, batch = expdesign$BATCH)
range(adjFPKM, na.rm = T)






reFPKM = rescale(adjFPKM, c(0, 4762.843417))









design = 

adjFPKM = removeBatchEffect(mergeFPKM, batch = c("EXP1", "EXP2", "EXP3", "EXP1", "EXP2", "EXP3"), design =  )
adjFPKM = adjFPKM+5



RES$ADJ = removeBatchEffect(RES$RNAAge, batch = c("EXP1", "EXP2", "EXP3", "EXP1", "EXP2", "EXP3"))

RES

RES$NORM = RES$RNAAge-mean(RES$RNAAge[4:6])



help(count2FPKM)





rownames(keepcounts) = keep2_annot$ORTHO
keepcounts = keepcounts[which(!duplicated(rownames(keepcounts))),]



dim(keepcounts)



library(RNAAgeCalc)
predict_age(
  stype = "all",
  exprdata = cpm(mergecounts),
  tissue = "muscle",
  exprtype = "counts",
  signature = "DESeq2",
  idtype = "SYMBOL"
) -> RES
plot(RES$RNAAge)

assay(current.vst) <- removeBatchEffect(x = assay(current.vst), batch = current.annot[colnames(current.vst),]$SEX, design = batchcov)






# Define Constrasts ------------------------------
contlist = list()
contlist[[1]] =  makeContrasts(EDA2R-GFP,levels = dds.design)
names(contlist) = c("oxEDA2R.vs.GFP")

#Reduce dataset ------------------------------
dds.sub = dds[which(!is.na(gene_annot[rownames(dds),]$SYMBOL)),]
temp = rowMeans(assay(dds.vst)[rownames(dds.sub),])
temp = temp[order(temp, decreasing = T)]
keep = names(temp[which(!duplicated(gene_annot[names(temp),]$SYMBOL))])
dds.sub = dds.sub[keep,]

# Perform DExpression ------------------------------
comlist = list()
for(i in 1:length(contlist)){
  print(i)
  temp = results(dds.sub, contrast =  contlist[[i]], independentFiltering = TRUE, alpha = 0.1, parallel = T, BPPARAM = MulticoreParam(4))
  temp = data.frame(gene_annot[rownames(temp),],temp, stringsAsFactors = F)
  temp = temp[order(temp$pvalue),]
  comlist[[i]] = temp
}
names(comlist) = names(contlist)

# load genes sets ------------------------------
library(msigdbr)
# Hallmarks
h_gene_sets = msigdbr(species = "Mus musculus", category = "H")
H <- list()
for(i in 1:length(table(h_gene_sets$gs_name))){
  print(i)
  current = names(table(h_gene_sets$gs_name))[i]
  H[[i]] = h_gene_sets[which(h_gene_sets$gs_name==current),]$gene_symbol
}
names(H) = names(table(h_gene_sets$gs_name))

# KEGG
kegg_gene_sets = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
K <- list()
for(i in 1:length(table(kegg_gene_sets$gs_name))){
  print(i)
  current = names(table(kegg_gene_sets$gs_name))[i]
  K[[i]] = kegg_gene_sets[which(kegg_gene_sets$gs_name==current),]$gene_symbol
}
names(K) = names(table(kegg_gene_sets$gs_name))

# REACTOME
R.Hs = qusage::read.gmt("~/Downloads/ReactomePathways.gmt")
  humGenes = unique(paste(unlist(R.Hs)))
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  genesV3 <- genesV2[which(!duplicated(genesV2[, 1])),]
  rownames(genesV3) = genesV3$HGNC.symbol
  # Print the first 6 genes found to the screen
  return(genesV3)
}
reactome.annot <- convert_human_to_mouse(humGenes)
rownames(reactome.annot) = reactome.annot[,1]
R.mm = R.Hs
for(i in 1:length(R.Hs)){
  print(i)
  temp = length(reactome.annot[rownames(reactome.annot) %in% R.Hs[[i]],])
  if(temp!=2){
    R.mm[[i]] = reactome.annot[rownames(reactome.annot) %in% R.Hs[[i]],][,2]
  }
  if(temp==2){
    R.mm[[i]] = reactome.annot[rownames(reactome.annot) %in% R.Hs[[i]],][2]
  }
}
R = R.mm

convert_mouse_to_human <- function(gene_list) { 
  output = c()
  mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  for(gene in gene_list) {
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "mouse, laboratory"))[['DB.Class.Key']]
    if( !identical(class_key, integer(0)) ) {
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes) {
        output = rbind(c(gene, human_gene), output)
      }
    }
  }
  return (output)
}

convert_human_to_mouse <- function(gene_list) {
  output = c()
  mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  for(gene in gene_list) {
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "human"))[['DB.Class.Key']]
    if( !identical(class_key, integer(0)) ) {
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(human_gene in human_genes) {
        output = rbind(c(gene, human_gene), output)
      }
    }
  }
  return (output)
}


# Import Gset order ------------------------------
hallmarks.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/hallmarks.R"
kegg_metabolism.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/kegg_metabolism.R"
kegg_gip.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/kegg_gip.R"
kegg_cprocess.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/kegg_cprocess.R"
kegg_orgsys.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/kegg_orgsys.R"
reactome_cellcycle.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_cellcycle.R"
reactome_cellcom.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_cellcomm.R"
reactome_develop.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_develop.R"
reactome_signaling.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_signaling.R"
reactome_damage.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_damage.R"
reactome_transcription.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_transcription.R"
reactome_transcriptionTargets.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_transcription_targets.R"
reactome_metab_carbo.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_metab_carbo.R"
reactome_metab_lipids.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_metab_lipids.R"
reactome_metab_energy.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_metab_energy.R"
reactome_metab_other.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_metab_other.R"

reactome_immune_innate.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_immune_innate.R"
reactome_immune_adaptive.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_immune_adaptive.R"
reactome_immune_cytokines.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_immune_cytokines.R"
reactome_apoptosis.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_apoptosis.R"
reactome_protlocalization.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_protlocalization.R"
reactome_vmtransport.sorted = "/Users/marcobolis/Downloads/CUSTOMSET/reactome_vmtransport.R"
reactome_sig_rtk1.sorted =  "/Users/marcobolis/Downloads/CUSTOMSET/reactome_sig_rtk1.R"
reactome_sig_rtk2.sorted =  "/Users/marcobolis/Downloads/CUSTOMSET/reactome_sig_rtk2.R"
reactome_sig_other1.sorted =  "/Users/marcobolis/Downloads/CUSTOMSET/reactome_sig_other1.R"
reactome_sig_other2.sorted =  "/Users/marcobolis/Downloads/CUSTOMSET/reactome_sig_other2.R"

CS = c(H, K)
inflammatory_custom = 
 c(
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_ALLOGRAFT_REJECTION",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_COAGULATION",
    "Complement cascade",
    "Initial triggering of complement",
    "Creation of C4 and C2 activators",
    "Alternative complement activation",
    "Activation of C3 and C5",
    "Terminal pathway of complement",
    "Regulation of Complement cascade"
)



# Import custom functions ------------------------------
source("/Users/marcobolis/Downloads/do_PlotGSEAdot.R")
source("/Users/marcobolis/Downloads/do_PlotGSEAdot2.R")

# Plot DE ------------------------------
vlist = c(2.5)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = H, gene_annot = gene_annot, orderlist = hallmarks.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = K, gene_annot = gene_annot, orderlist = kegg_metabolism.sorted, line = vlist)
#p1 <- do_PlotGSEAdot2(comlist = comlist, GS = K, gene_annot = gene_annot, orderlist = kegg_gip.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = K, gene_annot = gene_annot, orderlist = kegg_cprocess.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = K, gene_annot = gene_annot, orderlist = kegg_orgsys.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_cellcom.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_develop.sorted, line = vlist)

p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_sig_rtk1.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_sig_rtk2.sorted, line = vlist)

p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_sig_other1.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_sig_other2.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_signaling.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_damage.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_apoptosis.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_transcription.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_transcriptionTargets.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_metab_carbo.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_metab_lipids.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_metab_energy.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_metab_other.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_immune_innate.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_immune_adaptive.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_immune_cytokines.sorted, line = vlist)

p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_protlocalization.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comlist, GS = R, gene_annot = gene_annot, orderlist = reactome_vmtransport.sorted, line = vlist)
p1

CS = c(H,R)
p1 <- do_PlotGSEAdot2(comlist = comtemp, GS = CS, gene_annot = gene_annot, orderlist = custom, line = vlist)
p1

t1 <- do_PlotGSEAdot2(comlist = comtemp, GS = CS, gene_annot = gene_annot, orderlist = custom, line = vlist, returnText = T)

X1 = X # muscle human

comlist$oxEDA2R.vs.GFP$stat = res.orig[rownames(comlist$oxEDA2R.vs.GFP),]$stat
comlist$oxEDA2R.vs.GFP$stat = res.ashr[rownames(comlist$oxEDA2R.vs.GFP),]$log2FoldChange

t2 <- do_PlotGSEAdot3(comlist = comlist, GS = CS, gene_annot = gene_annot, orderlist = custom, line = vlist, returnText = F)

X2 = X

X = rbind(X1,X2)

P <- ggplot(X, aes(x = GROUP, y = fct_reorder(GS, ORDER))) + 
  geom_point(na.rm = T, aes(size = SCORE, color = I(COLOR))) +
  theme_bw(base_size = 8) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  #scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), limits=c(-14, 14), oob = scales::squish()) +
  ylab(NULL) +
  xlab(NULL) + 
  ggtitle("Hallmark pathway enrichment") +
  facet_grid(GS_GROUP ~ ., scale="free_y", space = "free_y") +
  theme(strip.text.y = element_text(angle = 0)) 
if(!is.null(line)){
  for(i in 1:length(line)){
    xi = line[i]
    P = P + geom_vline(xintercept = xi, linetype="dashed", color = "black", size=0.25)
  }
}

V = list() # top 500 up e top 500 down
V$age_up = agesig_up
V$age_dn = agesig_dn

convert

temp_eda2r = comlist$oxEDA2R.vs.GFP
temp_eda2r = temp_eda2r[which(!is.na(temp_eda2r$SYMBOL)),]
hsconv = convert_mouse_to_human(temp_eda2r$SYMBOL)
rownames(hsconv) = hsconv[,1]
temp_eda2r = temp_eda2r[which(temp_eda2r$SYMBOL %in% rownames(hsconv)),]
temp_eda2r$ortho = hsconv[temp_eda2r$SYMBOL,2]
temp_eda2r = temp_eda2r[which(!is.na(temp_eda2r$padj)),]
temp_eda2r = temp_eda2r[which(!is.na(temp_eda2r$ortho)),]

x = temp_eda2r$stat
names(x) = temp_eda2r$ortho
x = x[which(!is.na(x))]
idx.V <- ids2indices(V,id=names(x))
v_res = cameraPR.default(statistic = x, idx.V, use.ranks = F)


ggsave('plot.svg', P, width = 4, height = 3)
readtext("plot.svg")-> temp
write_clip(temp$text, object_type = "character")


res = comlist$oxEDA2R.vs.GFP

res.orig = results(dds.sub, contrast =  contlist[[1]], independentFiltering = TRUE, alpha = 0.1, parallel = T, BPPARAM = MulticoreParam(4))
res.ashr <- lfcShrink(res = res.orig, dds = dds.sub,  contrast = setNames(object = c(1, -1, 0, 0), resultsNames(dds.sub)), type="ashr")
res.apeglm <- lfcShrink(res = res.orig, dds = dds.sub,  contrast = setNames(object = c(1, -1, 0, 0), resultsNames(dds.sub)), type="apeglm")


temp = data.frame(gene_annot[rownames(res.ashr),],res.ashr, stringsAsFactors = F)
temp = temp[order(temp$pvalue),]


lab_italics <- paste0("italic('", rownames(res), "')")
selectLab_italics = paste0(
"italic('",
c('TET1','TET3','DNMT3A', 'DNMT3B','TDG','KDM1A','HDAC2','MEST','MYCN', 'MYBL'),
"')")

res = temp
res_noNans = temp[which(!is.na(temp$padj)),]
lab_italics <- paste0("italic('", rownames(res_noNans), "')")
p1 = EnhancedVolcano(res_noNans,xlim = c(-2,2), 
  lab = res_noNans$SYMBOL,
  x = 'log2FoldChange',
  y = 'padj',
  #selectLab = selectLab_italics,
  xlab = bquote(~Log[2]~ 'fold change'),
  pCutoff = 0.05,
  FCcutoff = 0.8,
  pointSize = 2,
  labSize = 4.0,
  labCol = 'black',
  #labFace = 'bold',
  #boxedLabels = TRUE,
  #parseLabels = TRUE,
  max.overlaps = 100,
  col = c('black', 'pink', 'purple', 'red3'),
  colAlpha = 5/5,
  legendPosition = 'bottom',
  legendLabSize = 14,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 1.0,
  colConnectors = 'black') #+ coord_flip()

ggsave('plot.svg', p1, width = 6, height = 6)
readtext("plot.svg")-> temp
write_clip(temp$text, object_type = "character")




# dotplot

library(ggplot2)

# Your provided matrix
dp <- matrix(c(5.7839607, 0.845474, 0.2045163, 6.250684, 10.8092191, 35.228787,
               0.5799195, 2.248266, 0.1010544, 4.597957, 13.4678749, 14.353339),
             nrow = 6, byrow = TRUE,
             dimnames = list(c("MET", "SEX", "AGE", "BMI", "SKING", "CRP"),
                             c("NTR", "NESDA")))
dp = dp[rev(c("SEX", "AGE", "BMI","SKING", "MET","CRP")),]

# Convert the matrix to a data frame
dp_df <- as.data.frame(as.table(dp))
colnames(dp_df) <- c("Predictor", "Dataset", "Value")


# Create the dot plot
p1 = ggplot(dp_df, aes(x = Dataset, y = Predictor, size = Value, fill = Value > 10)) +
    geom_point(shape = 21, alpha = 0.7) +
    scale_size_continuous(range = c(5, 15)) +
    scale_fill_manual(values = c("grey", "red")) +
    theme_minimal() +
    labs(x = "Dataset", y = "Predictor", size = "-10*log10(pvalue)")
ggsave('plot.svg', p1, width = 4, height = 3)
readtext("plot.svg")-> temp
write_clip(temp$text, object_type = "character")


# INFLAMMATION SPECIES
hsa = read.table("~/Desktop/REVX/hsa_muscle_corr.txt", sep = "\t", header = T, row.names = 1)
mmu = read.table("~/Desktop/REVX/mmu_muscle_corr.txt", sep = "\t", header = T, row.names = 1)
rno = read.table("~/Desktop/REVX/rno_muscle_corr.txt", sep = "\t", header = T, row.names = 1)
eda = res

mmu = mmu[order(abs(mmu$MUSCLE), decreasing = T),]
mmu$rank = seq(1:nrow(mmu))
colnames(mmu) = c("SYMBOL", "log2FoldChange", "pvalue")
mmu = mmu[which(!is.na(mmu$SYMBOL)),]
mmu$stat = mmu$log2FoldChange
mmu$padj = mmu$pvalue

hsa = hsa[order(abs(hsa$MUSCLE), decreasing = T),]
hsa$rank = seq(1:nrow(hsa))
colnames(hsa) = c("SYMBOL", "log2FoldChange", "pvalue")
hsa = hsa[which(!is.na(hsa$SYMBOL)),]
hsa$stat = hsa$log2FoldChange
hsa$padj = hsa$pvalue

rno = rno[order(abs(rno$MUSCLE), decreasing = T),]
rno$rank = seq(1:nrow(rno))

colnames(rno) = c("SYMBOL", "log2FoldChange", "pvalue")
rno = rno[which(!is.na(rno$SYMBOL)),]
rno$stat = rno$log2FoldChange
rno$padj = rno$pvalue

head(comlist$oxEDA2R.vs.GFP)

comtemp = list()
comtemp$hsa = hsa

comtemp$mmu = mmu

comtemp$rno = rno

p1 <- do_PlotGSEAdot2(comlist = comtemp, GS = H, gene_annot = gene_annot, orderlist = hallmarks.sorted, line = vlist)
p1 <- do_PlotGSEAdot2(comlist = comtemp, GS = R, gene_annot = gene_annot, orderlist = hallmarks.sorted, line = vlist)







## Create merged Seurat Object ------------------------------------------------------------------------
for(i in 1:length(samples)){
  print(i)
  currName = names(samples)[i]
  if(i==1){seed = CreateSeuratObject(as.data.frame(samples[[i]]@assays$RNA@counts), project = currName)}
  if(i==2){
    rest = list()
    rest[[i-1]] = CreateSeuratObject(as.data.frame(samples[[i]]@assays$RNA@counts), project = currName)
  }
  if(i>2){rest[[i-1]] = CreateSeuratObject(as.data.frame(samples[[i]]@assays$RNA@counts), project = currName)}
}
merged_seurat <- merge(x = seed, y = rest, add.cell.ids = names(samples) )

  

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = merged_seurat, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 10

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

merged_seurat@meta.data$BIOREP = rep(c("B1","B1","B1","B2","B2","B2","B3","B3","B3"),2)


# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = merged_seurat@meta.data)

# Normalize the counts
filtered_seurat <- NormalizeData(filtered_seurat)

# Identify the most variable genes
filtered_seurat <- FindVariableFeatures(filtered_seurat, 
                                        selection.method = "vst",
                                        nfeatures = 3000, 
                                        verbose = FALSE)

# Scale the counts
filtered_seurat <- ScaleData(filtered_seurat, vars.to.regress = "BIOREP")

# Perform PCA
filtered_seurat <- RunPCA(filtered_seurat, npcs = 16)
filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:10, n.neighbors = 15L)

# Plot UMAP                             
DimPlot(filtered_seurat)

DimPlot(filtered_seurat, reduction = "pca",  group.by = "orig.ident")

DimPlot(filtered_seurat, split.by = "orig.ident", group.by = "orig.ident")  

# Plot the elbow plot
ElbowPlot(object = filtered_seurat, ndims = 40)

# Determine percent of variation associated with each PC
pct <- filtered_seurat[["pca"]]@stdev / sum(filtered_seurat[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs

# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

# Determine the K-nearest neighbor graph
filtered_seurat <- FindNeighbors(object = filtered_seurat, dims = 1:pcs)
filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:pcs)

# Determine the clusters for various resolutions                                
filtered_seurat <- FindClusters(object = filtered_seurat, resolution = seq(from = 0.1, to = 1, by = 0.1))

clustree(filtered_seurat, prefix = "RNA_snn_res.") #stop

filtered_seurat@meta.data = cbind(filtered_seurat@meta.data, mannot[rownames(filtered_seurat@meta.data),])

DimPlot(filtered_seurat, reduction = "umap", pt.size = 0.5, label=TRUE, group.by = "RNA_snn_res.0.9") # A MANO
FeaturePlot(filtered_seurat, feature = "GROUP",  reduction = "umap" );p.umap 
DimPlot(filtered_seurat, reduction = "umap", pt.size = 0.5, label=TRUE, group.by = "GROUP") # A MANO



