setwd("/Users/marcobolis/Desktop/Projects/INFLAMMAGING/EDA2R-BRIEF/Progeria/")
rawcounts = read.table("./GSE165409_LMNA_rawCounts_18samp.txt", header = T, row.names = 1, sep = "\t")

plot(colSums(rawcounts))

# Load libraries ------------------------------
neededlibraries <- c("edgeR", "RColorBrewer", "forcats", "parallel", "doParallel", "BiocParallel", "DESeq2", "ggplot2", "ggrepel","org.Hs.eg.db", "GSVA", "monocle", "qusage", "ggplot2", "clipr", "readtext", "scales")
lapply(neededlibraries, require, character.only = TRUE)

# Create annotation file ------------------------------
write_clip(colnames(rawcounts))
annotations = data.frame(row.names = colnames(rawcounts), GROUP = c(rep("LMNA_Y",6), rep("WT_Y",6), rep("WT_O",6)), OUTLIER = c(rep("NO",12),rep("YES",6)))


# Load dataset ------------------------------
SPECIES = "MOUSE"
if(SPECIES=="MOUSE"){initial = "ENSM"}
if(SPECIES=="HUMAN"){initial = "ENSG"}


# Remove outliers ------------------------------
#annotations = annotations[which(annotations$OUTLIER=="NO"),]
rawcounts = rawcounts[,rownames(annotations)]

# retrieve gene annotations ------------------------------
if(SPECIES=="MOUSE"){
  library(org.Mm.eg.db)
  Db = org.Mm.eg.db
}
if(SPECIES=="HUMAN"){
  library(org.Hs.eg.db)
  Db = org.Hs.eg.db
}
symbols <- mapIds(Db,
                  keys = row.names(rawcounts),
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first"
)
entrez <- mapIds(Db,
                 keys = row.names(rawcounts),
                 column = "ENTREZID",
                 keytype = "ENSEMBL",
                 multiVals = "first")
gene_annot = data.frame(SYMBOL=symbols, ENTREZ=entrez)
rownames(gene_annot) <- rownames(rawcounts)

# Create DESeq2 Object ------------------------------
samples = annotations
dds <- DESeqDataSetFromMatrix(countData = rawcounts, colData = samples, design = ~ 0 + GROUP)
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
colTips= c( rep("white", ncol(dds.assay)  ) )
names(colTips) = colnames(dds.assay)
colTips[rownames(samples)[which(samples$GROUP=="LMNA_Y")]] = "blue"
colTips[rownames(samples)[which(samples$GROUP=="WT_Y")]] = "darkgreen"
colTips[rownames(samples)[which(samples$GROUP=="WT_O")]] = "red"
plot(as.phylo(hc), type = "phylogram", cex = 0.8, edge.color = "black", edge.width = 2, edge.lty = 1, tip.color = colTips)
# save 5x7

# Define design ------------------------------
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, parallel = T, BPPARAM = MulticoreParam(4))
dds.design = model.matrix(~ 0 + GROUP, data=samples)
colnames(dds.design) = gsub(x = colnames(dds.design), pattern = "GROUP", replacement = "")

# Define Constrasts ------------------------------
contlist = list()
contlist[[1]] =  makeContrasts(LMNA_Y-WT_Y,levels = dds.design)
contlist[[2]] =  makeContrasts(WT_O-WT_Y,levels = dds.design)
contlist[[3]] =  makeContrasts(LMNA_Y-WT_O,levels = dds.design)
names(contlist) = c("LMNA_Y.vs.WT_Y","WT_O.vs.WT_Y","LMNA_Y.vs.WT_O")

# Reduce dataset ------------------------------
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
reactome.annot <- convertHumanGeneList(humGenes)
R.mm = R.Hs
for(i in 1:length(R.Hs)){
  print(i)
  R.mm[[i]] = reactome.annot[rownames(reactome.annot) %in% R.Hs[[i]],]$MGI.symbol
}
R = R.mm

# Import Gset order ------------------------------
hallmarks.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/hallmarks.R"
kegg_metabolism.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/kegg_metabolism.R"
kegg_gip.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/kegg_gip.R"
kegg_cprocess.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/kegg_cprocess.R"
kegg_orgsys.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/kegg_orgsys.R"
reactome_cellcycle.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_cellcycle.R"
reactome_cellcom.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_cellcomm.R"
reactome_develop.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_develop.R"
reactome_signaling.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_signaling.R"
reactome_damage.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_damage.R"
reactome_transcription.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_transcription.R"
reactome_transcriptionTargets.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_transcription_targets.R"
reactome_metab_carbo.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_metab_carbo.R"
reactome_metab_lipids.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_metab_lipids.R"
reactome_metab_energy.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_metab_energy.R"
reactome_metab_other.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_metab_other.R"

reactome_immune_innate.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_immune_innate.R"
reactome_immune_adaptive.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_immune_adaptive.R"
reactome_immune_cytokines.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_immune_cytokines.R"
reactome_apoptosis.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_apoptosis.R"
reactome_protlocalization.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_protlocalization.R"
reactome_vmtransport.sorted = "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_vmtransport.R"
reactome_sig_rtk1.sorted =  "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_sig_rtk1.R"
reactome_sig_rtk2.sorted =  "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_sig_rtk2.R"
reactome_sig_other1.sorted =  "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_sig_other1.R"
reactome_sig_other2.sorted =  "/Users/marcobolis/Desktop/Projects/REPO/CUSTOMSET/reactome_sig_other2.R"


# Import custom functions ------------------------------
source("/Users/marcobolis/Desktop/Projects/REPO/do_PlotGSEAdot.R")
source("/Users/marcobolis/Desktop/Projects/REPO/do_PlotGSEAdot2.R")

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

ggsave('plot.svg', p1, width = 6.5, height = 9.7)
readtext("plot.svg")-> temp
write_clip(temp$text, object_type = "character")

# Counts per Million plot 
cpmx = dds.cpm["ENSMUSG00000034457",1:12]
cpmbplot = data.frame(GROUP = c(rep("HGPS",6),rep("WT",6)), EXP = c(as.numeric(as.character(dds.cpm["ENSMUSG00000034457",1:12]))))
cpmbplot$GROUP = factor(cpmbplot$GROUP, levels = c("WT", "HGPS"))

p_gene = ggplot(cpmbplot, aes(x=GROUP, y=EXP, fill=GROUP)) +
  #stat_compare_means(comparisons = my_comparisons) +
  #scale_shape_manual(name = "df", values = c(21)) +
  #scale_y_continuous(limits = c(40,180)) + 
  geom_violin(aes(fill=GROUP), width = 0.4, trim = F) +
  geom_jitter( aes( fill=GROUP), shape = 21,  size=0.5, alpha=1,  position = position_jitter(width = .04)) +
  #stat_compare_means(comparisons = my_comparisons,  color = "black", bracket.size = 0.5, ) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=-90,vjust = 0))
p_gene

ggsave('plot.svg', p_gene, width = 4, height = 3.8)
readtext("plot.svg")-> temp
write_clip(temp$text, object_type = "character")

# VolcanoPlot


res = comlist$LMNA_Y.vs.WT_Y
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


