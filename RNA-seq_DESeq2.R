######  title: "RNA-seq DESeq2"
######  author: Feng_yuanqing

library("DESeq2")
library("tidyverse")

# Here I use OCA2 CRISPRi as an example. 

#---------------------DESeq2 prep Counts----------------------#
# estimation of size factors: estimateSizeFactors
# estimation of dispersion: estimateDispersions
# Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest

# read in data (RNA read count)
data <- read.table(file.choose(),header = TRUE, row.names =  "id")
# change digits to integers, remove gene with mean count < 10
data2 <-(round(data)) %>% dplyr::filter(rowMeans(.)>10)
head(data2)
# set sample groups,Assign conditions
mysample<-gsub(".sg1sg2.[1-3]$|5.[1-3]", "", colnames(data2))
# condition
condition <- factor(mysample)
# Create a coldata frame and the DESeqDataSet
coldata <- data.frame(row.names=colnames(data2), condition)
dds <- DESeqDataSetFromMatrix(countData=data2, colData=coldata, design=~condition)
# dds
## perform differential calculation
dds2 <- DESeq(dds)
## Check the size factors
sizeFactors(dds2)

#——------------------- normalize read counts between groups -----------------
# normalized count with size factors, used for between sample comparison
# counts divided by sample-specific size factors determined by median Counts of gene counts relative to geometric mean per gene # NOTE: DESeq2 doesn’t actually use normalized counts, rather it uses the raw counts and models the normalization inside the Generalized Linear Model (GLM). These normalized counts will be useful for downstream visualization of results, but cannot be used as input to DESeq2 or any other tools that peform differential expression analysis which use the negative binomial model.
normalized_count <- counts(dds2, normalized=T) %>% 
  data.frame() %>% 
  round(digits = 2) %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble() 
write.table(normalized_count, file = "all_kallisto_normalized.count.tsv", row.names = F, sep = "\t", quote = F)

## Total number of normalized counts per sample
# colSums(normalized_count)
## Total number of raw counts per sample
# colSums(counts(dds2))

# Plot dispersion
png("DESeq2_all_dispersion.png", 600, 600, pointsize=20)
plotDispEsts(dds2, main="Dispersion_plot_all_samples")
dev.off()

#——-------------- DEG between groups -----------------
# colnames(normalized_count)
# Set condition
res1 <- results(dds2, contrast=c("condition","OCA2","plko"))
## MA Plot to show fold change. An MA-plot is a plot of log-intensity Countss (M-values) versus log-intensity averages (A-values). 
png("DESeq2_OCA2_CRISPRi_MAplot.png", 600, 600, pointsize=20)
plotMA(res1,  main="OCA2_enhancer_inhibition", ylim=c(-6,6),alpha = 0.05)
dev.off()
## order the data by p-value
resOrdered <- res1[order(res1$pvalue),] %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
head(resOrdered)
#summary(resOrdered)
write.table(resOrdered, file = "DESeq2_OCA2_CRISPRi_Pvalue.tsv", row.names = F, sep = "\t", quote = F)

#-------------------select significant genes ---------------
sum(resOrdered$padj < 0.1, na.rm=TRUE)
siggenes <- resOrdered  %>% filter(padj < 0.1 & abs(log2FoldChange) > 0.2)
# resOrdered[row.names(resOrdered)== "OCA2",] # p-value
# filter(normalized_count, gene=="OCA2") # normalized counts
Sig_gene_count <- normalized_count %>% filter(normalized_count$gene %in% siggenes$gene)
# reads count of significant genes
write.table(Sig_gene_count, file = "DESeq2_OCA2_CRISPRi_Siggene_counts.tsv", row.names = F, sep = "\t", quote = F)

#---------------- heatmap of selected genes---------------
#BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")
library("RColorBrewer")
# select top20 gene 
top20genes<-head(siggenes$gene,20)
top20count_matrix <- normalized_count %>% dplyr::slice(match(top20genes, normalized_count$gene)) %>% column_to_rownames( var = "gene") %>% dplyr::select(matches("OCA2|plko5")) %>% as.matrix() 
#  plot heatmaps of log(counts)
annotation_col = data.frame(Group = factor(c(rep("OCA2_sgRNA", 3),rep("Control_sgRNA", 3))))
ann_colors = list(Group = c(OCA2_sgRNA = "#984ea3", Control_sgRNA = "#a6d854"))
pdf("DESeq2_OCA2_CRISPRi_sigGene_top20_heatmap.pdf", width = 8, height = 8)
ht <- ComplexHeatmap::pheatmap(log10(top20count_matrix),annotation_col = annotation_col, annotation_colors = ann_colors, 
column_split = annotation_col$Group,cluster_cols = T, cluster_rows = F,heatmap_legend_param = list(title = "log10(count)"),main="Normalized counts of top 20 genes",scale="none",cellwidth=20,cellheight=15,angle_col="45",color = colorRampPalette(c("#4393c3", "#d1e5f0","#fddbc7","#d6604d"))(10) )
dev.off();ht

# ----------------------------------------Enhanced Volcano plot using a package-------------------------------- ---------- 
library("EnhancedVolcano")
library("scales")
# replace p=0 with 1e-100
resOrdered1<- resOrdered %>% mutate(padj = replace(padj, padj<1e-100, 1e-50)) %>% filter(pvalue<0.05) %>% na.omit()
# gene to highlight
mysig_gene<-head(resOrdered1$gene,10)
# set lim 
xlim=c(min(resOrdered1$log2FoldChange)-0.1,max(resOrdered1$log2FoldChange)+0.1)
ylim=c(0, max(-log10(resOrdered1$padj)) + 1)
# plot 
pv<-EnhancedVolcano(resOrdered1,lab = resOrdered1$gene, selectLab=mysig_gene, xlim=xlim,ylim=ylim, x = 'log2FoldChange', y = 'padj', ylab = bquote(~-Log[10] ~ italic(P)), pointSize = 2,labSize =4,title = "DESeq2_OCA2_EN_CRISPRi (MNT1)",subtitle = NULL,titleLabSize = 14, legendPosition = "top",legendLabSize = 14,cutoffLineWidth = 0.2, borderWidth=0.4, col = c('grey30', 'forestgreen', 'royalblue', 'red2'),colAlpha = 0.7,drawConnectors = TRUE,widthConnectors = 0.5, pCutoff = 0.01, FCcutoff = 0.5,caption = 'log2FoldChange cutoff, 0.5; padj cutoff, 0.01');pv
ggsave("DESeq2_OCA2_EnCRISPRi_Vplot.pdf",plot = pv, width =12, height = 16, units = "cm")


#---------------------Volcano plot of DEG --------------------------------
Vplot<-function(x){
  x$logp <- -log10(x$padj)
  snp_loc <- rownames(x)
  rownames(x) <- seq(1:length(x$log2FoldChange))
  siggene_colors= factor(ifelse(x$logp>3, 2, 1))
  mylabelgene<-ifelse(x$logp>10,snp_loc,'') # change the logp cutoff and xlim, ylim
  p<-ggplot(data = x, aes(x = log2FoldChange, y = logp)) + 
    geom_point(aes(color=siggene_colors),size=1,alpha=1) + 
    theme_bw() + labs(x ="Log2(Fold Change)", y="Log10(Adjust P)") + 
    xlim(-4, 4) + ylim(0, 40) + 
    theme(legend.title=NULL, legend.text=element_text(size=12), axis.text=element_text(size=12), 
          axis.title=element_text(size=12.5), plot.margin=margin(t=20, r=20, b=20, l=20, "pt")) + 
    geom_text_repel(aes(label=mylabelgene),box.padding = 0.5, max.overlaps = Inf,size=3, show.legend= F) +
    scale_colour_manual(values = c("grey80", "red"),labels = c("NS", "Log10(Adjust P)>3")) +
    geom_abline(intercept = 3, slope = 0, linetype = 2,color="grey") 
  p## could change the xlim and ylim scale to avoid lost extrem values.
  return(p)
}
# read in the dataframe from DESeq2
mysig_gene<-siggenes%>% column_to_rownames( var = "gene") %>% as.data.frame() 
p<-Vplot(mysig_gene);p
ggsave("DESeq2_OCA2_CRISPRi_sigGene_Vplot.pdf",plot = p, width =16, height = 16, units = "cm")
#ggplotly(p)

# ------ Volcano plot using a package ---------- 
#library("EnhancedVolcano")
#EnhancedVolcano(mysig_gene,lab=rownames(mysig_gene), x="log2FoldChange", y="padj", title = NULL,titleLabSize = 14, subtitle=NULL,axisLabSize=14, pCutoff = 10e-4, FCcutoff =0.5, xlim = c(-4, 4), ylim = c(0, 40), colAlpha = 0.8,shape = 19,legendLabSize =10, captionLabSize = 10, cutoffLineType = 'longdash', cutoffLineWidth = 0.2, drawConnectors = F)
#ggsave("Vocano.pdf", width =12.5, height = 15, units = "cm")

#----------------------- Scatter plot----------------------
Scatt_plot <- function(Counts){
  Ref=rowMeans(Counts[,1:3])  # Case
  Alt=rowMeans(Counts[,4:6])  # Control
  LFC=as.data.frame(Alt-Ref) # log2 Fold change
  Name=rownames(LFC)
  mylables<-ifelse(abs(LFC)>1,Name,'') # could modify lfc cutoff
  Counts= mutate(Counts,sig=ifelse(abs(LFC)>1,2,1)) # could modify lfc cutoff
  p <- ggplot(Counts, aes(x=Ref, y=Alt,color=factor(sig))) +
    theme_bw() + geom_point(size=1,alpha=1) +
    labs(y ="Log2(Count) - Control", x="Log2(Count) - OCA2 CRISPRi") + 
    geom_abline(intercept = 0, slope = 1, linetype = 2, color="grey") + coord_fixed() + 
    theme(legend.title=NULL, legend.text=element_text(size=12), axis.text=element_text(size=12), 
          axis.title=element_text(size=12.5),plot.margin=margin(t=20, r=20, b=20, l=20, "pt")) + 
    geom_text_repel(aes(label=mylables),box.padding = 0.5, max.overlaps = Inf,size=3, show.legend = F )
  p
  return(p)
}
# Select read counts of significant genes
sCounts<-Sig_gene_count %>% column_to_rownames(var="gene") %>%  dplyr::select(matches("OCA2|plko5")) %>% log2() %>% as.data.frame()
# scatter plot and label sig genes
p<-Scatt_plot(sCounts)+scale_color_manual(values = c("grey80", "red"),labels = c("NS", "Log2|FC|>1"));p
ggsave("DESeq2_OCA2_CRISPRi_sigGene_Splot.pdf", width =16, height = 16, units = "cm")


##------------ correlation plot between all samples --------------------
library("debrowser")
# modified from debrowser::all2all(sCounts)
allcor <- function (data, cex=2) {
  pcor2 <- function(x, y, ...) debrowser::panel.cor(x, y, cex.cor = cex, prefix = "")
  nr <- nrow(data)
  if (nr > 2000) 
    nr <- 2000 # only plot 2000 genes
  pairs(log2(data[1:nr, ]), cex = 0.25, diag.panel = panel.hist, lower.panel = pcor2,main = "Log2(Normalized counts)")
}
## Calculate correlation of normalized counts
allcor(sCounts) 
# save the picture as it is in Rstudio
dev.print(pdf, "DESeq2_OCA2_CRISPRi_sigGene_corr.pdf") 
dev.off()

