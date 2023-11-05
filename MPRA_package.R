# -----------------MPRA Data analysis-----------------------
# Feng_yuanqing
# This script is used to compare the allele activity using MPRA barcode counts data from DNA and RNA.

#1. Install mpra and its dependencies by running the lines below in R or RStudio.
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("mpra")
library("mpra")
library("debrowser")
library("tidyverse")
library("ggrepel")
library("ggpubr")
library("ggExtra")
library("reshape2")

## --------------------------------- -----read in a data frame-----------------------
## read in a data frame that have MPRA counts data, with col1:Oligo_name; col2: Sample1;  col3: Sample2; .... MPRA_DNA_RNA_barcode.count.tsv
data <- read.table(file.choose(),header =T)
head(data) 

## ------------------------extract Ref and Alt Allele separately------------------------
alt <- data %>% filter(str_detect(Enhancer_name,"^Alt"))
alt$oligo <- gsub("Alt_","",alt$Enhancer_name, perl=TRUE)
ref <- data %>% filter(str_detect(Enhancer_name,"^Ref")) 
ref$oligo <- gsub("Ref_","",ref$Enhancer_name, perl=TRUE)

## ------------------------extract SNPs have both Ref and Alt Allele------------------------
alt2 <- dplyr::filter(alt, alt$oligo %in% ref$oligo) 
ref2 <- dplyr::filter(ref, ref$oligo %in% alt$oligo) 
col_name <- c("ref_1","ref_2","ref_3", "alt_1","alt_2","alt_3")
WM88_DNA <- data.frame(cbind(ref2[,4:6], alt2[,4:6]),row.names = alt2$oligo); colnames(WM88_DNA) <- col_name
colSums(WM88_DNA, na.rm = TRUE)
WM88_RNA <- data.frame(cbind(ref2[,7:9], alt2[,7:9]),row.names = alt2$oligo); colnames(WM88_RNA) <- col_name
colSums(WM88_RNA, na.rm = TRUE)
MNT1_DNA <- data.frame(cbind(ref2[,10:12], alt2[,10:12]),row.names = alt2$oligo); colnames(MNT1_DNA) <- col_name
colSums(MNT1_DNA, na.rm = TRUE)
MNT1_RNA <- data.frame(cbind(ref2[,13:15], alt2[,13:15]),row.names = alt2$oligo); colnames(MNT1_RNA) <- col_name
colSums(MNT1_RNA, na.rm = TRUE)

##---------------------------- construct a MPRASet ------------------------
mpra_MNT <- MPRASet(DNA = MNT1_DNA, RNA = MNT1_RNA, eid = alt2$oligo, eseq = NULL, barcode = NULL)
mpra_WM88 <- MPRASet(DNA = WM88_DNA, RNA = WM88_RNA, eid = alt2$oligo, eseq = NULL, barcode = NULL)

##----------------------------- normalization of the barcode counts bwtween samples------------------
## Calculate the RNA/DNA ratio of normalized counts (both DNA and RNA are sceled to have 10 million counts )
Ratio_mnt1 <- getRNA(normalize_counts(mpra_MNT)) / getDNA(normalize_counts(mpra_MNT))
Ratio_wm88 <- getRNA(normalize_counts(mpra_WM88)) / getDNA(normalize_counts(mpra_WM88))
write.table(Ratio_mnt1, file = "MNT1_Ref_vs_Alt_ratio.txt", row.names = T, sep = "\t", quote = F)
write.table(Ratio_wm88, file = "WM88_Ref_vs_Alt_ratio.txt", row.names = T, sep = "\t", quote = F)
# Should normalized the DNA or RNA counts independently using limma 

#------------------------ MPRASet limma_voom_eBayes DEG Analysis for MNT1------------------------
detach("package:tidyverse", unload=TRUE)
mpra_MNT <- MPRASet(DNA = MNT1_DNA, RNA = MNT1_RNA, eid = alt2$oligo, eseq = NULL, barcode = NULL)
# the desgin of the MPRASet, 3 ref, 3 alt.
design <- data.frame(intcpt = 1, alleleB = grepl("alt", colnames(mpra_MNT)))
# block is a vector that is used when the columns of the MPRAset object are paired.
block_vector <- rep(1:3, 2)
# Linear models for differential analysis of MPRA data
mpralm_MNT_fit <- mpralm(object = mpra_MNT, design = design, aggregate = "none", normalize = TRUE, block = block_vector, model_type = "corr_groups", plot = TRUE)
toptab_MNT <- topTable(mpralm_MNT_fit, coef = 2, number = Inf)
head(toptab_MNT) 
toptab_MNT_ranked <- toptab_MNT %>% rownames_to_column() %>% arrange(adj.P.Val) 
## Write the results into a txt file
write.table(toptab_MNT_ranked, file = "MNT1_Ref_vs_Alt_mpralm_sigdiff.tsv", row.names = F, sep = "\t", quote = F)

#------------------------ MPRASet limma_voom_eBayes DEG Analysis for WM88 ------------------------
mpra_WM88 <- MPRASet(DNA = WM88_DNA, RNA = WM88_RNA, eid = alt2$oligo, eseq = NULL, barcode = NULL)
design <- data.frame(intcpt = 1, alleleB = grepl("alt", colnames(mpra_WM88)))
block_vector <- rep(1:3, 2)
mpralm_WM88_fit <- mpralm(object = mpra_WM88, design = design, aggregate = "none", normalize = TRUE, block = block_vector, model_type = "corr_groups", plot = TRUE) ## here using normalized counts, All columns are scaled to have 10 million counts.
toptab_WM88 <- topTable(mpralm_WM88_fit, coef = 2, number = Inf)
head(toptab_WM88)
toptab_WM88_ranked <- toptab_WM88 %>% rownames_to_column() %>% arrange(adj.P.Val) 
## Write the results into a txt file
write.table(toptab_WM88_ranked, file = "WM88_Ref_vs_Alt_mpralm_sigdiff.tsv", row.names = F, sep = "\t", quote = F)

# ----------------------------------------FC P-values plot using a package-------------------------------- ---------- 
library("EnhancedVolcano")
library("scales")
# wm88 select top snp to highlight
selectLab <- toptab_WM88 %>% filter(adj.P.Val<1e-20 & abs(logFC)>0.25 ) %>% rownames()
summary(-log10(toptab_WM88$adj.P.Val))
summary(toptab_WM88$logFC)
# plotting ggplot
p<-EnhancedVolcano(toptab_WM88_M, lab=rownames(toptab_WM88), selectLab=selectLab, x="logFC", y="adj.P.Val",hlineCol="steelblue",hline= c(1e-20), borderWidth = 0.25, labSize = 4, title = NULL,titleLabSize = 16, subtitle=NULL,axisLabSize=16, pCutoff =1e-2, FCcutoff =0.2, xlim = c(-1.2, 1.2), ylim = c(0, 35), pointSize=1, colAlpha = 0.8,shape = 19,legendLabSize =15, captionLabSize = 15, cutoffLineType = 'longdash', cutoffLineWidth = 0.25, drawConnectors = T)
p 
ggsave("Volcano_WM88.pdf", width =5, height = 6, units = "in")

# MNT1
selectLab2 <- toptab_MNT %>% filter(adj.P.Val<1e-20 & abs(logFC)>0.25 ) %>% rownames()
summary(-log10(toptab_MNT$adj.P.Val))
summary(toptab_MNT$logFC)
p2<-EnhancedVolcano(toptab_MNT, lab=rownames(toptab_MNT), selectLab=selectLab2, x="logFC", y="adj.P.Val",hlineCol="steelblue",hline= c(1e-20), borderWidth = 0.25, labSize = 4, title = NULL,titleLabSize = 16, subtitle=NULL,axisLabSize=16, pCutoff =1e-2, FCcutoff =0.2, xlim = c(-2.4, 1), ylim = c(0, 120), pointSize=1, colAlpha = 0.6,shape = 19,legendLabSize =15, captionLabSize = 15, cutoffLineType = 'longdash', cutoffLineWidth = 0.25, drawConnectors = T)
p2
ggsave("Volcano_MNT1.pdf", width =5, height = 6, units = "in")


#-------------- correlation of RNA/DNA ratio mnt and wm88 ------------ 
panel.hist <- function(x, ...){ 
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE,breaks = "FD" )
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan2", ...)}

allcor <- function (data, cex=2) {
  pcor2 <- function(x, y, ...) panel.cor(x, y, cex.cor = cex, prefix = "")
  nr <- nrow(data)
  if (nr > 2000) 
    nr <- 2000
  pairs(log10(data[1:nr, ]), cex = 0.25, diag.panel = panel.hist, lower.panel = pcor2)}
## Calculate the RNA/DNA ratio of normalized counts
Ratio_mnt1 <- getRNA(normalize_counts(mpra_MNT)) / getDNA(normalize_counts(mpra_MNT))
Ratio_wm88 <- getRNA(normalize_counts(mpra_WM88)) / getDNA(normalize_counts(mpra_WM88))
allcor(Ratio_mnt1)
allcor(Ratio_wm88)
colnames(Ratio_mnt1)<-c("MNT1_Ref1","MNT1_Ref2","MNT1_Ref3","MNT1_Alt1","MNT1_Alt2","MNT1_Alt3")
colnames(Ratio_wm88)<-c("WM88_Ref1","WM88_Ref2","WM88_Ref3","WM88_Alt1","WM88_Alt2","WM88_Alt3")
allratio<-cbind(Ratio_mnt1,Ratio_wm88)
allcor(allratio)

#------------------- PCA plot of RNA/DNA ratios ref.alt in mnt and wm88 ----------------
PCA_plot <- function(x) {
  pca_data <- debrowser::run_pca(x)
  # x is a read counts data.frame.
  p_data <-  debrowser::prepPCADat(pca_data)
  pcx=1;pcy=2
  xaxis <- sprintf("PC%d (%.2f%%)", pcx, round(pca_data$explained[pcx] * 100, 2))
  yaxis <- sprintf("PC%d (%.2f%%)", pcy, round(pca_data$explained[pcy] * 100, 2))
  plot1 <- ggplot(data = p_data, aes(x = x, y = y)) + geom_point(aes(color = color), size = 4) + labs(x = xaxis, y = yaxis) + theme( legend.text=element_text(size=12), axis.text=element_text(size=12), axis.title=element_text(size=14),plot.margin=margin(t=20, r=20, b=20, l=20, "pt")) + geom_text(aes(label = samples), vjust = 0, nudge_y = -1) + theme_light() 
  return(plot1)
}
PCA_plot(getDNA(normalize_counts(mpra_WM88)))
PCA_plot(Ratio_wm88)
PCA_plot(allratio)

##------------Compare Fold change and pvalues between two cell lines using 4 Quadrants --------------------
# sorting the data by oligo name
wm_rank <- toptab_WM88 %>% rownames_to_column() %>% arrange(rowname) 
mnt_rank <- toptab_MNT %>% rownames_to_column() %>% arrange(rowname) 

# build a factor for intersected significant snps and uniq ones. 
sig_mnt <- dplyr::filter(mnt_rank,adj.P.Val<0.05) 
sig_wm <- dplyr::filter(wm_rank,adj.P.Val<0.05)
both <- dplyr::filter(sig_mnt,sig_mnt$rowname %in% sig_wm$rowname) %>% pull(rowname)
mnt_only <- dplyr::filter(sig_mnt,!(sig_mnt$rowname %in% sig_wm$rowname)) %>% pull(rowname)
wm_only <- dplyr::filter(sig_wm,!(sig_wm$rowname %in% sig_mnt$rowname)) %>% pull(rowname)
snp_rank <-as.data.frame(wm_rank$rowname)
snp_rank= mutate(snp_rank,b =ifelse(snp_rank[,1] %in% both ,1,0))
snp_rank= mutate(snp_rank,mnt =ifelse(snp_rank[,1] %in% mnt_only ,2,0))
snp_rank= mutate(snp_rank,wm =ifelse(snp_rank[,1] %in% wm_only ,3,0))
snp_rank= mutate(snp_rank,sum = rowSums(snp_rank[,2:4]))
head(snp_rank)

# combine the fold change data from MNT and WM88, With factor info.
fc2 <- cbind.data.frame(mnt_rank$logFC,wm_rank$logFC,snp_rank$sum)
colnames(fc2)=c("MNT_log2FC","WM88_log2FC","factor")
# simple plot
p1<-ggscatter(fc2, x = "MNT_log2FC", y = "WM88_log2FC", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", main="MNT1_vs_WM88",font.main=c(12, "bold.italic", "blue4"),  xlab="MNT_log2FC ", ylab="WM88_log2FC", shape = 20, size=2, color ="dodgerblue4", cor.coeff.args = list(method = "pearson", label.x = -2, label.sep = "\n"))
p1 + geom_hline(yintercept=0,linetype = 2,color="grey") + geom_vline(xintercept=0,linetype = 2,color="grey") 
# plot with classification factors. 
# using summary(fc2) to determine the axis scale.
p3 <- ggplot(fc2, aes(x=MNT_log2FC, y=WM88_log2FC,color=factor(fc2$factor))) +
  theme_bw() + geom_point(size=0.5,alpha=1) + 
  scale_color_manual(values=c("grey","#1B9E77","#D95F02","#7570B3"),labels = c("NS","both","MNT_only","WM88_only"))+
  labs(x ="MNT_log2FC", y="WM88_log2FC") + 
  geom_hline(yintercept=0,linetype = 2,color="grey") + 
  geom_vline(xintercept=0,linetype = 2,color="grey") +
  xlim(-2.5, 1) + ylim(-0.75, 0.75) + 
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title=element_text(size=14),plot.margin=margin(t=20, r=20, b=20, l=20, "pt")) 
p3
## add the percentage at four Quadrants
library("scales")
total <- length(both)+length(mnt_only)+length(wm_only)
per1 <- filter(fc2,fc2$MNT_log2FC>0 &fc2$WM88_log2FC>0) %>% filter(factor!=0); pp<-length(per1$factor)/total
per2 <- filter(fc2,fc2$MNT_log2FC>0 &fc2$WM88_log2FC<0) %>% filter(factor!=0); pn<-length(per2$factor)/total
per3 <- filter(fc2,fc2$MNT_log2FC<0 &fc2$WM88_log2FC>0) %>% filter(factor!=0); np<-length(per3$factor)/total
per4 <- filter(fc2,fc2$MNT_log2FC<0 &fc2$WM88_log2FC<0) %>% filter(factor!=0); nn<-length(per4$factor)/total
p3+ annotate("text", x=0.75, y=0.75, label=percent(pp), size=5)+ 
  annotate("text", x=0.75, y=-0.75, label=percent(pn), size=5)+
  annotate("text", x=-2, y=0.75, label=percent(np), size=5)+
  annotate("text", x=-2, y=-0.75, label=percent(nn), size=5)

##----------------------------- --------- Density plot of the barcode counts ------------------------------ --------------
##data transformation, mark each row by samples using melt, for density plot
dens_tf <- function(x) {
  ID <- rownames(x)
  dat1 <- cbind(ID, x)
  colnames(dat1) <- c("ID", colnames(x))
  data <- dat1
  head(data)
  mdata <- melt(as.data.frame(data), "ID")
  data2 <-data.frame(samples=mdata$variable, logcount= log10(as.numeric(mdata$value)))
  return(data2)}
# raw count data
tfdata <- dens_tf(MNT1_DNA)
tfdata <- dens_tf(MNT1_RNA)
tfdata <- dens_tf(WM88_DNA)
tfdata <- dens_tf(WM88_RNA)
# normalized count data. All columns are scaled to have 10 million counts.
tfdata <- dens_tf(getDNA(normalize_counts(mpra_MNT)))
tfdata <- dens_tf(getRNA(normalize_counts(mpra_MNT)))
tfdata <- dens_tf(getDNA(normalize_counts(mpra_WM88)))
tfdata <- dens_tf(getRNA(normalize_counts(mpra_WM88)))
# Density plot using ggplot
p<-ggplot(data = tfdata, aes(x = logcount,fill = samples, colour = samples)) + geom_density( alpha = 0.5) + labs(x = "log10(Barcode_Counts)", y = "Density") + theme_bw() + theme(legend.position="top",legend.text=element_text(size=12),axis.text=element_text(size=14), axis.title=element_text(size=14,face=NULL, color = NULL), plot.margin=margin(t=20, r=20, b=20, l=20, "pt")) ;p
ggsave("WM88_RNA_cpm_norm.pdf", width =14, height = 14, units = "cm")
## plot using Plot_ly, click "show in a new window", and save as an svg.
ggplotly(p, width =800, height = 660) %>% config(toImageButtonOptions = list(format = "svg"))

##------------------------------- ----------- IQR plot using ggplot ------------------------------ ----------
cond <- rep(c("ref", "alt"), each = length(tfdata$samples)/2)
p<-ggplot(data = tfdata, aes(x =samples,y= logcount,fill= cond)) + geom_boxplot(outlier.colour="coral2", outlier.size = 1) + labs(y = "log10(Barcode_Counts)", x=NULL) + theme_bw() + theme(legend.position="top",legend.text=element_text(size=12),axis.text=element_text(size=14), axis.title=element_text(size=14,face=NULL, color = NULL), plot.margin=margin(t=20, r=20, b=20, l=30, "pt")) +scale_fill_brewer(palette="Accent")
ggsave("IQR_MNT1_RNA_norm.pdf", width =12, height = 10, units = "cm");p

## IQR plot using Plot_ly
plot_ly(tfdata, x = ~samples, y = ~logcount, type = "box", width = 500, height = 300, marker = list(color = "rgb(8,81,156)", outliercolor = "rgba(219, 64, 82, 0.6)", line = list(outliercolor = "rgba(219, 64, 82, 1.0)", outlierwidth = 2))) %>% plotly::layout( xaxis = list(title = "samples"), yaxis = list(title = "log10(Barcode_Counts)")) %>% config(toImageButtonOptions = list(format = "svg"))

##------------------------------------ PCA Plot----------------------------------------------
PCA_plot <- function(x) {
pca_data <- debrowser::run_pca(x)
# x is a read counts data.frame.
p_data <-  debrowser::prepPCADat(pca_data)
shape2 <- as.factor(cbind("Alt","Alt","Alt","Ref","Ref","Ref"))
pcx=1;pcy=2
xaxis <- sprintf("PC%d (%.2f%%)", pcx, round(pca_data$explained[pcx] * 100, 2))
yaxis <- sprintf("PC%d (%.2f%%)", pcy, round(pca_data$explained[pcy] * 100, 2))
plot1 <- ggplot(data = p_data, aes(x = x, y = y)) + geom_point(aes(color = color,shape =shape2), size = 4) + labs(x = xaxis, y = yaxis) + theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10, "pt")) + geom_text(aes(label = samples), vjust = 0, nudge_y = 1) 
return(plot1)
}
p <-PCA_plot(getDNA(normalize_counts(mpra_WM88)))
ggsave("PCA_WM88_DNA_norm.pdf", width =12, height = 10, units = "cm");p

#plot interactive PCA
ggplotly(p, width = 600, height = 500) %>% config(toImageButtonOptions = list(format = "svg"))


# ---------------------------------------- Volcano plot using a package-------------------------------- ---------- 
library("EnhancedVolcano")
library("scales")
# wm88 select top snp to highlight
selectLab <- toptab_WM88 %>% filter(adj.P.Val<1e-20 & abs(logFC)>0.25 ) %>% rownames()
summary(-log10(toptab_WM88$adj.P.Val))
summary(toptab_WM88$logFC)
# plotting ggplot
p<-EnhancedVolcano(toptab_WM88_M, lab=rownames(toptab_WM88), selectLab=selectLab, x="logFC", y="adj.P.Val",hlineCol="steelblue",hline= c(1e-20), borderWidth = 0.25, labSize = 4, title = NULL,titleLabSize = 16, subtitle=NULL,axisLabSize=16, pCutoff =1e-2, FCcutoff =0.2, xlim = c(-1.2, 1.2), ylim = c(0, 35), pointSize=1, colAlpha = 0.8,shape = 19,legendLabSize =15, captionLabSize = 15, cutoffLineType = 'longdash', cutoffLineWidth = 0.25, drawConnectors = T)
p 
ggsave("Volcano_WM88.pdf", width =5, height = 6, units = "in")

# MNT1
selectLab2 <- toptab_MNT %>% filter(adj.P.Val<1e-20 & abs(logFC)>0.25 ) %>% rownames()
summary(-log10(toptab_MNT$adj.P.Val))
summary(toptab_MNT$logFC)
p2<-EnhancedVolcano(toptab_MNT, lab=rownames(toptab_MNT), selectLab=selectLab2, x="logFC", y="adj.P.Val",hlineCol="steelblue",hline= c(1e-20), borderWidth = 0.25, labSize = 4, title = NULL,titleLabSize = 16, subtitle=NULL,axisLabSize=16, pCutoff =1e-2, FCcutoff =0.2, xlim = c(-2.4, 1), ylim = c(0, 120), pointSize=1, colAlpha = 0.6,shape = 19,legendLabSize =15, captionLabSize = 15, cutoffLineType = 'longdash', cutoffLineWidth = 0.25, drawConnectors = T)
p2
ggsave("Volcano_MNT1.pdf", width =5, height = 6, units = "in")

#------------------------------------ ----------- Scatter plot----------------------------- ----------------
## geom_text_repel is used to avoid overlap.
Scatt_plot <- function(Ratio,x=12,y=12,z=1){
  Ref=rowMeans(Ratio[,1:3]) 
  Alt=rowMeans(Ratio[,4:6])
  LFC=as.data.frame(log2(Ref/Alt))
  Name=rownames(FC)
  Ratio= mutate(Ratio,sig=ifelse(abs(LFC)>z ,2,1))
  p <- ggplot(Ratio, aes(x=Ref, y=Alt,color=factor(Ratio$sig))) +
    theme_bw() + geom_point(size=1,alpha=1) + 
    labs(x ="Average RNA/DNA Ratio of Ref", y="Average RNA/DNA Ratio of Alt") + 
    geom_abline(intercept = 0, slope = 1, linetype = 2, color="grey") + 
    xlim(0, x) + ylim(0, y)+ coord_fixed() + 
    theme(legend.title=NULL, legend.text=element_text(size=12), axis.text=element_text(size=12), 
          axis.title=element_text(size=12.5),plot.margin=margin(t=20, r=20, b=20, l=20, "pt")) + 
    geom_text_repel(aes(label=ifelse(abs(LFC)>z ,Name,''),hjust=-0.1, vjust=0.1),size=3, show.legend= F )
return(p)
}

Ratio <- getRNA(normalize_counts(mpra_MNT)) / getDNA(normalize_counts(mpra_MNT))
p<-Scatt_plot(Ratio,x=12,y=12,z=1) # x:xlim, y:ylim, z:fold_change.
p + scale_colour_manual(values = c("grey80", "red"),labels = c("NS", "Log2|FC|>1") )
ggsave("Scatter_MNT.pdf", width =16, height = 16, units = "cm")

Ratio <- getRNA(normalize_counts(mpra_WM88)) / getDNA(normalize_counts(mpra_WM88))
p<-Scatt_plot(Ratio,x=5.2,y=5.2,z=0.5) # x:xlim, y:ylim, z:fold_change.
p + scale_colour_manual(values = c("grey80", "red"),labels = c("NS", "Log2|FC|>0.5") )
ggsave("Scatter_WM88.pdf", width =16, height = 16, units = "cm")

##------------Bar plot of GO or KEGG --------------------
#Read in the GO enrichment dataframe or the mouse phenotype enrichment dataframe from GREAT
data <- read.table(file.choose(),header= T,sep="\t")
head(data) 
GO_rank <- cbind.data.frame(data[,1],data[,9])
colnames(GO_rank) <- c("GO_Term","qvalue")
GO_rank$qvalue <- -log10(GO_rank$qvalue)
head(GO_rank)
GO_rank <- arrange(GO_rank, desc(qvalue))
ggplot(GO_rank, aes(x= reorder(GO_Term, qvalue), y=qvalue)) + geom_bar(stat="identity",fill="#08519c", width=.5, position="dodge" ) + coord_flip() + ggtitle("GO_term enrichment by GREAT") + labs(y="-log10(Binomial q-Value)", x="GO_Term") + theme_classic() + theme(legend.title=NULL, axis.text=element_text(size=14), axis.title=element_text(size=14), plot.margin=margin(t=20, r=20, b=20, l=20, "pt")) 
ggsave("Scatter_WM88.pdf", width =8.5, height = 7.2, units = "in")

#Read in the mouse phenotype enrichment dataframe from GREAT
data <- read.table(file.choose(), sep="\t")
head(data) 
GO_rank <- cbind.data.frame(data[,1],data[,4])
head(GO_rank)
colnames(GO_rank) <- c("GO_Term","qvalue")
GO_rank$qvalue <- -log10(GO_rank$qvalue)
head(GO_rank)
GO_rank <- arrange(GO_rank, desc(qvalue))
ggplot(GO_rank, aes(x= reorder(GO_Term, qvalue), y=qvalue)) + geom_bar(stat="identity",fill="#08519c", width=.5, position="dodge" ) + coord_flip() + ggtitle("Mouse Phenotype Enrichment") + labs(y="-log10(Binomial q-Value)", x="GO_Term") + theme_classic() + theme(legend.title=NULL, axis.text=element_text(size=14), axis.title=element_text(size=14), plot.margin=margin(t=20, r=20, b=20, l=20, "pt")) 
ggsave("Scatter_WM88.pdf", width =8.5, height = 7.2, units = "in")


#-------------------------------- ---------- Volcano plot ----------------------------- ----------------
# read in the dataframe of toptable /Users/fengyq/Documents/2018_Pigmentation/2019_MPRA_data/Yuanqing_0815/0_MPRA_Rpackage_analysis/DEG/MNT1_Ref_vs_Alt.txt
toptab <- read.table(file.choose(),header =T) 
head(toptab)
Vplot<-function(toptab){
  toptab$logp <- -log10(toptab$adj.P.Val)
  snp_loc <- rownames(toptab)
  rownames(toptab) <- seq(1:length(toptab$logFC))
  toptab$color= ifelse(toptab$logp>3, 2, 1)
  p<-ggplot(data = toptab, aes(x = logFC, y = logp)) + 
    geom_point(aes(color=factor(toptab$color)),size=1,alpha=1) + 
    theme_bw() + labs(x ="Log2(Fold Change)", y="Log10(Adjust P)") + 
    xlim(-2.5, 2) + ylim(0, 9) +
    theme(legend.title=NULL, legend.text=element_text(size=12), axis.text=element_text(size=12), 
          axis.title=element_text(size=12.5), plot.margin=margin(t=20, r=20, b=20, l=20, "pt")) + 
    geom_text_repel(aes(label=ifelse(toptab$logp>3,snp_loc,''),hjust=-0.1, vjust=0.1),size=3, show.legend= F) +
    scale_colour_manual(values = c("grey80", "red"),labels = c("NS", "Log10(Adjust P)>3")) +
    geom_abline(intercept = 3, slope = 0, linetype = 2,color="grey") 
  return(p)
}
## could change the xlim and ylim scale to avoid lost extrem values.
p<-Vplot(toptab_MNT)
p
Vplot(toptab_MNT); ggsave("Volcano_MNT-1.pdf", width =16, height = 16, units = "cm")
Vplot(toptab_WM88); ggsave("Volcano_WM88-1.pdf", width =16, height = 16, units = "cm")
ggplotly(p)

