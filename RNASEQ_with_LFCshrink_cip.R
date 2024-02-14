setwd('/***')
name = 'Donor_cip'
library(dplyr)
library(stringr)
library(ggplot2)
library(DESeq2)

##Block 0 - Load salmon sf quantifications into data frame counts_initial
counts_initial <- read.table("RNASEQ_cip.sf")
names(counts_initial) <- counts_initial[1,]
counts_initial <- counts_initial[-1,]
rownames(counts_initial) <- counts_initial[,1]
counts_initial <- counts_initial[,-1]
counts_initial <- counts_initial %>% mutate_if(is.character, as.numeric)


##Block 1 - Load Counts and Filter Low Coverage
counts=counts_initial[rowSums(counts_initial)>20,]
samples = data.frame(label=colnames(counts))
samples$generation <- as.factor(str_extract(samples$label, "^*[0-9a-z]+(?=_)"))
samples$replicate <- as.factor(sub(".*_", "", samples$label))
samples <- samples[order(samples$generation, samples$replicate), ]
rownames(samples) = samples$label
counts=counts[,samples$label]
remove(counts_initial)
#------------------------------------------------------
##Block 2 - Deseq2 Computations
dds <- DESeqDataSetFromMatrix(countData=round(counts),colData=samples,design = ~ replicate + generation)

dds$generation <- relevel(dds$generation, "7421")

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, fitType = 'local')
dds <- nbinomWaldTest(dds, betaPrior=FALSE)

res <- results(dds,contrast=c('generation', "7421ac", "7421"))
resLFC <- lfcShrink(dds, coef="generation_7421ac_vs_7421", type="apeglm")

res$target=rownames(res)
resLFC$target=rownames(resLFC)

summary(res)
writeLines(capture.output(summary(res)), paste0("2.0_",name,"_DESEQstats.txt"))

##Block 2 Metrics
pdf(paste0("2.1_",name,"_DispEst.pdf"),width=6,height=4,colormodel='rgb',paper = 'A4')
plotDispEsts(dds,main = "Dispersion Estimates")
dev.off()

pdf(paste0("2.2_",name,"_MA.pdf"),width=8,height=6,colormodel='rgb',paper = 'A4')
plotMA(dds,alpha=0.05,main="LFC",colNonSig="gray60",colLine="gray40",ylim=c(-4,4))
dev.off()

forpca <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")
pdf(paste0("2.3_",name,"_PCA.pdf"),width=6,height=6,colormodel='rgb',paper = 'A4')
plotPCA(forpca, intgroup = "label",ntop = 500)
dev.off()
remove(counts,dds,forpca,samples)


data <- as.data.frame(res[,c('target','baseMean','log2FoldChange','padj')])
data <- data[!is.na(data$log2FoldChange),]
data <- data[!is.na(data$padj),]
data=data[order(data$log2FoldChange),]
write.csv(data,paste0('2.4_',name,'_Gene_test.csv'))
dataFC=data[data$log2FoldChange>1 | data$log2FoldChange < -1,]
write.csv(dataFC,paste0('2.5_',name,'_GeneFC_test.csv'))

dataLFC <- as.data.frame(resLFC[,c('target','baseMean','log2FoldChange','padj')])
dataLFC <- dataLFC[!is.na(dataLFC$log2FoldChange),]
dataLFC <- dataLFC[!is.na(dataLFC$padj),]
dataLFC=dataLFC[order(dataLFC$log2FoldChange),]
write.csv(dataLFC,paste0('2.4_',name,'_Gene_LFCshrink.csv'))
dataFC_LFC=dataLFC[dataLFC$log2FoldChange>1 | dataLFC$log2FoldChange < -1,]
write.csv(dataFC_LFC,paste0('2.5_',name,'_GeneFC_LFCshrink.csv'))

