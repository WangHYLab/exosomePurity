#-----------------------upload data

normal_miExpr1 <- read.csv("./testdata/01cancer_healthy_cellline_derived_count_exosomes/healthy_cellline_derived_exosomes.csv")
normal_miExpr<-normal_miExpr1
normal_miExpr<-normal_miExpr[,-1] 
  
cell_line_cpm_miExpr1 <- read.csv("./testdata/01cancer_healthy_cellline_derived_count_exosomes/cancer_cellline_derived_exosomes1.csv")
colnames(cell_line_cpm_miExpr1)[1]<- "gene_id"
cell_line_cpm_miExpr2 <- read.csv("./testdata/01cancer_healthy_cellline_derived_count_exosomes/cancer_cellline_derived_exosomes2.csv")
colnames(cell_line_cpm_miExpr2)[1]<- "gene_id"
cell_line_cpm_miExpr5 <- read.csv("./testdata/01cancer_healthy_cellline_derived_count_exosomes/cancer_cellline_derived_exosomes3.csv")
colnames(cell_line_cpm_miExpr5)[1]<- "gene_id"
  
cell_line_cpm_miExpr <- merge(cell_line_cpm_miExpr1,cell_line_cpm_miExpr2,by="gene_id",all=FALSE)
cell_line_cpm_miExpr <- merge(cell_line_cpm_miExpr,cell_line_cpm_miExpr5,by="gene_id",all=FALSE)

tumor_miExpr<-cell_line_cpm_miExpr

colnames(normal_miExpr)[1]<- "gene_id"
colnames(tumor_miExpr)[1]<- "gene_id"
  
countData <- merge(normal_miExpr,tumor_miExpr,all=FALSE,by="gene_id")
rownames(countData)<- countData$gene_id
countData<- countData[,-1]
countData[is.na(countData)] <- 0

#----------------------DEG analysis
library("DESeq2")
condition<- factor(c(rep("normal",31),rep("cancer",25)))
colData <- data.frame(row.names =  colnames(countData), condition)
colData

library(RUVSeq)
filter <- apply(countData, 1, function(x) length(x[x>5])>=20)
filtered <- countData[filter,]
set<-newSeqExpressionSet(as.matrix(filtered),
                         phenoData = data.frame(condition, row.names=colnames(filtered)))
library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[condition])
plotPCA(set, col=colors[condition], cex=1.2)
design <- model.matrix(~condition, data=pData(set))
y <- DGEList(counts=counts(set), group=condition)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:300]))]

set2 <- RUVg(set, empirical, k=1)
pData(set2)
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[condition])
plotPCA(set2, col=colors[condition], cex=1.2)

dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = pData(set2),
                              design = ~ W_1+condition)

dds2 <- DESeq(dds)
resultsNames(dds2)
res1 <- results(dds2,contrast=c("condition", "normal", "cancer"))
table(res1$padj<0.01)
res1 <- res1[order(res1$pvalue),]
summary(res1)
diff_gene_deseq <-subset(res1, padj < 0.01 & abs(log2FoldChange) > 1  )
diff_gene_deseq2 <- data.frame(diff_gene_deseq)
diff_gene_deseq2$gene_id <- rownames(diff_gene_deseq2)
diff_gene2 <- rownames(diff_gene_deseq2)
diff_gene2 = data.frame(diff_gene2)
colnames(diff_gene2)<- "gene_id"
save(diff_gene2,file = "./resultdata/01signature_DEGS_data/CRC_0.01_160.Rdata")

