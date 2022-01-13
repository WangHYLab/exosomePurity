#-----------------------data
tumor_miExpr<-read.csv(file= "./testdata/cancer_cellline_derived_exosomes.csv")
tumor_miExpr<-tumor_miExpr[,-1]
normal_miExpr<- read.csv(file = './testdata/healthy_cellline_derived_exosomes.csv')
normal_miExpr<-normal_miExpr[,-1]
countData <- merge(normal_miExpr,tumor_miExpr,all=FALSE,by="gene_id")
rownames(countData)<- countData$gene_id
countData<- countData[,-1]
countData[is.na(countData)] <- 0

#----------------------DEGs
library("DESeq2")
condition<- factor(c(rep("normal",41),rep("cancer",21)))
colData <- data.frame(row.names =  colnames(countData), condition)
colData
dds <- DESeqDataSetFromMatrix(countData=countData , colData=colData,design = ~condition )
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

CRC_diff_1 <- diff_gene2
save(CRC_diff_1,file = "./resultdata/CRC_0.01_200.Rdata")
