{
  load(file = './resultdata/CRC_0.01_200.Rdata')
  diff<-as.data.frame(CRC_diff_1)
  
  cell_line_cpm_miExpr <- read.csv("./testdata/cell_line_cpm_miExpr_var.csv")
  cell_line_cpm_miExpr<-cell_line_cpm_miExpr[,-1]
  normal_cpm_miExpr <- read.csv("./testdata/normal_cpm_miExpr_var.csv")
  normal_cpm_miExpr<-normal_cpm_miExpr[,-1]
}


fun_var<-function(ref,miEXpr){

diff<-ref
cpm_miExpr<-miEXpr
colnames(cpm_miExpr)[1]<- "gene_id"

library(preprocessCore)
temp <- cpm_miExpr
gene_id <- temp$gene_id
rownames(temp)<- temp$gene_id
temp <- temp[,-1]
temp <- data.matrix(temp)
temp <- log2(temp+1) 
boxplot(temp)

normalize_temp=normalize.quantiles(temp)
boxplot(normalize_temp)
normalize_temp <- data.frame(normalize_temp)
normalize_temp$gene_id<- gene_id

hist(normalize_temp$X5,probability = T ,main = " Histogram")
lines(density(normalize_temp$X5), col="red", lwd=2)
#

diff_temp <- merge(normalize_temp,diff,by='gene_id')
rownames(diff_temp) <- diff_temp$gene_id
diff_temp <- diff_temp[,-1]

diff_var <- matrix(ncol  = 1,nrow = length(diff_temp[,1]))
j=1
for(i in (1:length(diff_temp[,1]))){
  test <-as.numeric( as.character( diff_temp[i,]))
  diff_var[j]=var(test)
  j= j+1
}

diff_var <- data.frame(diff_var)
rownames(diff_var)<- rownames(diff_temp)
diff_var$gene_id <- rownames(diff_temp)
summary(diff_var)
diff_var <- diff_var[order(-diff_var$diff_var),]
table(diff_var$diff_var<2)

diff_var<- diff_var[diff_var$diff_var<2,]

#filter by expression
rownames(normalize_temp) <- normalize_temp$gene_id
normalize_temp<-normalize_temp[,-dim(normalize_temp)[2]]
filter<-rowSums(normalize_temp) >= 20
filter<-normalize_temp[filter,]
gene_var<-diff_var$gene_id
gene_filter<-rownames(filter)
final_gene<-intersect(gene_var,gene_filter)
return(final_gene)
}


normal<-fun_var(diff,normal_cpm_miExpr)
cancer<-fun_var(diff,cell_line_cpm_miExpr)
common<-intersect(normal,cancer)
reference_gene<-as.data.frame(common)
colnames(reference_gene)<-'gene_id'
save(reference_gene,file='./resultdata/CRC_var2_h20_51refgenes.Rdata')
