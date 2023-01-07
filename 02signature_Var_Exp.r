#-----------------------upload data
{
  
  load(file = './resultdata/01signature_DEGS_data/CRC_0.01_160.Rdata')
  diff<-as.data.frame(diff_gene2)
  normal_miExpr1 <- read.csv("./testdata/02cancer_healthy_cellline_derived_cpm_exosomes/healthy_cellline_derived_cpm_exosomes.csv")
  normal_cpm_miExpr<-normal_miExpr1
  normal_cpm_miExpr<-normal_cpm_miExpr[,-1]
  
  cell_line_cpm_miExpr1 <- read.csv("./testdata/02cancer_healthy_cellline_derived_cpm_exosomes/cancer_cellline_derived_cpm_exosomes1.csv")
  colnames(cell_line_cpm_miExpr1)[1]<- "gene_id"
  cell_line_cpm_miExpr2 <- read.csv("./testdata/02cancer_healthy_cellline_derived_cpm_exosomes/cancer_cellline_derived_cpm_exosomes2.csv")
  colnames(cell_line_cpm_miExpr2)[1]<- "gene_id"
  cell_line_cpm_miExpr5 <- read.csv("./testdata/02cancer_healthy_cellline_derived_cpm_exosomes/cancer_cellline_derived_cpm_exosomes3.csv")
  colnames(cell_line_cpm_miExpr5)[1]<- "gene_id"
  
  cell_line_cpm_miExpr <- merge(cell_line_cpm_miExpr1,cell_line_cpm_miExpr2,by="gene_id",all=FALSE)
  cell_line_cpm_miExpr <- merge(cell_line_cpm_miExpr,cell_line_cpm_miExpr5,by="gene_id",all=FALSE)
}

#----------------------var
fun_var<-function(deg,exp){
  
  diff<-deg
  cpm_miExpr<-exp
  colnames(cpm_miExpr)[1]<- "gene_id"
  
  library(preprocessCore)
  temp <- cpm_miExpr
  gene_id <- temp$gene_id
  rownames(temp)<- temp$gene_id
  temp <- temp[,-1]
  temp <- data.matrix(temp)
  temp <- log2(temp+1) 
  boxplot(temp)
  #---yes
  normalize_temp=normalize.quantiles(temp)
  boxplot(normalize_temp)
  normalize_temp <- data.frame(normalize_temp)
  normalize_temp$gene_id<- gene_id
  
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

save(common,file='./resultdata/02signature_Var_data/CRC_var2_h20_48refgenes.Rdata')



