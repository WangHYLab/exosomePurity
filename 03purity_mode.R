#-----------------------upload data
{
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


zfGenes_cpm <- merge(cell_line_cpm_miExpr,normal_cpm_miExpr,by="gene_id",all = FALSE)
}

{

  median_count <- function(exosome_miExpr,col2name){
    median_tpm <- exosome_miExpr[,-1]
    median_matrix <- cbind(exosome_miExpr$gene_id, apply(median_tpm,1,median))
    median_matrix <- as.data.frame(median_matrix)
    median_matrix$V1 <- exosome_miExpr$gene_id
    colnames(median_matrix)[1]<- "gene_id"
    colnames(median_matrix)[2]<- col2name
    return(median_matrix)
  }
  
  require(quadprog)
  
  normal_cpm_median <- median_count(normal_cpm_miExpr,"normal_exo")
  cell_line_cpm_median <- median_count(cell_line_cpm_miExpr,"cell_line_exo")
  rawdata<-data.frame( merge(normal_cpm_median,cell_line_cpm_median,by="gene_id",all=T))
}


#-----------------------serum data
{
  
  tumor_mix<-read.csv("./testdata/03cancer_serum_cpm_exosomes/cancer_serum_cpm_exosomes.csv")
  tumor_mix<-tumor_mix[,-1]
  
}

mix_counts <-merge( tumor_mix,zfGenes_cpm,by="gene_id",all=FALSE)#
rownames(mix_counts)<- mix_counts[,1]
gene_id <-  mix_counts[,1]
mix_counts <-data.frame( mix_counts[,-1])
sample_name <- colnames(mix_counts)

library('preprocessCore')
mix_counts <- data.matrix(mix_counts)
mix_counts=normalize.quantiles(mix_counts)
mix_counts <- data.frame(mix_counts)
rownames(mix_counts)<- gene_id
colnames(mix_counts)<- sample_name
mix_counts <- log2(mix_counts+1)
boxplot(mix_counts)


{
  mix_counts$gene_id <-rownames(mix_counts)
  load(file = './resultdata/02signature_Var_data/CRC_var2_h20_48refgenes.Rdata')
  Lung_diff <- data.frame(common)
  colnames(Lung_diff)[1]<- "gene_id"
  mix_counts_diff <-merge(mix_counts,Lung_diff,by="gene_id",all=FALSE)
  Lung_reference <- data.frame( mix_counts_diff$gene_id)
  rownames(mix_counts_diff)<- mix_counts_diff$gene_id
  
  a=ls()
  rm(list=a[which(a!="mix_counts_diff"&a!="tumor_mix"&a!="cell_line_cpm_miExpr"&a!="normal_cpm_miExpr")])
  
  median_count_ref <- function(exosome_miExpr,col2name){
    median_tpm <- exosome_miExpr
    exosome_miExpr$gene_id <-rownames(exosome_miExpr)
    median_matrix <- as.data.frame( apply(median_tpm,1,median))
    median_matrix <- as.data.frame(median_matrix)
    median_matrix$V1 <- exosome_miExpr$gene_id
    colnames(median_matrix)[2]<- "gene_id"
    colnames(median_matrix)[1]<- col2name
    return(median_matrix)}
  
  m<-length(tumor_mix)
  n<-length(cell_line_cpm_miExpr)
  k<-length(normal_cpm_miExpr)
  
  
  
  cell_line_miExpr_median <- median_count_ref(mix_counts_diff[,c((m+1):(m+n-1))],"cell_line_exo")
  normal_miExpr_median <- median_count_ref(mix_counts_diff[,c((m+n):(m+n+k-2))],"normal_exo")
  
  simulate_miExpr<- mix_counts_diff[,c(1:m)]
  cellline_purity_miExpr<- mix_counts_diff[,c(1,(m+1):(m+n-1))]
  normal_purity_miExpr <- mix_counts_diff[,c(1,(m+n):(m+n+k-2))]
}


{
  library(quadprog)
  single_sample_purity_simulate<- function(exosome_miExpr,normal_miExpr_median,cell_line_miExpr_median){
    
    for(i in 2:length(exosome_miExpr)){
      
      temp <- exosome_miExpr[,i]
      
      temp <- data.frame(exosome_miExpr$gene_id,temp)
      colnames(temp)[2]<-"tumor_exo"
      colnames(temp)[1]<- "gene_id"
      
      temp_median<- merge(temp,normal_miExpr_median,by="gene_id",all = FALSE)
      temp_median<- merge(temp_median,cell_line_miExpr_median,by="gene_id",all = FALSE)
      
      
      temp_median <- as.data.frame(temp_median)
      normal <- as.matrix ( as.numeric( as.character(temp_median$normal_exo))) 
      cell_line <- as.matrix ( as.numeric( as.character(temp_median$cell_line_exo))) 
      X <- cbind(normal,cell_line)
      y <- temp_median$tumor_exo

      Dmat <- crossprod(X) 
      
      dvec <- as.vector(t(as.matrix(y)) %*% X) 
      bvec <- c(-1,1,rep(0, (dim(X)[2]))) 
      Amat <- rbind(-1,1,diag(dim(X)[2])) 
      
      w <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat), bvec=bvec, 
                    meq=0)$solution
      
      if(i==2){
        result1 <- data.frame(w)
        colnames(result1)[1]<- colnames(exosome_miExpr)[i]
      }else{
        result<- data.frame(w)
        colnames(result)[1]<- colnames(exosome_miExpr)[i]
        result1 <-as.data.frame(cbind(result1,result))
      }
    }
    return(result1)
  }
  
  simulate_purity <- single_sample_purity_simulate(simulate_miExpr,normal_miExpr_median,cell_line_miExpr_median)
  normal_purity <- single_sample_purity_simulate(normal_purity_miExpr,normal_miExpr_median,cell_line_miExpr_median)
  celline_purity <- single_sample_purity_simulate(cellline_purity_miExpr,normal_miExpr_median,cell_line_miExpr_median)
  
}


  write.csv(simulate_purity,file="./resultdata/03cancer_serum_purity/serum_putity.csv")







