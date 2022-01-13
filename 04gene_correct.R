tumor_cpm_miExpr <- read.csv("./testdata/serum_cancer_cpm.csv")
normal_cpm_miExpr<- read.csv(file= "./testdata/serum_normal_cpm.csv")
colnames(normal_cpm_miExpr)[1] <- "gene_id"
X<-normal_cpm_miExpr
Y<-tumor_cpm_miExpr
rownames(X)<- X$gene_id
rownames(Y)<- Y$gene_id
X <- X[,-1]
Y<- Y[,-1]

rawdata = na.omit(cbind(Y,X[rownames(Y),]))
rawdata <- data.frame(rawdata)
filter_mean <- apply(rawdata, 1, mean) 
summary(filter_mean)
filter <- apply(rawdata, 1, function(x) sum(x)>100)
filtered <- rawdata[filter,]
rawdata <- filtered

{
  library('preprocessCore')
  gene_id <-  rownames(rawdata)
  sample_name <- colnames(rawdata)
  dat.transform <- data.matrix(rawdata)
  dat.transform=normalize.quantiles(dat.transform)
  befor_dat.transform <- dat.transform
  rownames(befor_dat.transform)<-gene_id
  dat.transform <- data.frame(dat.transform)
  rownames(dat.transform)<- gene_id
  colnames(dat.transform)<- sample_name
  dat.transform <- log2(dat.transform+1)
}

purity <-  read.csv("./testdata/crc_tumor_putity.csv",header = T)[,-1]
purity <- t(purity)


  mat = cbind(1, c(purity[,2], rep(0, ncol(X))))
  mat.trans = t(mat)
  invXX = solve(mat.trans%*%mat)
  H = invXX %*% mat.trans
  coefs = H %*% t(dat.transform)
  Ypred = t(mat %*% coefs)
  resi = dat.transform - Ypred
  ncase = ncol(Y)
  ncntl = ncol(X)
  s2.case = rowSums(resi[,1:ncase]^2) / (ncase - ncol(mat))
  s2.cntl = rowSums(resi[,ncase+1:ncntl]^2) / (ncntl - ncol(mat))
  
  shrinker = function(vv) {
    tmp = log(vv+1)
    exp((tmp+mean(tmp))/2)
  }
  
  s2.case = shrinker(s2.case)
  s2.cntl = shrinker(s2.cntl)
  ix = s2.case<s2.cntl
  s2.case[ix] = s2.cntl[ix] + 0.001
  H1 = H[,1:ncase]; H2 = H[,ncase+1:ncntl]
  se.new = sqrt((H1%*%t(H1))[2,2]*s2.case + (H2%*%t(H2))[2,2]*s2.cntl)
  stats = coefs[2,]/se.new#t
  df = ncol(X) + ncol(Y) - ncol(mat)
  pval = 2*pt(-abs(stats), df=df)
  out = as.data.frame(stats)
  out$pval = pval
  out$qval = p.adjust(out$pval,method = "BH")
  table(out$qval<0.05)
  
  lcf=function(Ydata,Xdata){
    #median
    Ymean <- apply(Ydata,1,median)
    Xmean <- apply(Xdata,1,median)
    log2((Ymean)/(Xmean))
    
  }
  
  befor_dat.transform <- data.frame(befor_dat.transform)
  out$log2FC <- lcf(befor_dat.transform[1:30],befor_dat.transform[31:120])
  befor_dat.transform$gene_id <- rownames(out)
  table(abs( out$log2FC)>1)
  out<-na.omit(out)
  diffgene<- out[out$qval<0.05&abs( out$log2FC)>1,]
  save(diffgene,file='./resultdata/crc_correct_113.Rdata')
  
