#-----------------------upload data
normal_cpm_miExpr<- read.csv("./testdata/04correct_DEG_gene/healthy_cpm_exosomes.csv")
colnames(normal_cpm_miExpr)[1]<- "gene_id"
cell_line_cpm_miExpr<- read.csv("./testdata/04correct_DEG_gene/cancer_cpm_exosomes.csv")[,-1]
colnames(cell_line_cpm_miExpr)[1]<- "gene_id"
tumor_cpm_miExpr<-cell_line_cpm_miExpr 

X<-normal_cpm_miExpr
Y<-tumor_cpm_miExpr
rownames(X)<- X$gene_id
rownames(Y)<- Y$gene_id
X <- X[,-1]
Y<- Y[,-1]

ny<-length(Y)
nx<-length(X)

rawdata = na.omit(cbind(Y,X[rownames(Y),]))
rawdata <- data.frame(rawdata)
filter_mean <- apply(rawdata, 1, mean) 
summary(filter_mean)
filter <- apply(rawdata, 1, function(x) sum(x)>100)
filtered <- rawdata[filter,]
rawdata <- filtered

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

#--------------------------------upload serum purity
purity <-  read.csv("./resultdata/03cancer_serum_purity/serum_putity.csv",header = T)[,-1]
purity <- t(purity)

#-------------------------------correct
mat = cbind(1, c(purity[,2], rep(0, ncol(X))))#W
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

stats = coefs[2,]/se.new
df = ncol(X) + ncol(Y) - ncol(mat)
pval = 2*pt(-abs(stats), df=df)
out = as.data.frame(stats)
out$pval = pval
out$qval = p.adjust(out$pval,method = "BH")
table(out$qval<0.05)


lcf=function(Ydata,Xdata){
  Ymean <- apply(Ydata,1,mean)
  Xmean <- apply(Xdata,1,mean)
  log2((Ymean)/(Xmean))
  
}


befor_dat.transform <- data.frame(befor_dat.transform)
out$log2FC <- lcf(befor_dat.transform[,c(1:ny)],befor_dat.transform[,c((ny+1):(ny+nx))])
befor_dat.transform$gene_id <- rownames(out)
table(abs( out$log2FC)>1)
out<-na.omit(out)
diffgene<- out[out$qval<0.05&abs( out$log2FC)>1,]
diffgene<- out[out$pval<0.05&abs( out$log2FC)>1,]
diff<-diffgene

save(diff,file = './resultdata/04correct_DEG_gene/01after_crc_p0.05_71.Rdata')

load('./testdata/04correct_DEG_gene/01before_crc_0.05_188.Rdata')

diff_common<-intersect(rownames(diffgene),diff_gene2$gene_id)
n<-match(diff_common,rownames(diffgene))
diff_new<-rownames(diffgene)[-n]


save(diff_new,file = './resultdata/04correct_DEG_gene/01new_crc_p0.05_27.Rdata')
