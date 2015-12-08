
setwd("~/Dropbox/github/SKAT2/inst")
library(SKAT2)
data(mouse)
n=1825
p=30
pheno=mouse.pheno[1:n,]
Z1=mouse.X[1:n,1:p]
X=cbind(model.matrix(~pheno$GENDER)[,-1],pheno$CageDensity)
X1=cbind(1,model.matrix(~pheno$GENDER)[,-1],pheno$CageDensity)
y=pheno$Obesity.BMI
#G=tcrossprod(scale(mouse.X[1:n,],T,F))
#eigenG=getEigenZd(Kd=G)

#save(G,file="G.rda")
#save(eigenG,file="eigenG.rda")

load("G.rda")
load("eigenG.rda")

Zt=colmult(X,Z1)$Z


Rprof("SKATmem.out",memory.profiling = TRUE)
out=testWindow(y,X=X1,Zt=Zt,W=list(Z1),eigenG=eigenG)
Rprof(NULL)
summaryRprof("SKATmem.out",memory="both")


##test on hpcc
library(BEDMatrix)
setwd("~/BC1958")
pheno=read.table("dataClean.fam",sep="",head=F)
##simulate phenotype and a covariate
sex=factor(pheno[,5])
set.seed(12345)
X2=cbind(model.matrix(~sex)[,-1],rnorm(length(sex)))
X=cbind(1,X2)
y=rnorm(length(sex))
geno=BEDMatrix('dataClean.bed', n=2994,p=657756)
eigenGA=readRDS("eigenGA.rds")
Z1=geno[,1:30]
Zt=colmult(X2,Z1)$Z

#fit.optim 0.58 s for SKAT2
library(SKAT2)
Rprof("SKATmem.out",memory.profiling = TRUE)
out=testWindow(y,X=X,Zt=Zt,W=list(Z1),eigenG=eigenGA)
Rprof(NULL)
summaryRprof("SKATmem.out",memory="both")

##fit.optim 6.34 for RSKAT2
library(RSKAT2)
Rprof("SKATmemR.out",memory.profiling = TRUE)
out=testWindow(y,X=X,Zt=Zt,W=list(Z1),eigenG=eigenGA)
Rprof(NULL)
summaryRprof("SKATmemR.out",memory="both")

