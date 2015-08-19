
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


