library(SKAT2)
data(mouse)
n=1825
p=30
pheno=mouse.pheno[1:n,]
Z1=mouse.X[1:n,1:p]
X=cbind(model.matrix(~pheno$GENDER)[,-1],pheno$CageDensity)
X1=cbind(1,model.matrix(~pheno$GENDER)[,-1],pheno$CageDensity)
y=pheno$Obesity.BMI
G=tcrossprod(scale(mouse.X[1:n,],T,F))
eigenG=getEigenZd(Kd=G)
W=colmult(X,Z1)$Z


Rprof("SKATmem.out",memory.profiling = TRUE)
out=testWindow(y,X=X1,Zt=Z1,W=list(W),eigenG=eigenG)
Rprof(NULL)
summaryRprof("SKATmem.out",memory="both")


Rprof("GibbsI.out",memory.profiling = TRUE)
FW(y,VAR,ENV)
Rprof(NULL)
summaryRprof("GibbsI.out",memory="both")

Rprof("GibbsOLS.out",memory.profiling = TRUE)
FW(y,VAR,ENV,method="OLS")
Rprof(NULL)
summaryRprof("GibbsI.out",memory="both")