cd ~/Dropbox/github/SKAT2/src
R CMD SHLIB matrix.c getDL.c -o getDL.so
R CMD SHLIB qfc.cpp -o qfc.so
dir="~/Dropbox/github/SKAT2"
setwd(file.path(dir,"src"))
dyn.load("getDL.so")
dyn.load("qfc.so")
library(minqa)
sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
       if(trace) cat(nm,":")
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
 }

sourceDir(file.path(dir,"R"))

load(file.path(dir,"data/mouse.RData"))
n=500
p=20
pheno=mouse.pheno[1:n,]
Z1=mouse.X[1:n,1:p]
X=cbind(model.matrix(~pheno$GENDER)[,-1],pheno$CageDensity)
X1=cbind(1,model.matrix(~pheno$GENDER)[,-1],pheno$CageDensity)
y=pheno$Obesity.BMI
G=tcrossprod(scale(mouse.X[1:n,],T,F))
eigenG=getEigenZd(Kd=G)
W=colmult(X,Z1)$Z

X=X1
Zt=Z1
tauRel=NULL
eigenZd=eigenG
U1=eigenZd$U1
d1=eigenZd$d1
windowtest=c("SKAT","Score")
tU1X=crossprod(U1,X)
tU1y=crossprod(U1,y)
tXX=crossprod(X)
tXy=crossprod(X,y)
tyy=crossprod(y,y)
tU1W=crossprod(U1,W)
tXW=crossprod(X,W)
tWW=crossprod(W)
tWy=crossprod(W,y)
tZtZt=crossprod(Zt)
tU1Zt=crossprod(U1,Zt)
tXZt=crossprod(X,Zt)
tyZt=crossprod(y,Zt)
tWZt=crossprod(W,Zt)
logVar=T
optimizer="bobyqa"
n=length(y)
    out=getDL(var_e=0.5,taud=0.5,d1=d1,n=n,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,tauw=0.5,kw=ncol(W),tU1W=tU1W,tXW=tXW,tWW=tWW,tWy=tWy,getQ=T,getS=T,tZtZt=tZtZt,tU1Zt=tU1Zt,tXZt=tXZt,tyZt=tyZt,tWZt=tWZt,getNeg2Log=T)
    