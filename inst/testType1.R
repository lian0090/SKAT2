library(SKAT2)
data(mouse)
##only fixed effects as NULL model 
y=1+rnorm(nrow(mouse.X),0,1)
H0=fitNULL(y~1)
setsize=10
nmarkers=ncol(mouse.X)
nsets=ceiling(nmarkers/setsize)
pvalues=NULL
for(i in 1:nsets){
  print(i)
  Xcols=c((setsize*(i-1)): min(i*setsize,nmarkers))
  X0=mouse.X[,Xcols]
  pvalues=rbind(pvalues,testZ(H0,Zt=X0,method=c("Score","SKAT"))$p.value)
}

##ranom effects in  NULL model 
data(mouse)
y=1+rnorm(nrow(mouse.X),0,1)
H0=fitNULL(y~1)
setsize=10
nmarkers=ncol(mouse.X)
nsets=ceiling(nmarkers/setsize)
pvalues=NULL
for(i in 1:nsets){
  print(i)
  Xcols=c((setsize*(i-1)): min(i*setsize,nmarkers))
  X0=mouse.X[,Xcols]
  pvalues=rbind(pvalues,testZ(H0,Zt=X0,method=c("Score","SKAT"))$p.value)
}



apply(pvalues,2,function(a)sum(a<=0.05)/length(a))

