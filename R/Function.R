##functions
#Impute marker genotypes by column mean
meanImpute=function(X){
	X=apply(X,2,function(a){if(any(is.na(a))){a[which(is.na(a))]=mean(a,na.rm=T)};return(a)})
	return(X)
	}
	


#calculate matrix trace
tr=function(X){
	out=sum(diag(X))
	return(out)
}	


#get loglikelihood for Var
getLoglik=function(Var,y,X,W=NULL,kw=NULL,eigenZd,logVar,tauRel=NULL,REML=T){
if(is.null(names(Var))){stop("Var must have names")}	
n=length(y)
U1=eigenZd$U1
d1=eigenZd$d1
tU1X=crossprod(U1,X)
tU1y=crossprod(U1,y)
tXX=crossprod(X)
tXy=crossprod(X,y)
tyy=sum(y^2)
  if(!is.null(W)){
    if(is.null(kw))stop("kw must be speficied for W")
    nw=length(kw)
    if(sum(kw)!=ncol(W))stop("sum of kw should be equal to the number of columns in W")
    tauw=rep(0,nw)
    tU1W=crossprod(U1,W)
    tXW=crossprod(X,W)
    tWW=crossprod(W,W)
    tWy=crossprod(W,y)
  }else {
    nw=0
    tauw=NULL
    tU1W=NULL
    tXW=NULL
    tWW=NULL
    tWy=NULL
  }

out=neg2Log(Var=Var,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,d1=d1,n=n,tU1W=tU1W,tXW=tXW,tWW=tWW,tWy=tWy,kw=kw,tauRel=tauRel,logVar=logVar,REML=REML)
out=-1/2*out
return(out)
}

