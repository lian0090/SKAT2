##functions
getEigenG=function(G=NULL,Zg=NULL,precision=1e-5){
  out=list()
  if(!is.null(G) & !is.null(Zg)) stop("Only use one of G or Zg")
  if(is.null(G) & is.null(Zg)) stop("Provide G or Zg")
  if(!is.null(Zg)){
  	Zg=as.matrix(Zg)
  	if(any(is.na(Zg))){
  		Zg=meanImpute(Zg)
  	}
  	if(nrow(Zg)<=ncol(Zg)) {
  		G=tcrossprod(Zg)
  	    eigenG=eigen(G,symmetric=T)
    	U1=eigenG$vectors
    	d1=eigenG$values	
  	}else{
  	svdZg=svd(Zg,nv=0)
    U1=svdZg$u
    d1=svdZg$d^2 ###Do not forget d^2!!!
  	}
  }else{
    eigenG=eigen(G,symmetric=T)
    U1=eigenG$vectors
    d1=eigenG$values	
  }
  
  wh0=which(d1<precision)
  if(length(wh0>0)){
  	d1=d1[-wh0] 
  	U1=U1[,-wh0]
 	}
  
  out$d1=d1
  out$U1=U1
  class(out)=c("eigenG")
  return(out)
}


#Impute marker genotypes by column mean
meanImpute=function (X) 
{
  X=as.matrix(X)
  X = apply(X, 2, function(a) {
    if (any(is.na(a))) {
      a[which(is.na(a))] = mean(a, na.rm = T)
    }
    return(a)
  })
  return(X)
}

	

Meff=function(X){
	X=as.matrix(X)
	#X is the incidence matrix for the variables to be tested
	if(ncol(X)==1){return(1)}
	corX=cor(X)
	lambda=eigen(corX,only.values=T)$values
	return(Meff.lambda(lambda))
}

Meff.lambda=function(lambda){
	##round lambda first, so that if the floor will work correctly in the case of numerical precision. 
	#for example: a=rnorm(100);b=a+1;d=a+3;eigen(cor(cbind(a,b,d))...
	lambda=round(lambda,digits=4)
	x1=as.numeric(lambda>=1)+(lambda-floor(lambda))
	return(sum(x1))
	}



#calculate matrix trace
tr=function(X){
	out=sum(diag(X))
	return(out)
}	


	
	



fitNULL=function(null.formula,data=NULL){ 
   RandomTerm=findCallPattern(null.formula,c(".R",".G",".eigenG"))
   if(!is.null(RandomTerm)){
     FaST=getFaST(null.formula,data=data)
     out=fitNULL.FaST(FaST)
     return(out)
   }else{
     out=lm(null.formula,data=data)
     out$FaST=getFaST(null.formula,data=data)
     out$Var=summary(out)$sigma^2
     names(out$Var)="var_e"
     out$loglik=logLik(out,REML=T)
     class(out)="lm"
     return(out) 
   }
  }

#FaST changed, but Var not changing, for P3D methods
updatefitNULLFaST=function(fit0,FaST){
  if(class(fit0)=="lm") stop("cannot update fit0 inherited from lm")
  fit0$FaST=FaST
  #fit0$loglikML=-1/2*getDL(fit0$Var, FaST, getNeg2Log = T, REML = F)$neg2logLik
  if(FaST$method=="FaST"|FaST$method=="brute"){
  fit0$loglikREML=-1/2*getDL(fit0$Var, FaST, getNeg2Log = T, REML = T)$neg2logLik
  }
}

neg2Log=function(Var,FaST,logVar,REML){
  #logVar:whether the input Var is in log scale.
  if(logVar)Var=exp(Var)
  getDL(Var, FaST, getNeg2Log = T, REML = REML)$neg2logLik
}

fitNULL.FaST=function(FaST){
  method=FaST$method
  if(method=="FaST"|method=="brute"){
    n.Var=FaST$n.randomTerm+1
    Var=rep(0.5,n.Var)
    tmpfit <- bobyqa(par = Var, fn = neg2Log, FaST=FaST,logVar = T, REML=T)
    Var = exp(tmpfit$par)
    names(Var)=c("var_e",FaST$namesRandomTerm)
    tmpfit$value = tmpfit$fval 
    out = list()
    out$counts = tmpfit$feval
    out$convergence = tmpfit$ierr
    out$message = tmpfit$msg    
    if (is.na(tmpfit$value)) {
      stop("objective function returned NA, please check input values")
    }
    out$Var = Var
    out$loglikREML = -1/2 * tmpfit$value
    out$FaST=FaST
    class(out)="lmm" 
    return(out)
  }else stop("only FaST and brute are supported in fitNULL.FaST")  
}








