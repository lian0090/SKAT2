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

fit.optim = function(FaST, VarRel=NULL, logVar = T, optimizer = "bobyqa") {
	VarReduced=reduceInitVar(FaST$nw,VarRel)
	
	if (optimizer == "optim") {
		fit <- optim(par = VarReduced, fn = neg2LogVarReduced, FaST=FaST,VarRel=VarRel, logVar = logVar, REML=T)
	}
	#The default fitting algorithm in lmer function
	
    if (optimizer == "bobyqa") {
		tmpfit <- bobyqa(par = VarReduced, fn = neg2LogVarReduced, FaST=FaST, VarRel=VarRel, logVar = logVar, REML=T)
		fit = list()
		fit$par = tmpfit$par
		fit$value = tmpfit$fval
		fit$counts = tmpfit$feval
		fit$convergence = tmpfit$ierr
		fit$message = tmpfit$msg
	}
	if (is.na(fit$value)) {
		stop("objective function returned NA, please check input values")
	}
	namesPar=names(VarReduced)
	names(fit$par) = namesPar
	Var = fit$par
	Var = expandVar(Var, logVar, VarRel)
	fit$Var = Var
	fit$loglik = -1/2 * fit$value
	if (logVar == T) {
		fit$par = exp(fit$par)
	}
	return(fit)

}


fitNULL=function(null.formula,data=NULL,fit=T){
   FaST=getFaST(null.formula,data=data)
   out=fitNULL.FaST(FaST,fit=fit)
   return(out)
  }
   
fitNULL.FaST=function(FaST,fit=T){
   if(fit){
   y=FaST$y
   X=FaST$X
   if (is.null(FaST$U1)) {
     out=lm(y~-1+X)
	   out$FaST=FaST
	   out$Var=summary(out)$sigma^2
	   names(out$Var)="VarE"
	   #out$loglik=logLik(fit0,REML=T)
	   class(out)="lm"
	   return(out)
	} else {
	out = fit.optim(FaST, VarRel=NULL,logVar = T,  optimizer = "bobyqa" )
	out$FaST=FaST	
	class(out) = c("lmm")
	}	
	}else{
	 	out=list(FaST=FaST)
	 	if(is.null(FaST$U1)){
	 		class(out)="lm.FaST"
	 	}else{
	 		class(out)="lmm.FaST"
	 	}
	}
	return(out)
}




neg2LogVarReduced = function(VarReduced, FaST, VarRel=NULL, logVar = T, REML = T) {

	Var=expandVar(VarReduced, logVar, VarRel)

	out <- getDL(Var, FaST, getNeg2Log = T, REML = REML)$neg2logLik

	return(out)
}


#interface for getting loglikelihood for Var
getLoglik = function(Var, y, X, Z, REML = T) {	
	FaST=getFaST(y=y,X=X,Z=Z)
	
	out = getDL(Var,FaST, getNeg2Log = T, REML = REML)$neg2logLik
	out = -1/2 * out
	return(out)
}


