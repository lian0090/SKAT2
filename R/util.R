##functions
getEigenG=function(G=NULL,Zg=NULL,precision=1e-5){
  out=list()
  if(!is.null(G) & !is.null(Zg)) stop("Only use one of G or Zg")
  if(is.null(G) & is.null(Zg)) stop("Provide G or Zg")
  if(!is.null(Zg)){
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
  class(out)=c("list","eigenG")
  return(out)
}


#Impute marker genotypes by column mean
meanImpute=function(X){
	X=apply(X,2,function(a){if(any(is.na(a))){a[which(is.na(a))]=mean(a,na.rm=T)};return(a)})
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

fit.optim = function(par, fn, logVar = T, tauRel = NULL, optimizer = "bobyqa", ...) {
	namesPar = names(par)
	if (is.null(namesPar)) {
		stop("par must have names")
	}
	if (optimizer == "optim") {
		fit <- optim(par = par, fn = fn, logVar = logVar, tauRel = tauRel, ...)
	}
	#The default fitting algorithm in lmer function
	
    if (optimizer == "bobyqa") {
		tmpfit <- bobyqa(par = par, fn = fn, logVar = logVar, tauRel = tauRel, ...)
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
	names(fit$par) = namesPar
	Var = fit$par
	Var = get_tau(Var, logVar, tauRel)
	for (i in 1:length(Var)) {
		#var_e,taud,tauw
		assign(names(Var)[i], Var[[i]])
	}
	fit$outVar = Var
	fit$loglik = -1/2 * fit$value
	if (logVar == T) {
		fit$par = exp(fit$par)
	}
	return(fit)

}

neg2Log = function(Var, tU1y, tU1X, tXX, tXy, tyy, d1, n, tU1W = NULL, tXW = NULL, tWW = NULL, tWy = NULL, kw = NULL, logVar = T, tauRel = NULL, 
	REML = T) {

	#d1 and U1 from d1=svd(Zd)$d^2, U1=svd(Zd)$u 
	Var = get_tau(Var, logVar, tauRel)
	for (i in 1:length(Var)) {
		assign(names(Var)[i], Var[[i]])
	}


	out <- getDL(var_e = var_e, taud = taud, d1 = d1, n = n, tU1y = tU1y, tU1X = tU1X, tXX = tXX, tXy = tXy, tyy = tyy, tauw = tauw, kw = kw, 
		tU1W = tU1W, tXW = tXW, tWW = tWW, tWy = tWy, getNeg2Log = T, REML = REML)$neg2logLik

	return(out)
}


#interface for getting loglikelihood for Var
getLoglik = function(Var, y, X, Z = list(), REML = T) {
	if (length(Z) < 1) {
		stop("must supply at least one random component")
	} else {
		if (!("eigenG" %in% class(Z[[1]]))) {
			cat("It is better to supply the first element of Z as an eigenG object to save computation\n")
			eigenG = getEigenG(Zg = Z[[1]])
		} else {
			eigenG = Z[[1]]
		}
		if (length(Z) == 1) {
			Zw = NULL
		} else {
			Zw = Z[-1]
		}
	}

	n = length(y)
	U1 = eigenG$U1
	d1 = eigenG$d1

	#order of elements of Var: var_e, var_d, all the var_w
	if (!is.null(Zw)) {

		if (!is.list(Zw)) {
			stop("Zw must be a list of all other random effect incidence matrix")
		}
		kw = sapply(Zw, ncol)
		nw = length(kw)
		Zw = do.call(cbind, Zw)
		if (sum(kw) != ncol(Zw)) 
			stop("sum of kw should be equal to the number of columns in Zw")
		tauw = rep(0, nw)
		tU1W = crossprod(U1, Zw)
		tXW = crossprod(X, Zw)
		tWW = crossprod(Zw, Zw)
		tWy = crossprod(Zw, y)
	} else {
		nw = 0
		kw = NULL
		tauw = NULL
		tU1W = NULL
		tXW = NULL
		tWW = NULL
		tWy = NULL
	}

	nVar = nw + 2
	if (length(Var) != nVar) 
		stop("number of elements in Var must be equal to the number of random components")
	if (nw > 0) {
		names.tauw = paste("tauw", c(1:nw), sep = "")
	} else names.tauw = NULL
	names(Var) = c("var_e", "taud", names.tauw)

	tU1X = crossprod(U1, X)
	tU1y = crossprod(U1, y)
	tXX = crossprod(X)
	tXy = crossprod(X, y)
	tyy = sum(y^2)


	out = neg2Log(Var = Var, tU1y = tU1y, tU1X = tU1X, tXX = tXX, tXy = tXy, tyy = tyy, d1 = d1, n = n, tU1W = tU1W, tXW = tXW, tWW = tWW, tWy = tWy, 
		kw = kw, tauRel = NULL, logVar = F, REML = REML)
	out = -1/2 * out
	return(out)
}


