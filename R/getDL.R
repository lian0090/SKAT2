
#Z=sweep(op(X),v,MARGINopX,"*")
sweep_prod=function(X,v,transX,MARGINopx){
 return(.Call("Rsweep_prod",X,v,as.integer(transX),MARGINopx) )
}


##does not depend on the real names of Var. 
getDL = function(Var, FaST, getQ = F, getS = F, getNeg2Log = T, REML = T, getAlphaHat=F) {
  if(FaST$method=="FaST"){
  var_e=Var[1]
  var_d=Var[FaST$whichFaST+1]
  
  if (FaST$nw == 0) {
		VarW = NULL
	} else {
		VarW = Var[-c(1,FaST$whichFaST+1)]
	}
	##.Call will not be able to take NULL values. 
	FaSTnames=c("d1","n","tU1y","tU1X","tXX","tXy","tyy","kw","nw","tU1W","tXW","tWW","tWy","tZtZt","tU1Zt","tXZt","tyZt","tWZt")
	for(i in 1:length(FaSTnames)){
		namei=FaSTnames[i]
		if(is.null(FaST[[namei]])){
			FaST[[namei]]=NA
		}
	}
	
	out <- .Call("C_getDL", var_e, var_d, FaST$d1, FaST$n, FaST$tU1y, FaST$tU1X, FaST$tXX, FaST$tXy, FaST$tyy, VarW, FaST$kw, FaST$nw, FaST$tU1W, FaST$tXW, FaST$tWW, FaST$tWy, FaST$tZtZt, FaST$tU1Zt, FaST$tXZt, FaST$tyZt, 
		FaST$tWZt, as.integer(getQ), as.integer(getS), as.integer(getNeg2Log), as.integer(REML),as.integer(getAlphaHat))
	return(out)
  }
  if(FaST$method=="brute"){
  	out=vector(mode="list",length=7)
    names(out)=c("neg2logLik","Q","lambda","S","sdS","hat_alpha","invtXVinvX")
    y=FaST$y
    X=FaST$X
    listG=FaST$listG
    Ztest=FaST$Ztest
	if (length(Var) != (length(listG) + 1)) 
		stop("number of variance must be equal to 1 plus the number of G in listG\n")
	var_e = Var[1]
	var_z = Var[-1]
	n = nrow(as.matrix(X))
	nG = length(listG)
	V = diag(var_e, n)
	#if we use LL'=V, the calculations might be more efficient
	for (i in 1:nG) {
		V = V + listG[[i]] * var_z[i]
	}
	Vinv = .Call("C_matinv", V)
	tXVinv=crossprod(X,Vinv)
	tXVinvX = tXVinv %*% X
	invtXVinvX=.Call("C_matinv",tXVinvX)
	P = Vinv - t(tXVinv) %*% invtXinvX %*% tXVinv
    if(getAlphaHat)	{
    	tXVinvy=tXVinv %*% y
    	out$hat_alpha=invtXVinvX%*%tXVinvy
    	out$invtXVinvX=invtXVinvX
    }
	
	if(getNeg2Log){
	Py=P%*%y
	neg2logLik1=determinant(V, logarithm = T)$modulus
	neg2logLik3=as.numeric(t(y) %*% Py)

	if(REML){
	neg2logLik2=determinant(tXVinvX, logarithm = T)$modulus
	neg2logLik= neg2logLik1 + neg2logLik2  + neg2logLik3 		}else{
	neg2logLik= neg2logLik1+neg2logLik3+n*log(2*pi)		
		}
	out$neg2loglik=neg2logLik	
		}
	if(getQ|getS){
		
		tZtP=crossprod(Ztest,P)
		RQ=tZtP%*%y
		Q = sum(RQ^2)/2
		tZtPZt = (tZtP%*% Ztest)
		if(getQ){
		out$Q=Q
		out$lambda = (eigen(tZtPZt, symmetric = TRUE, only.values = TRUE)$values/2)
		}
		if(getS){
		out$S = Q - sum(diag(tZtPZt))/2
        out$sdS=sqrt(sum(tZtPZt^2)/2)
		}
	}	
   return(out)
  	
  }
  
}






