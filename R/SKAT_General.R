fit.optim=function(par,fn,logVar=T,tauRel=NULL, optimizer="bobyqa",...){
	namesPar=names(par)
if(is.null(namesPar)){stop("par must have names")}
  if(optimizer=="optim"){
  fit<-optim(par=par,fn=fn,logVar=logVar,tauRel=tauRel, ...)
  }
  #The default fitting algorithm in lmer function
  #
  if(optimizer=="bobyqa"){
  tmpfit<-bobyqa(par=par,fn=fn,logVar=logVar,tauRel=tauRel, ...) 
  fit=list()
  fit$par=tmpfit$par
  fit$value=tmpfit$fval
  fit$counts=tmpfit$feval
  fit$convergence=tmpfit$ierr
  fit$message=tmpfit$msg
  }
  if(is.na(fit$value)){
  stop("objective function returned NA, please check input values")
  }
  names(fit$par)=namesPar
  Var=fit$par
  Var=get_tau(Var,logVar,tauRel)
  for(i in 1:length(Var)){
  	#var_e,taud,tauw
  	assign(names(Var)[i],Var[[i]])
  } 
  fit$outVar=Var
  fit$loglik=-1/2*fit$value
  if(logVar==T){fit$par=exp(fit$par)}
  return(fit)

}

 get_tau=function(Var,logVar,tauRel){
    
  if(logVar==T){
    Var=exp(Var)
  }

  for(i in 1:length(Var)){
  	assign(names(Var)[i],Var[i])
  }
  namesPar=names(Var)
  if(any(!gsub("\\d*","",namesPar) %in% c("var_e","taud","tauw"))){
  	stop("Var names must be var_e, taud, or tauwD")
  }
 
  if(is.null(tauRel)){
  	names.tauw=grep("tauw",namesPar,value=T)
  	  	} else{
  	for(i in 1:length(tauRel)){
  		eval(parse(text=tauRel[[i]]))
  		
  	}  		
   	split.tau=strsplit(tauRel,split="=")
   	names.tauw=grep("tauw",c(unlist(split.tau),namesPar),value=T)
   	names.tauw=gsub(".*(tauw\\d+).*","\\1",names.tauw)
  	names.tauw=unique(names.tauw)
  }
  if(length(names.tauw)>0){	
  index.tauw=as.numeric(gsub("tauw(\\d+)","\\1",names.tauw))
  names.tauw=names.tauw[order(index.tauw)]
  
  tauw=vector()
  for(i in 1:length(names.tauw)){
  	tauw[i]=get(names.tauw[i])
  }
  names(tauw)=names(tauw)
  }else tauw=NULL
  Var=list(var_e=var_e,taud=taud,tauw=tauw)
  return(Var)
  }
 getDL.XYZ=function(var_e,taud,tauw=NULL,eigenZd,X,y,W=NULL,kw=NULL,Zt=NULL){
	n=length(y)
 	U1=eigenZd$U1
	d1=eigenZd$d1
	tXX=crossprod(X)
	tU1y=crossprod(U1,y)
 	tU1X=crossprod(U1,X)
 	tXy=crossprod(X,y)
 	tyy=sum(y^2) 	
 	if(!is.null(W)){
 		tU1W=crossprod(U1,W)
 		tXW=crossprod(X,W)
 		tWW=crossprod(W)
 		tWy=crossprod(W,y)
 	}else{
 		tU1W=tXW=tWW=tWy=NULL
 		
 	}
 	if(!is.null(Zt)){
 		tZtZt=crossprod(Zt)
 		tU1Zt=crossprod(U1,Zt)
 		tXZt=crossprod(X,Zt)
 		tyZt=crossprod(y,Zt)
 		if(!is.null(W)){
 			tWZt=crossprod(W,Zt)
 		}else{tWZt=NULL}	
 	}else{
 		tZtZt=tU1Zt=tXZt=tyZt=NULL
 	}
 	 out=getDL(var_e,taud,d1=d1,n=n,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,tauw=tauw,kw=kw,tU1W=tU1W,tXW=tXW,tWW=tWW,tWy=tWy,tZtZt=tZtZt,tU1Zt=tU1Zt,tXZt=tXZt,tyZt=tyZt,tWZt=tWZt,getQ=F,getS=F,getNeg2Log=T,REML=T)
 	 
 return(out)	
 }
 
neg2Log=function(Var,tU1y,tU1X,tXX,tXy,tyy,d1,n,tU1W=NULL,tXW=NULL,tWW=NULL,tWy=NULL,kw=NULL,logVar=T,tauRel=NULL,REML=T){
  
  #d1 and U1 from d1=svd(Zd)$d^2, U1=svd(Zd)$u 
  Var=get_tau(Var,logVar,tauRel)
  for(i in 1:length(Var)){
  	assign(names(Var)[i],Var[[i]])
  } 
 

  out<-getDL(var_e=var_e,taud=taud,d1=d1,n=n,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,tauw=tauw,kw=kw,tU1W=tU1W,tXW=tXW,tWW=tWW,tWy=tWy,getNeg2Log=T,REML=REML)$neg2logLik
  
  return(out)
 }


 
getDL=function(var_e,taud,d1,n,tU1y,tU1X,tXX,tXy,tyy,tauw=NULL,kw=NULL,tU1W=NULL,tXW=NULL,tWW=NULL,tWy=NULL,tZtZt=NULL,tU1Zt=NULL,tXZt=NULL,tyZt=NULL,tWZt=NULL,getQ=F,getS=F,getNeg2Log=T,REML=T)
{	if(is.null(tauw))tauw=NA
	if(is.null(kw))kw=NA
	if(is.null(tU1W))tU1W=NA
	if(is.null(tXW))tXW=NA
	if(is.null(tWW))tWW=NA
	if(is.null(tWy))tWy=NA
	if(is.null(tZtZt))tZtZt=NA
	if(is.null(tU1Zt))tU1Zt=NA
	if(is.null(tXZt))tXZt=NA
	if(is.null(tyZt))tyZt=NA
	if(is.null(tWZt))tWZt=NA
	getQ=as.integer(getQ)
	getS=as.integer(getS)
	getNeg2Log=as.integer(getNeg2Log)
	REML=as.integer(REML)
	out<-.Call("C_getDL",var_e,taud, d1,n,tU1y,tU1X,tXX,tXy,tyy,tauw,kw,tU1W,tXW,tWW,tWy,tZtZt,tU1Zt,tXZt,tyZt,tWZt,getQ,getS,getNeg2Log,REML) 
	return(out)	
	
}
	
	
getEigenG=function(G,precision=1e-5){
  out=list()
  eigenG=eigen(G,symmetric=T)
  U1=eigenG$vectors
  d1=eigenG$values	
  
  wh0=which(d1<precision)
  if(length(wh0>0)){
  	d1=d1[-wh0] 
  	U1=U1[,-wh0]
 	}
  
  out$d1=d1
  out$U1=U1
  class(out)=c("list","EigenG")
  return(out)
}




getEigenZd=function(Kd=NULL,Zd=NULL,precision=1e-5){
  out=list()
  if(!is.null(Kd) & !is.null(Zd)) stop("Only use one of Kd or Zd")
  if(!is.null(Zd)){
  	if(any(is.na(Zd))){
  		Zd=meanImpute(Zd)
  	}
  	if(nrow(Zd)<=ncol(Zd)) {
  		Kd=tcrossprod(scale(Zd,T,F))
  	    eigenKd=eigen(Kd,symmetric=T)
    	U1=eigenKd$vectors
    	d1=eigenKd$values	
  	}else{
  	svdZd=svd(Zd,nv=0)
    U1=svdZd$u
    d1=svdZd$d^2 ###Do not forget d^2!!!
  	}
  }	else{
    eigenKd=eigen(Kd,symmetric=T)
    U1=eigenKd$vectors
    d1=eigenKd$values	
  }
  
  wh0=which(d1<precision)
  if(length(wh0>0)){
  	d1=d1[-wh0] 
  	U1=U1[,-wh0]
 	}
  
  out$d1=d1
  out$U1=U1
 
  return(out)
}

      	

testZ=function(y,X,W=NULL,tauRel=NULL,Zt,eigenZd,windowtest,tU1X=NULL,tU1y=NULL,tXX=NULL,tXy=NULL,tyy=NULL,logVar=T,optimizer="bobyqa"){
	
    #check input
    X=as.matrix(X)
	Zt=as.matrix(Zt)
	if(any(is.na(y))){
  		#optim function will report not being able to evalue function at intial values when there is NA
  		if((!is.null(tU1X))|(!is.null(tU1y))|!is.null(tXX)|(!is.null(tXy))|(!is.null(tyy))){
  			stop("there should be no missing values when tU1X, tU1y, tXX, tXy or tyy is supplied")
  		}
	whNAy=which(is.na(y))
	y=y[-whNAy]
	X=X[-whNAy,,drop=F]
	Zt=Zt[-whNAy,,drop=F]
    eigenZd$U1=eigenZd$U1[-whNAy,]
    }else{
      	whNAy=NULL
    }
	out=list()
  	#Null model with no random effects
  	#classical SKAT test
 	 if(!is.null(windowtest)){
 		if(is.null(eigenZd)){
  			mod=lm(y~-1+X)
  			resid=residuals(mod)
  			s2 = summary(mod)$sigma**2
  			Q=sum(crossprod(resid,Zt)^2)/s2/2
  			tXZt=crossprod(X,Zt)
  			W.1=(crossprod(Zt) - t(tXZt) %*% solve(crossprod(X)) %*% tXZt)/2
  			lambda=eigen(W.1,symmetric=TRUE, only.values = TRUE)$values
			lambda1=lambda
			IDX1<-which(lambda >= 0)
			# eigenvalue bigger than mean(lambda1[IDX1])/100000 
			IDX2<-which(lambda1 > mean(lambda1[IDX1])/100000)
			lambda<-lambda1[IDX2]
			varS=(sum(W.1^2)/(s2^2))/2
			S=(Q-(sum(diag(W.1))))/s2
			out$p.SKAT<-Get_PValue.Lambda(lambda, Q)
			out$Q=Q
			out$S=S   
			out$p.Score=pnorm(S,mean=0,sd=sqrt(varS),lower.tail=F)
  			return(out)
  		}
  }
  
  #logVar, paramterize variance components with log when using REML to restrict variance component to be larger than 0.
  	
  d1=eigenZd$d1
  U1=eigenZd$U1
  
  if(is.null(tU1X)) tU1X=crossprod(U1,X)
  if(is.null(tU1y)) tU1y=crossprod(U1,y)
  
  n=length(y)
  nd=length(d1)
  
  
  if(nd<n){
  if(is.null(tXy)) tXy=crossprod(X,y)
  if(is.null(tXX)) tXX=crossprod(X,X)
  if(is.null(tyy)) tyy=sum(y^2)
  }
  
  if(!is.null(W)){
  	
  	if(!is.list(W)){stop("W must be a list of all other random effect incidence matrix")}
  	kw=sapply(W,ncol)
    nw=length(kw)
    W=do.call(cbind,W)
    if(sum(kw)!=ncol(W))stop("sum of kw should be equal to the number of columns in W")
    #remove NA values
    if(length(whNAy)>0){W=W[-whNAy,,drop=F]} 
    tauw=rep(0,nw)
    tU1W=crossprod(U1,W)
    tXW=crossprod(X,W)
    tWW=crossprod(W,W)
    tWy=crossprod(W,y)
    if(!is.null(windowtest)){tWZt=crossprod(W,Zt)}else{tWZt=NULL}
  }else {
    nw=0
    kw=NULL
    tauw=NULL
    tU1W=NULL
    tXW=NULL
    tWW=NULL
    tWy=NULL
    tWZt=NULL
  }
  
  if(!is.null(windowtest)){
  	tU1Zt=Crossprod(U1,Zt)
  	tyZt=crossprod(y,Zt)
  	tXZt=crossprod(X,Zt)
  	tZtZt=crossprod(Zt)	
  }
  
  ##test with low rank Zh	
  namesPar=c("var_e","taud")
  if(nw>0){
  	namesPar=c("var_e","taud",paste("tauw",c(1:nw),sep=""))
  	if(!is.null(tauRel)){
  		splittau=strsplit(tauRel,split="=")
  		unlist.splittau=unlist(splittau)
  		#only the independent Var need to be positive if using logVar
  		#exp transformation will only be applied to par, not the dependent var
  		#if(any(grepl("-[[:punct:][:digit:]]*tau",unlist.splittau))){if(logVar==T)stop("logVar must be set to False when the tau values are of different signs")}
  		n.Rel=length(splittau)
  		if(length(unlist.splittau)>2*length(splittau)) stop("Every equation must contain only one relationship")
  		#strip off the digits, decimal point (.) and +,-,*,/
  		if(any(!(gsub("([[:digit:][:punct:]e]+)","",unlist.splittau)%in% c("taud","tauw","")))){
  			stop("in tauRel, must only specify taud or tauwD, where D is any integer less or equal to the number of W matrix")
  			}
  		left.tau=gsub(".*(tau[wd]{1}\\d*).*","\\1",sapply(splittau,function(a)a[1]))
        right.tau=gsub(".*(tau[wd]{1}\\d*).*","\\1",sapply(splittau,function(a)a[2]))
        if(any(duplicated(left.tau))){stop("duplicated terms on the lest side of equations")}
        if(length(intersect(right.tau,left.tau))>0){stop("one variable can only appear on one side of equation")}
  		
  		for(i in 1:length(splittau)){
  			splittau.i=gsub(".*(tau[dw]{1}\\d*).*","\\1",splittau[[i]])
  			
  			n.split=length(splittau.i)
  			##only keep the last tau from relationship equation
  			namesPar=setdiff(namesPar,splittau.i[-n.split])
  		}
  		}
  	}
  par=rep(0.5,length(namesPar))
  names(par)=namesPar	
 
  fit0=fit.optim(par=par,fn=neg2Log,logVar=logVar,d1=d1,n=n,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,tU1W=tU1W,tXW=tXW,tWW=tWW,tWy=tWy,kw=kw,tauRel=tauRel,optimizer=optimizer)
  out$fit0=fit0
    
  ##SKAT test or LR test
  
  if("SKAT" %in% windowtest | "Score" %in% windowtest){
    var_e=fit0$outVar$var_e
    taud=fit0$outVar$taud
    if(nw>0)tauw=fit0$outVar$tauw
    getQ=("SKAT" %in% windowtest)
    getS= ("Score" %in% windowtest)
    Qdis=getDL(var_e,taud=taud,d1=d1,n=n,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,tauw=tauw,kw=kw,tU1W=tU1W,tXW=tXW,tWW=tWW,tWy=tWy,tZtZt=tZtZt,tU1Zt=tU1Zt,tXZt=tXZt,tyZt=tyZt,tWZt=tWZt,getQ=getQ,getS=getS,getNeg2Log=F,REML=T)	
    
    if("SKAT" %in% windowtest){
      Q=Qdis$Q
      lambda=Qdis$lambda      
      lambda1=lambda
      IDX1<-which(lambda >= 0)
      #eigenvalue bigger than mean(lambda1[IDX1])/100000 
      IDX2<-which(lambda1 > mean(lambda1[IDX1])/100000)
      lambda<-lambda1[IDX2]
      p.value<-Get_PValue.Lambda(lambda,Q) 
      out$p.SKAT=p.value
      out$Q=Q
    }
    #Score test
    if("Score" %in% windowtest){
      S=Qdis$S
      sdS=Qdis$sdS	
      out$Score=S
      out$sdScore=sdS
      out$p.Score=pnorm(S/sdS,lower.tail=F)
    }
    
  }
 return(out)  
}

