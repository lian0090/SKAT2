fit.optim=function(par,fn,logVar=T,tauRel=NULL, ...){
	namesPar=names(par)
if(is.null(namesPar)){stop("par must have names")}
  #fit<-optim(par=par,fn=fn,logVar=logVar,tauRel=tauRel, ...)
  #The default fitting algorithm in lmer function
  fit<-bobyqa(par=par,fn=fn,logVar=logVar,tauRel=tauRel, ...) 
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

getDL=function(var_e,taud,d1,n,tU1y,tU1X,tXX=NULL,tXy=NULL,tyy=NULL,tauw=NULL,kw=NULL,tU1W=NULL,tXW=NULL,tWW=NULL,tWy=NULL,getQ=F,getS=F,get.tU1ehat=T,tZtZt=NULL,tU1Zt=NULL,tXZt=NULL,tyZt=NULL,tWZt=NULL,get_tSNP=F){
  out=list()
  kd=length(d1)
  d_sharp=1/(d1*taud+var_e)
  out$d_sharp=d_sharp
   
  if(kd<n){
    d_tau=d_sharp-1/var_e
    tXU1d_tau=sweep(t(tU1X),2,d_tau,"*")
    tXVdX=tXU1d_tau%*%tU1X+tXX/var_e
    tXVdy=tXU1d_tau%*%tU1y+tXy/var_e
    out$d_tau=d_tau
  }else{
    tXU1d_sharp=sweep(t(tU1X),2,d_sharp,"*")
    tXVdX=tXU1d_sharp%*%tU1X
    tXVdy=tXU1d_sharp%*%tU1y
  }
  if(!is.null(tU1W)){
    if(kd<n){	
      tWU1d_tau=sweep(t(tU1W),2,d_tau,"*")	
      tXVdW=tXU1d_tau%*%tU1W+tXW/var_e
      tWVdW=tWU1d_tau%*%tU1W+tWW/var_e
      tWVdy=tWU1d_tau%*%tU1y+tWy/var_e
      
    }else{
      tWU1d_sharp=sweep(t(tU1W),2,d_sharp,"*")	
      tXVdW=tXU1d_sharp%*%tU1W
      tWVdW=tWU1d_sharp%*%tU1W
      tWVdy=tWU1d_sharp%*%tU1y
    }
    Gamma=rep(tauw,kw) ##kw is needed for constructing Gamma
   
    Gamma_tWVdW=sweep(tWVdW,1,Gamma,"*") #!sweep at dim 1
    Vgamma=Gamma_tWVdW
    diag(Vgamma)=diag(Vgamma+1)
    Cgamma=sweep(solve(Vgamma),2,Gamma,"*")
    tXVdW_Cgamma_tWVdX=tXVdW%*%Cgamma%*%t(tXVdW)
    tXVdW_Cgamma_tWVdy=tXVdW%*%Cgamma%*%tWVdy
    #print(tXVdW_Cgamma_tWVdX)
    #cat("\n")
    tXVinvX=tXVdX-tXVdW_Cgamma_tWVdX
    tXVinvy=tXVdy-tXVdW_Cgamma_tWVdy
    out$tWVdy=tWVdy
    out$tXVdW=tXVdW
    out$Vgamma=Vgamma
    out$Cgamma=Cgamma	

  }else{
    tXVinvX=tXVdX
    tXVinvy=tXVdy
  }
  
  invtXVinvX=solve(tXVinvX)
  hat_alpha=invtXVinvX%*%tXVinvy
  if(get.tU1ehat==T){
  tU1_ehat=tU1y-as.vector(tU1X%*%hat_alpha)
  out$tU1_ehat=tU1_ehat
  }
  out$hat_alpha=hat_alpha
  out$tXVinvX=tXVinvX
  
  if(getQ==T|getS==T){
    
    if(kd<n){
      
      tXVdZt=tXU1d_tau%*%tU1Zt+tXZt/var_e			
      tyVdZt=t(tU1y*d_tau)%*%tU1Zt+tyZt/var_e
      tZtU1d_tau=sweep(t(tU1Zt),2,d_tau,"*")
      tZtVdZt=tZtU1d_tau%*%tU1Zt+tZtZt/var_e
      if(!is.null(tauw)){tWVdZt=tWU1d_tau%*%tU1Zt+tWZt/var_e}
    }else{
      tXVdZt=tXU1d_sharp%*%tU1Zt
      tyVdZt=t(tU1y*d_sharp)%*%tU1Zt	
      tZtU1d_sharp=sweep(t(tU1Zt),2,d_sharp,"*")
      tZtVdZt=tZtU1d_sharp%*%tU1Zt
      if(!is.null(tauw))tWVdZt=tWU1d_sharp%*%tU1Zt
    }
    tehatVdZt=tyVdZt-t(hat_alpha)%*%tXVdZt		
    
    if(!is.null(tauw)){	
      tehatVdW=t(tWVdy)-t(hat_alpha)%*%tXVdW				
      LQ=tehatVdZt-tehatVdW%*%Cgamma%*%tWVdZt
      tZtVinvZt=tZtVdZt-t(tWVdZt)%*%Cgamma%*%tWVdZt
      tXVinvZt=tXVdZt-tXVdW%*%Cgamma%*%tWVdZt
    }else{
      LQ=tehatVdZt
      tZtVinvZt=tZtVdZt
      tXVinvZt=tXVdZt
    }
    
    tZtPZt=tZtVinvZt-t(tXVinvZt)%*%invtXVinvX%*%tXVinvZt
    Q=1/2*sum(LQ^2)
    if(getQ==T) {
      lambda=eigen(tZtPZt,only.values=T,symmetric=T)$values/2
      out$Q=Q
      out$lambda=lambda
    }
    
    if(getS==T){
      out$sdS=sqrt(sum(tZtPZt^2)/2)
      out$S=Q-sum(diag(tZtPZt))/2
      
    }
    
  }
  if(get_tSNP==T){
  	##has not finished yet. 
    	
    }
  return(out)
}




neg2Log=function(Var,tU1y,tU1X,tXX,tXy,tyy,d1,n,tU1W=NULL,tXW=NULL,tWW=NULL,tWy=NULL,kw=NULL,logVar=T,tauRel=NULL,REML=T){
  
  #d1 and U1 from d1=svd(Zd)$d^2, U1=svd(Zd)$u 
  Var=get_tau(Var,logVar,tauRel)
  for(i in 1:length(Var)){
  	assign(names(Var)[i],Var[[i]])
  } 
 
  
  kd=length(d1)
  terms=getDL(var_e=var_e,taud=taud,d1=d1,n=n,tauw=tauw,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,tU1W=tU1W,tXW=tXW,tWW=tWW,tWy=tWy,kw=kw,getQ=F)
  tXVinvX=terms$tXVinvX
  tU1_ehat=terms$tU1_ehat
  hat_alpha=terms$hat_alpha
  #Zd low rank (kd<n)
  if(kd<n){
    d_tau=terms$d_tau	
    logDetVd=sum(log(d1*taud+var_e))+(n-kd)*log(var_e) 
    tehat_Vd_ehat=sum(tU1_ehat^2*d_tau)+(tyy+t(hat_alpha)%*%tXX%*%hat_alpha-2*sum(hat_alpha*tXy))/var_e
  }else{
    d_sharp=terms$d_sharp
    logDetVd=sum(log(d1*taud+var_e)) 
    tehat_Vd_ehat=sum(tU1_ehat^2*d_sharp)
  }
  
  if(is.null(tU1W)){	
    neg2logLik1=logDetVd
    neg2logLik3=tehat_Vd_ehat
  }else{
    tWVdy=terms$tWVdy
    tXVdW=terms$tXVdW
    Vgamma=terms$Vgamma
    Cgamma=terms$Cgamma	
    tW_Vd_ehat=tWVdy-t(tXVdW)%*%hat_alpha	
    logDetVgamma=determinant(Vgamma,logarithm=T)$modulus
    neg2logLik1=logDetVd+logDetVgamma
    tehat_Vd_W_Cgamma_tW_Vd_ehat=t(tW_Vd_ehat)%*%Cgamma%*%tW_Vd_ehat
    neg2logLik3=tehat_Vd_ehat-tehat_Vd_W_Cgamma_tW_Vd_ehat
  } 	
  
  neg2logLik2=determinant(tXVinvX,logarithm=T)$modulus
 if(REML==T){
   out<- sum(neg2logLik1,neg2logLik2,neg2logLik3)
   #compatible with emma
   #kx=ncol(tXX)
   #out<- sum(neg2logLik1,neg2logLik2,neg2logLik3)+(n-kx)*log(2*pi)-log(det(tXX))
  }else{
   out<- sum(neg2logLik1,neg2logLik3)+n*log(2*pi)
   }
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

#A wrapper for testZ
testWindow=function(y,X,eigenG,Zt,W=NULL,removeZtFromG=F){
	
	#W should be a list of all other incidence matrix for random effects 
	#eigenG is produced by eigenZd
	if(removeZtFromG==T){
	W=c(W,list(Zt))
	nw=length(W)
	out=testZ(y,X,eigenZd=eigenG,Zt=Zt,W=W,tauRel=paste("tauw",nw,"=-taud",sep=""),windowtest=c("Score","SKAT"))
	}else{
		
	out=testZ(y,X,eigenZd=eigenG,Zt=Zt,W=W,tauRel=NULL,windowtest=c("Score","SKAT"))	
	}
	return(out)
}

testZ=function(y,X,W=NULL,tauRel=NULL,Zt,eigenZd,windowtest,nperm=0,tU1X=NULL,tU1y=NULL,tXX=NULL,tXy=NULL,tyy=NULL,logVar=T){
  if(any(is.na(y))){
  	#optim function will report not being able to evalue function at intial values when there is NA
  	stop("there should be no missing values")
  	}
  	
  out=list()
  #Null model with no random effects
  #classical SKAT test
  if(!is.null(windowtest)){
  if(is.null(eigenZd)){
  mod=lm(y~-1+X)
  resid=residuals(mod)
  s2 = summary(mod)$sigma**2
  Q=sum((t(resid)%*%(Zt))^2)/s2/2
  W.1=t(Zt) %*% Zt - (t(Zt) %*%X)%*%solve(t(X)%*%X)%*% (t(X) %*% Zt )
	lambda=eigen(W.1/2,symmetric=TRUE, only.values = TRUE)$values
	lambda1=lambda
	IDX1<-which(lambda >= 0)
	# eigenvalue bigger than mean(lambda1[IDX1])/100000 
	IDX2<-which(lambda1 > mean(lambda1[IDX1])/100000)
	lambda<-lambda1[IDX2]
	out$p.SKAT<-Get_PValue.Lambda(lambda, Q)   
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
  if(is.null(tXX)) tXX=crossprod(X)
  if(is.null(tyy)) tyy=sum(y^2)
  }
  
  if(!is.null(W)){
  	if(!is.list(W)){stop("W must be a list of all other random effect incidence matrix")}
   kw=sapply(W,ncol)
   nw=length(kw)
   W=do.call(cbind,W)
    if(sum(kw)!=ncol(W))stop("sum of kw should be equal to the number of columns in W")
    tauw=rep(0,nw)
    tU1W=crossprod(U1,W)
    tXW=crossprod(X,W)
    tWW=crossprod(W,W)
    tWy=crossprod(W,y)
    if(!is.null(windowtest)){tWZt=crossprod(W,Zt)}else{tWZt=NULL}
  }else {
    nw=0
    tauw=NULL
    tU1W=NULL
    tXW=NULL
    tWW=NULL
    tWy=NULL
    tWZt=NULL
  }
  
  if(!is.null(windowtest)){
  tU1Zt=crossprod(U1,Zt)
  tZty=crossprod(Zt,y)
  tyZt=t(tZty)
  tXZt=crossprod(X,Zt)
  tZtZt=crossprod(Zt,Zt)	
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
  fit0=fit.optim(par=par,fn=neg2Log,logVar=logVar,d1=d1,n=n,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,tU1W=tU1W,tXW=tXW,tWW=tWW,tWy=tWy,kw=kw,tauRel=tauRel)
  out$fit0=fit0
    
  ##SKAT test or LR test
  
  if("SKAT" %in% windowtest | "Score" %in% windowtest){
    var_e=fit0$outVar$var_e
    taud=fit0$outVar$taud
    if(nw>0)tauw=fit0$outVar$tauw
    getQ=("SKAT" %in% windowtest)
    getS= ("Score" %in% windowtest)
    Qdis=getDL(var_e,taud=taud,tauw=tauw,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,d1=d1,n=n,kw=kw,tU1W=tU1W,tXW=tXW,tWW=tWW,tWy=tWy,getQ=getQ,getS=getS,tZtZt=tZtZt,tU1Zt=tU1Zt,tXZt=tXZt,tyZt=tyZt,tWZt=tWZt)	
    
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
  #LR test only
  if("LR" %in% windowtest){	
    if(is.null(W)){
      tWZt=NULL
      kw_H1=ncol(Zt)
      tXW_H1=tXZt
      tU1W_H1=tU1Zt
      tWW_H1=tZtZt
      tWy_H1=tZty
    }else{
      tWZt=crossprod(W,Zt)
      kw_H1=c(kw,ncol(Zt))
      tXW_H1=cbind(tXW,tXZt)
      tU1W_H1=cbind(tU1W,tU1Zt)
      tWW_H1=rbind(cbind(tWW,tWZt),cbind(t(tWZt),tZtZt))
      tWy_H1=rbind(tWy,tZty)
    }		
    namesPar_H1=c(namesPar,paste("tauw",nw+1,sep=""))
    parH1=rep(0.5,length(namesPar_H1))
    names(parH1)=namesPar_H1
    fit1<-fit.optim(par=parH1,fn=neg2Log,logVar=logVar,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,tU1W=tU1W_H1,tXW=tXW_H1,tWW=tWW_H1,tWy=tWy_H1,d1=d1,n=n,kw=kw_H1,tauRel=tauRel)
    neg2_log0=fit0$value
    neg2_log1=fit1$value
    LR=neg2_log0-neg2_log1
    out$fit1=fit1
    out$LR=LR
    out$p.LR=pchisq(LR,1,lower.tail=F)/2
       if(nperm>0){
      LR.perm=vector()
      taut.perm=vector()
      for(j in 1:nperm){
        Ztp=Zt[sample(1:nrow(Zt)),]
        tZtZtp=crossprod(Ztp)
        tZtyp=crossprod(Ztp,y)
        tyZtp=t(tyZtp)
        tU1Ztp=crossprod(U1,Ztp)
        tXZtp=crossprod(X,Ztp)
        if(is.null(W)){
          tWZtp=NULL
          kw_H1p=ncol(Ztp)
          tXW_H1p=tXZtp
          tU1W_H1=tU1Ztp
          tWW_H1p=tZtZtp
          tWy_H1p=tZtyp
        }else{
          tWZtp=crossprod(W,Ztp)
          kw_H1p=c(kw,ncol(Ztp))
          tXW_H1p=cbind(tXW,tXZtp)
          tU1W_H1p=cbind(tU1W,tU1Ztp)
          tWW_H1p=rbind(cbind(tWW,tWZtp),cbind(t(tWZt),tZtZtp))
          tWy_H1p=rbind(tWy,tZtyp)
        }
        
        fit1<-fit.optim(par=rep(0.5,length(namesPar)+1),fn=neg2Log,namesPar=c(namesPar,paste("tauw",nw+1,sep="")),logVar=T,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,tU1W=tU1W_H1p,tXW=tXW_H1p,tWW=tWW_H1p,tWy=tWy_H1p,d1=d1,n=n,kw=kw_H1,tauRel=tauRel)
      
        #if the log likelihood at 0 is larger than at the specified value, will put taut=0  
        newVar=fit$par
        n.par=length(fit$par)
        newVar[n.par]=0
              neg2_log1.0=neg2Log(Var=newVar,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,tU1W=tU1W_H1p,tXW=tXW_H1p,tWW=tWW_H1p,tWy=tWy_H1p,d1=d1,n=n,kw=kw_H1,tauRel=tauRel,logVar=F)
        if(neg2_log1.0<fit1$value){tau2.permj=0}else{tau2.permj=fit1$par[n.par]}	
        neg2_log0=fit0$value
        neg2_log1=fit1$value
        LR=neg2_log0-neg2_log1      
        LR.perm=c(LR.perm,LR)
        tau2.perm=c(tau2.perm,tau2.permj)
      }
      out$LR.perm=LR.perm
      out$tau2.perm=tau2.perm		
    }
  }
  
  return(out)  
}

pLR.Listgarten=function(LR.perm,tau2.perm,LR,topP=0.1){
  #using the 0.1 tail is indeed better. This gives all the weight of fitting to the top 10 percent
  n=length(LR.perm)
  ntopP=round(n*topP)
  pi=mean(tau2.perm==0)
  n0=ceiling(pi*n)
  ntop=min(ntopP,n-n0)
  
  topLR=sort(LR.perm,decreasing=T)[1:ntop]
  topID=n:(n-ntop+1)
  logp.expect=log(pi+(1-pi)*(topID-0.5-n0)/(n-n0)) #quantile of the non-0 values
  
  get_ad=function(ad=c(a,d),decreasingLR.perm,logp.expect){
    a=ad[1]
    d=ad[2]
    
    logp.observe=log(pi+(1-pi)*pchisq(decreasingLR.perm/a,df=d,lower.tail=T))
    return(mean((logp.expect-logp.observe)^2))
  }
  fit=optim(par=c(1,1),fn=get_ad,decreasingLR.perm=topLR,logp.expect=logp.expect)
  a=fit$par[1]
  d=fit$par[2]
  p.LR=pchisq(LR/a,df=d,lower.tail=F)*(1-pi)
  return(list(p.LR=p.LR,a=a,d=d,p=pi))
}



pLR.Greven=function(LR.perm,LR){
	#this does not work well, because, sometimes, E(t2) might be way larger than E(t)
  Et=mean(LR.perm)
  Et2=mean(LR.perm^2)
  p=1-3*Et^2/Et2
  a=Et2/(3*Et)
  p.LR=pchisq(LR/a,df=1,lower.tail=F)*(1-p)
  out=list(p.LR=p.LR,a=a,p=p)
  return(out)
  }