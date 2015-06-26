#t-test on individual markers
#It is not exactly the same as that from emma, due to the small difference in estimating variance components. 
#If I use the same variance component from emma to put into getDL, I will get exactly the same pvalue 
#Population structure previously determined. 
P3D.NULL=function(y,X0,eigenG){
	#X0: fixed effects not being tested
	X=as.matrix(X0)
	n=length(y)
	U1=eigenG$U1
	d1=eigenG$d1
	tXX=crossprod(X)
	tU1y=crossprod(U1,y)
 	tU1X=crossprod(U1,X)
 	tXy=crossprod(X,y)
 	tyy=sum(y^2)
    Var=c(0.5,0.5)
    names(Var)=c("var_e","taud")
	fit0<-fit.optim(par=Var,fn=neg2Log,logVar=T,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,d1=d1,n=n)
	Var=fit0$par
	names(Var)=c("var_e","var_g")
	return(Var)
}
singleSNP.P3D=function(y,X0,Xt,Var,eigenG,method="LR"){
 #X0: incidence matrix for fixed effects not being tested
 #Xt: incidence matrix for the marker or for the GxE  to be tested	
 #Var: population variance components: var_e and var_g 
 #methods: "t" for t-test, "LR" for likelihood ratio test, default is likelihood ratio test
    
    X0=as.matrix(X0)
    Xt=as.matrix(Xt)
  if(any(is.na(y))){
  	#optim function will report not being able to evalue function at intial values when there is NA
  	#stop("there should be no missing values")
	whNA=which(is.na(y))
	y=y[-whNA]
	X0=X0[-whNA,,drop=F]
	Xt=Xt[-whNA,,drop=F]
    eigenG$U1=eigenG$U1[-whNA,]
  	}
    
   
    X=cbind(X0,Xt)
    names(Var)=c("var_e","taud")
    out=list()
 	n.beta=ncol(X)		
 	test=c((ncol(X0)+1):n.beta)	
 	
 	if(length(test)>1){
 		if("t" %in% method){
 			stop("use LR test when there is more than one fix effect to be tested")
 		}
 	}
 	#use LR test if length.test>1
 	out$p.value=rep(0,length(method))
 	names(out$p.value)=method
 	for(i in 1:length(method)){
 		
 	if("LR" == method[i]){
 		trypvalue=try({
 		ln0=getLoglik(Var=Var,y,X=X0,eigenZd=eigenG,logVar=F,REML=F)
 	 	ln1=getLoglik(Var=Var,y,X=X,eigenZd=eigenG,logVar=F,REML=F)
 	 	Q=-2*(ln0-ln1)
 	 	p.value=pchisq(Q,df=length(test),lower.tail=F)
 	 	out$ML1=ln1
 	 	out$ML0=ln0
 	 	out$LR=Q
 	 	})
 	 	}else if ("t" == method[i]){
 	n=length(y)
	U1=eigenG$U1
	d1=eigenG$d1
	tXX=crossprod(X)
	tU1y=crossprod(U1,y)
 	tU1X=crossprod(U1,X)
 	tXy=crossprod(X,y)
 	tyy=sum(y^2)
 	var_e=Var[1]
 	taud=Var[2]
 	trypvalue=try({
 	outDL=getDL(var_e=var_e,taud=taud,d1=d1,n=n,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,get.tU1ehat=F)
 	beta=outDL$hat_alpha
 	vbeta=solve(outDL$tXVinvX)
 	tscore=beta[test]/sqrt(vbeta[test,test])
 	##note: the df for t-distribution is not corrected by Satterthwaite's method. Likelihood ratio test should be better.
 	p.value=2*pt(tscore,df=n-n.beta,lower.tail=F)
 	})
 	}
   
   if(inherits(trypvalue, "try-error")){
  		p.value=NA
  		}	
 	 out$p.value[i]=p.value
 	 }	
    return(out)
 }
 ##perform association mapping for provided markers while correcting for multiple test.
 GWAS.P3D=function(y,X0,Xt,Var,eigenG,multipleCorrection=T){
 	   	 Xt=as.matrix(Xt)
 	   	 X0=as.matrix(X0)
 	   	 if(multipleCorrection==T){
 	   	 	Me=Meff(Xt)
 	   	 	cat("effetive number of test is ",Me,"\n")
 	   	 	}else{
 	   	 		Me=1
 	   	 	}
 	   	 p.value=rep(NA,ncol(Xt))
 	   	 for(i in 1:ncol(Xt)){
 	   	 p.value[i]=singleSNP.P3D(y,X0=X0,Xt=Xt[,i],Var=Var,eigenG=eigenG,method="LR")$p.value*Me
 	   	 }
 	   	 
 	   	 return(list(Me=Me,p.value=p.value))
 	   	}

 	
 
 ###population parameter re-estimated for each marker
 singleSNP=function(y,X,eigenG){
 	n=length(y)
	U1=eigenG$U1
	d1=eigenG$d1
	tXX=crossprod(X)
	tU1y=crossprod(U1,y)
 	tU1X=crossprod(U1,X)
 	tXy=crossprod(X,y)
 	tyy=sum(y^2)
    Var=c(0.5,0.5)
    names(Var)=c("var_e","taud")
	fit0<-fit.optim(par=Var,fn=neg2Log,logVar=T,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,d1=d1,n=n)
	Var=fit0$par
    outDL=getDL(var_e=Var[1],taud=Var[2],d1=d1,n=n,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,get.tU1ehat=F)
 	beta=outDL$hat_alpha
 	vbeta=solve(outDL$tXVinvX)
 	n.beta=length(beta)
 	tscore=beta[n.beta]/sqrt(vbeta[n.beta,n.beta])
 	##note: the df for t-distribution is not corrected by Satterthwaite's method. Likelihood ratio test should be better.
 	p.value=2*pt(tscore,df=n-n.beta,lower.tail=F)
    return(p.value)
 }