#t-test on individual markers
#It is not exactly the same as that from emma, due to the small difference in estimating variance components. 
#If I use the same variance component from emma to put into getDL, I will get exactly the same pvalue 
#Population structure previously determined. 
P3D.NULL=function(y,X,eigenG){
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
	return(fit0$par)
}
tSNP.P3D=function(y,X,Var,eigenG,test=NULL,LRfix=T){
 #Var: population variance components: var_e and taud	
 #fit NULL model without SNP and SNP GxE effet	
 #Xf: fixed effect (not included in GxE)
 #Xe: fixed effect (included for GxE) 	
 #test: which fixed effect to be tested. The default is to test the last one of fixed effect
 	if(is.null(test)){
 	test=n.beta	
 	}
 	if(length(test>1)){
 		if(LRfix==F){
 			stop("use LR test when there is more than one fix effect to be tested")
 		}
 	}
 	#use LR test if length.test>1
 	if(LRfix==T){
 		ln0=getLoglik(Var=Var,y,X=X[,setdiff((1:ncol(X)),test)],eigenZd,logVar=F)
 	 	ln1=getLoglik(Var=Var,y,X=X,eigenZd=eigenG,logVar=F)
 	 	Q=-2*(ln0-ln1)
 	 	p.value=pchisq(Q,df=length(test),lower.tail=F)
 	 	}else{
 	n=length(y)
	U1=eigenG$U1
	d1=eigenG$d1
	tXX=crossprod(X)
	tU1y=crossprod(U1,y)
 	tU1X=crossprod(U1,X)
 	tXy=crossprod(X,y)
 	tyy=sum(y^2)
 	outDL=getDL(var_e=var_e,taud=taud,d1=d1,n=n,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,get.tU1ehat=F)
 	beta=outDL$hat_alpha
 	vbeta=solve(outDL$tXVinvX)
 	n.beta=length(beta) 		
 	tscore=beta[test]/sqrt(vbeta[test])
 	##note: the df for t-distribution is not corrected by Satterthwaite's method. Likelihood ratio test should be better.
 	p.value=2*pt(tscore,df=n-n.beta,lower.tail=F)
 	 	}
    return(p.value)
 }
 
 ###population parameter re-estimated for each marker
 tSNP=function(y,X,eigenG){
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