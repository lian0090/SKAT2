#t-test on individual markers
#It is not exactly the same as that from emma, due to the small difference in estimating variance components. 
#If I use the same variance component from emma to put into getDL, I will get exactly the same pvalue 
#Population structure previously determined. 
P3D.NULL=function(y,X0,G){
	out=list()
	if(class(G)=="matrix"){
		cat("To save computation time, it is better to supply eigenG instead of G \n");
		eigenG=getEigenG(G);
	}else if( !("eigenG" %in% class(G))){
		eigenG=G;
	}else{
		stop("G must be a matrix or eigenG")
	}
	#X0: fixed effects not being tested
	 if(any(is.na(y))){
  	#optim function will report not being able to evalue function at intial values when there is NA
  	#stop("there should be no missing values")
	whNA=which(is.na(y))
	y=y[-whNA]
	X0=X0[-whNA,,drop=F]
    eigenG$U1=eigenG$U1[-whNA,]
  	}
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
	out$y=y;
	out$X0=X0;
	out$eigenG=eigenG;
	out$Var=Var;
	class(out)=c("list","P3D0")
	return(Var)
}

##perform association mapping for provided markers while correcting for multiple test.
 GWAS.P3D=function(Xt,P3D0,multipleCorrection=T){
 	 	 if(! ("P3D0" %in% class(P3D0))) stop("P3D0 must be the model fitting from getP3D0")
 	 	 
 	 	 Xt=as.matrix(Xt)
 	   	 X0=as.matrix(P3D0$X0)
 	   	 y=P3D0$y
 	   	 Var=P3D0$Var
 	   	 eigenG=P3D0$eigenG
 	   	 if(multipleCorrection==T){
 	   	 	Me=Meff(Xt)
 	   	 	cat("effective number of test is ",Me,"\n")
 	   	 	}else{
 	   	 		Me=1
 	   	 	}
 	   	 p.value=rep(NA,ncol(Xt))
       	 
 	   	 p.value=apply(Xt,2,function(a){
 	   	 	singleSNP.P3D(y,X0=X0,Xt=a,Var=Var,eigenG=eigenG,method="LR")$p.value*Me})

 	     return(list(Me=Me,p.value=p.value))
 	   	}


singleSNP=function(y,X0,Xt,Var=NULL,G,method="LR",P3D=T,lm0=NULL){
	if(class(G)=="matrix"){
		cat("To save computation time, it is better to supply eigenG instead of G \n");
		eigenG=getEigenG(G);
	}else if( !("eigenG" %in% class(G))){
		eigenG=G;
	}else{
		stop("G must be a matrix or eigenG")
	}
	##Var is the variance for the NULL model, without fitting the test SNPs 
    X0=as.matrix(X0)
    Xt=as.matrix(Xt)
  if(any(is.na(y))){
  	#optim function will report not being able to evalue function at intial values when there is NA
  	#stop("there should be no missing values")
	whNA=which(is.na(y))
	y=y[-whNA]
	X0=X0[-whNA,,drop=F]
	Xt=Xt[-whNA,,drop=F]
   if(!is.null(eigenG)) {
   	eigenG$U1=eigenG$U1[-whNA,]
   	}
  	}
  	X=cbind(X0,Xt)
    out=list()
 	n.beta=ncol(X)		
 	test=c((ncol(X0)+1):n.beta)	
 	
 	if(length(test)>1){
 		if("t" %in% method){
 			stop("use LR test when there is more than one fix effect to be tested")
 		}
 	}
 	out$p.value=rep(0,length(method))
 	names(out$p.value)=method
    
    ##Method for not fitting genetic background
    if(is.null(eigenG)){
    	if("LR" %in% method){
    		if(is.null(lm0)){
    	  warning("it is better to supply lm0 when not fitting genetic background and using LR test ")
    	  lm0=lm(y~-1+X0) 
    	   }    	
    	}
    	lm1=lm(y~-1+X)
    	
    	for(i in 1:length(method)){
    		if(method[i]=="LR"){
    			Q=-2*(logLik(lm0)-logLik(lm1))
    			out$p.value[i]=pchisq(Q,df=length(test),lower.tail=F)
    		}
    		if(method[i]=="t"){
    			out$p.value[i]=summary(lm1)$coefficients[n.beta,4]
    		}
   	
    	}
    	
    }
    ##end methods for NULL eigenG
    
    ##Methods for fitting genetic background
    if(!is.null(eigenG)){
    #Var for NULL model
    
		if(is.null(Var)){
			if(P3D==T | "LR" %in% method){
		stop("must supply Var for NULL model if P3D is true or if using LR mothod")	
		}
		}else{
			names(Var)=c("var_e","taud")
		}
		if(P3D==T){
			Var1=Var
		}else{
	Var1=P3D.NULL(y,X0=X,eigenG)
	names(Var1)=c("var_e","taud")
		}
		   	
	for(i in 1:length(method)){
 	if("LR" == method[i]){
 		trypvalue=try({
 		ln0=getLoglik(Var=Var,y,X=X0,eigenZd=eigenG,logVar=F,REML=F)
 	 	ln1=getLoglik(Var=Var1,y,X=X,eigenZd=eigenG,logVar=F,REML=F)
 	 	Q=-2*(ln0-ln1)
 	 	p.value=pchisq(Q,df=length(test),lower.tail=F)
 	 	out$ML1=ln1
 	 	out$ML0=ln0
 	 	out$LR=Q
 	 	out$Var=Var
 	 	out$Var1=Var1
 	 	})
 	 	}else if ("t" == method[i]){
 	 	trypvalue=try({
 	    outDL=outDL=getDL.XYZ(var_e=Var1[1],taud=Var1[2],eigenZd=eigenG,X=X,y=y)
 	    beta=outDL$hat_alpha
 	    vbeta=solve(outDL$tXVinvX)
 	   tscore=abs(beta[test])/sqrt(vbeta[test,test])
 	  ##note: the df for t-distribution is not corrected by Satterthwaite's method. Likelihood ratio test should be better.
 	   p.value=2*pt(tscore,df=n-n.beta,lower.tail=F)
 	})
 	}
   if(inherits(trypvalue, "try-error")){
  		p.value=NA
  		}	
 	 out$p.value[i]=p.value
 	}		
 	}
 	
    return(out)
}


##backward compatability
singleSNP.P3D=function(y,X0,Xt,Var,eigenG,method="LR"){
 #X0: incidence matrix for fixed effects not being tested
 #Xt: incidence matrix for the marker or for the GxE  to be tested	
 #Var: population variance components: var_e and var_g 
 #methods: "t" for t-test, "LR" for likelihood ratio test, default is likelihood ratio test
 	#use LR test if length.test>1
 	out=singleSNP(y=y,X0=X0,Xt=Xt,Var=Var,eigenG=eigenG,method=method)
    return(out)
 }


 
 	
 
