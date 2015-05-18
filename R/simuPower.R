##simulate power and size for GxE
##plink returns G matrix as the mena variance of marker genotype. 
##if a subset of individual is selected, will eigenG work for a smaller number of individuals?
##only test the autosome SNPs
simuPower=function(geno,SNPstart,SNPend,nsets=NULL,eigenG,var_e,kg=0,ks=0.2,kx=0,winsize=30,seed=1,nQTL=0,Xf,Xe,SKAT=T,Score=T,LR=F,alpha=0.05,saveAt=NULL,singleSNPtest=F,rewrite=F,GxE=T){
  #geno:  matrix, or gds.class object, snps in columns and individual in rows.
  #SNPstart: start position of snp to be tested
  #SNPend: end position of snp to be tested
  #nsets: number of sets for simulation
  #eigenG: eigen decoposition for G matrix
  #winsize: windowsize 
  #nQTL: number of major QTLs in the genetic background, defaul is 0
  #Xf: fixed effect (not included in GxE)
  #Xe: fixed effect (included for GxE) 
  #alpha: test size
  #singleSNPtest: whether to get p-value by single SNP test

`[.gds.class` <- 
    function(x, samples,snps)
{
       if(missing(samples)){
       	sample.id=NULL
       }else{
       	if(is.numeric(samples)){
       		sample.id=read.gdsn(index.gdsn(x,"sample.id"))[samples]
       	}else if (is.character(samples)) sample.id=samples
       	}       	
       if(missing(snps)){
       	snp.id=NULL
       }else{
        if(is.numeric(snps)){	
       	snp.id=read.gdsn(index.gdsn(x,"snp.id"))[snps]
       	}else if (is.character(j)) snp.id=snps
       }
    return(snpgdsGetGeno(x,sample.id=sample.id,snp.id=snp.id))
}  
  ##begin subsetting populations
  set.seed(seed)
  if(is.null(saveAt)){
    saveAt=paste("ks",ks,"kx",kx,".dat",sep="")
  }
  
  n.windowtest=length(which(c(SKAT,Score,LR)==T))
  if(rewrite==F){
  if(file.exists(saveAt)){
    stop(saveAt, " exists in disk, please specify a new name for save file")
  }
  }
  
  cat(
      "#var_e=",var_e,"\n",
      "#kg=",kg,"\n",
      "#ks=",ks, "\n",
      "#kx=",kx, "\n",
      "#nsets=",nsets,"\n",
      "#winsize=", winsize,"\n",
      "#nQTL=",nQTL, "\n",
      "#alpha=",alpha,"\n",
      "#n.windowtest=",n.windowtest,"\n",file=saveAt,append=F)
      
  
  p=SNPend-SNPstart+1
  N=nrow(eigenG$U1)
  tnsets=ceiling(p/winsize)
  if(is.null(nsets)){
    nsets=tnsets
    sets=c(1:nsets)
  }else{
    sets=sample(1:tnsets,nsets,replace=F)
  }
    
  if(nQTL>0){
    bcQTLs=sample(1:p,nQTL)    
    Zg=genofile[,bcQTLs]
    ug=simuBeta(Zg,kg,var_e)$u
  }else{
  	ug=eigenG$U1%*%rnorm(length(eigenG$d1),mean=0,sd=sqrt(kg*var_e/mean(eigenG$d1)))
  }
 
  
  X=cbind(Xf,Xe)
  beta_x=rep(0.5/ncol(X),ncol(X))
     
  e=rnorm(N,0,sqrt(var_e))
  y0=X%*%beta_x+ug+e
 
  if(singleSNPtest==T){
  	##fit P3D.NULL
  	tSNP.fit0=P3D.NULL(y0,X,eigenG)
  }
  
  
  for(i in 1:nsets) {
  	    seti=sets[i]
  	    cat(i,"set:", seti,"\n")
  	    ptm=proc.time()[3]
    win.start=(seti-1)*winsize+snpstart
    win.count=min(winsize,p-winsize*(seti-1))
    win.end=win.start+win.count-1
    Zs=genofile[,win.start:win.end]	
    Zs=simuBeta(Zs,ks,var_e) 
    y=y0+Zs$u
    if(GxE==T){
    #GxE
    Zx=colmult(Xe,Zs$Z)
    Zx=simuBeta(Zx,kx,var_e)	
    y=y+Zx$u
    out=testZ(y=y,X=X,W=Zs$Z,kw=c(ncol(Zs$Z)),Zt=Zx$Z,eigenZd=eigenG,SKAT=SKAT,Score=Score,LR=LR)
    }else{
    out=testZ(y=y,X=X,Zt=Zs$Z,eigenZd=eigenG,SKAT=SKAT,Score=Score,LR=LR)		
    }
    ptm2=proc.time()[3]
    if(SKAT==T){
      used.timeSKAT=ptm2-ptm	
      cat(out$p.SKAT$p.value,"\t",file=saveAt,append=T)
      cat(i,ncol(Zs$Z),"used.timeSKAT:", used.timeSKAT,"\n")
    }
    if(Score==T){
      cat(out$p.Score,"\t",file=saveAt,append=T)
    }
    if(LR==T){
      cat(out$p.LR,"\t",file=saveAt,append=T)
    }
    if(singleSNPtest==T){
    ##if marker is non-polymorphic, fixed effect will not work!!	
  	for(j in 1:win.count){
  	Zsj=Zs$Z[,j,drop=F]
  	if(GxE==T){
  	#col.Zxj=((ncol(Xe)*(j-1)+1):(ncol(Xe)*j))
  	#Zxj=Zx$Z[,col.Zxj,drop=F]
  	Zxj=colmult(Xe,Zsj)
  	nZxj=ncol(Zxj)	
  	test=c((ncol(X)+ncol(Zsj))+(1:nZxj))
  	p.value=singleSNP.P3D(y,cbind(X,Zsj,Zxj),Var=tSNP.fit0,eigenG=eigenG,test=test)$p.value
  	}else{
  	test=c(ncol(X)+c(1:ncol(Zsj)))	
  	p.value=singleSNP.P3D(y,cbind(X,Zsj),Var=tSNP.fit0,eigenG=eigenG,test=test)$p.value	
  	}
  	cat(p.value,"\t",file=saveAt,append=T)
  }
  ptm3=proc.time()[3]
  used.timeSingleSNP=ptm3-ptm2
  cat(i,ncol(Zs$Z),"used.timeSingleSNP:",used.timeSingleSNP,"\n")
 }

    
    cat("\n",file=saveAt,append=T)
    
    }
  power=list()
  pvalues=matrix(scan(saveAt,comment="#"),nrow=nsets,byrow=T)
  p.window=pvalues[,c(1:n.windowtest),drop=F]
  power.window=apply(p.window,2,function(a)mean(a<alpha)) 
  power$power.window=power.window
  if(singleSNPtest==T){
  p.singleSNP=pvalues[,-c(1:n.windowtest),drop=F]
  power.singleSNP=mean(na.omit(as.vector(p.singleSNP))<alpha)
  power$power.singleSNP=power.singleSNP
  }
  return(power)
}

