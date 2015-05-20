getSetsSNPStartEnd=function(winsize,SNPstart=NULL,SNPend=NULL,chr=NULL,uniqchr=NULL,nsets=NULL,sets=NULL,returnChr=F){
	##get sets Chr information
	   if(is.null(uniqchr)|is.null(chr)){
	   if(is.null(SNPstart)| is.null(SNPend)){stop("Must specify SNPstart and SNPend when chr and uniqchr not specified")}
	   nchr=1	
	   totalSNP=SNPend-SNPstart+1
	   totalSets=ceiling(totalSNP/winsize)
	   SNPChrStart=1
	   SNPChrEnd=totalSNP
	   SetsChrStart=1
	   }else{	
	   nchr=length(uniqchr)
	   SNPChrStart=SNPChrEnd=SetsChrStart=rep(0,nchr)
	   uniqchr=sort(uniqchr)
	   SetsChrStart[1]=1
	   totalSets=0
	   for(i in 1:nchr){
	   	chr.i=uniqchr[i]
	   	whichchr=which(chr==chr.i)
		SNPChrStart[i]=min(whichchr)
	    p=length(whichchr)
	    SNPChrEnd[i]=SNPChrStart[i]+p-1
	    totalsetsPerChr=ceiling(p/winsize)
	    if(i<nchr){SetsChrStart[i+1]=SetsChrStart[i]+totalsetsPerChr}	    
	    totalSets=totalSets+totalsetsPerChr
	   }
	   }
		 ##sample sets	   
	 if(is.null(sets)){ 
	 if(is.null(nsets)){
    sets=c(1:totalSets)
     }else{
    sets=sample(1:totalSets,nsets,replace=F)
  	}
  	}
  ##get SetsSNPStart, SetsSNPEnd
	nsets=length(sets)
	SetsSNPStart=rep(0,nsets)
    SetsSNPEnd=rep(0,nsets)
	for(i in 1:nsets){
	chr.set=max(which(SetsChrStart<=sets[i]))
  	SetsSNPStart[i]=SNPChrStart[chr.set]+winsize*(sets[i]-SetsChrStart[chr.set])
  	SetsSNPEnd[i]=min(SetsSNPStart[i]+winsize,SNPChrEnd[chr.set])
	}
	out=list(sets=sets,nsets=nsets,SetsSNPStart=SetsSNPStart,SetsSNPEnd=SetsSNPEnd)
    if(returnChr==T){
    	out$SetsChrStart=SetsChrStart
    	out$SNPChrStart=SNPChrStart
    	out$SNPChrEnd=SNPChrEnd
    	out$totalSets=totalSets
    }
    return(out)	    
  }
  
simu_ug=function(eigenG,geno=NULL,nQTL=0){
  if(nQTL>0){
    bcQTLs=sample(SNPstart:SNPend,nQTL)    
    Zg=geno[,bcQTLs]
    ug=simuBeta(Z=Zg,k=kg,Type="Normal")$u
  }else{
  	ug=eigenG$U1%*%rnorm(length(eigenG$d1),mean=0,sd=sqrt(kg*var_e/mean(eigenG$d1)))
  }
  return(ug)
  }


Get_BetaMAF<-function(Type, MAF, MaxValue=1.6,Sign=0){

	n<-length(MAF)
	re<-rep(0,n)
	IDX<-which(MAF > 0)
	if(Type == "LogMAF"){
		re[IDX]<-abs(log10(MAF[IDX]))*0.4
	} else if (Type == "FixedMAF") {
		re[IDX]<-MaxValue
	}
	#Lian added this line
	re[which(re>MaxValue)]=MaxValue	

    	if(Sign > 0){
      		#temp.n<-round(n * Sign)
		temp.n<-floor(n * Sign)
		if(temp.n > 0){
      			temp.idx<-sample(1:n, temp.n)
      			re[temp.idx]<--re[temp.idx]
		}
    	} 
	return(re)
  
}

##MinMAF: minimum MAF frequency allowed for test. If MAF<MinMAF,this marker is not included in association test
simuBeta<-function(Z,k=NULL, Type="Normal", MAF=NULL, openLowerMAF=NULL,openUpperMAF=NULL,Causal.Ratio=1,Causal.MAF.Cutoff=0.03,Sign=0,scaleZ=F)
{
	if(!is.null(MAF)){
	if(!is.null(openLowerMAF) | !is.null(openUpperMAF)){
	testID=Get_testSNPs(MAF,openLowerMAF,openUpperMAF)
	Z=Z[,testID,drop=F]
	MAF=MAF[testID]
	}
	}
	beta=rep(0,ncol(Z))
	Z=meanImpute(Z)
	Z=scale(Z,T,scale=scaleZ)
    if(Type=="Normal"){
    sumvar=sum(apply(Z,2,var))
    beta=rnorm(ncol(Z),0,sqrt(k/sumvar))
    u=Z%*%beta
    	}else{
    causalID=Get_CausalSNPs(MAF,Causal.Ratio,Causal.MAF.Cutoff)
    if(length(causalID)>0){
    beta_causal=Get_BetaMAF(Type=Type,MAF=MAF[causalID],MaxValue=1.6,Sign=Sign)		
	u=Z[,causalID]%*%beta_causal
	beta[causalID]=beta_causal  
	}else {
		u=rep(0,nrow(Z))
    }
}
    return(list(Z=Z,beta=beta,u=u))
}
Get_testSNPs<-function(MAF,openLowerMAF=0,openUpperMAF=NULL){
	if(is.null(openUpperMAF)){
		IDX<-which(MAF > openLowerMAF)
		}else{
		IDX<-which(MAF > openLowerMAF & MAF < openUpperMAF)
	}
	return(IDX)
}

Get_CausalSNPs<-function(MAF, Causal.Ratio, Causal.MAF.Cutoff){

	IDX<-which(MAF < Causal.MAF.Cutoff)
	if(length(IDX) == 0){
		msg<-sprintf("No SNPs with MAF < %f",Causal.MAF.Cutoff)
		stop(msg)
	}
	

	N.causal<-round(Causal.Ratio * length(IDX))
	if(N.causal < 1){
		N.causal = 1
	}
	#print(N.causal)
	#print(Causal.Ratio)
	#print(length(IDX))
	re<-sort(sample(IDX,N.causal))
	return(re)
}
 
##simulate power and size for GxE
##plink returns G matrix as the mena variance of marker genotype. 
##if a subset of individual is selected, will eigenG work for a smaller number of individuals?
##only test the autosome SNPs
simuPower=function(geno,SNPstart,SNPend,nsets=NULL,sets=NULL,eigenG,kg=0,ks=0.2,kx=0,winsize=30,seed=1,nQTL=0,Xf,Xe,SKAT=T,Score=T,LR=F,alpha=0.05,saveAt=NULL,singleSNPtest=F,rewrite=F,GxE=T,betaType=c("Normal","LogMAF","FixedMAF")[1],MAF=NULL,Causal.MAF.Cutoff=0.03){
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
  
  if(betaType=="LogMAF"| betaType=="FixedMAF"){if(is.null(MAF)) stop("must specify MAF if betaType is LogMAF or FixedMAF")}
  
  cat(
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
  
  SetsSNPStartEnd=getSetsSNPStartEnd(winsize,p,chr=NULL,uniqchr=NULL,nsets=nsets)
  sets= SetsSNPStartEnd$sets
  nsets= SetsSNPStartEnd$nsets
 
    
  ug=simu_ug(eigenG=eigenG,geno=geno,nQTL=nQTL)
  	
  X=cbind(Xf,Xe)
  beta_x=rep(0.5/ncol(X),ncol(X))
     
  e=rnorm(N,0,1)
  y0=X%*%beta_x+ug+e
 
  if(singleSNPtest==T){
  	##fit P3D.NULL
  	tSNP.fit0=P3D.NULL(y0,X,eigenG)
  }
  
  
  for(i in 1:nsets) {
  	    seti=sets[i]
  	    set.seed(seti)
  	    cat(i,"set:", seti,"\n")
  	   
    win.start=SetsSNPStartEnd$SetsSNPStart[i]
    win.end= SetsSNPStartEnd$SetsSNPEnd[i]
    win.count=win.end-win.start+1
    Zs=geno[,win.start:win.end]	
    Zs=simuBeta(Z=Zs,k=ks,Type=betaType,MAF=MAF,Causal.MAF.Cutoff=Causal.MAF.Cutoff) 
    y=y0+Zs$u
    
    if(GxE==T){
    #GxE
    Zx=colmult(Xe,Zs$Z)
    Zx=simuBeta(Z=Zx,k=kx,Type="Normal")	
    y=y+Zx$u
    ptm=proc.time()[3]
    out=testZ(y=y,X=X,W=Zs$Z,kw=c(ncol(Zs$Z)),Zt=Zx$Z,eigenZd=eigenG,SKAT=SKAT,Score=Score,LR=LR)
    }else{
    ptm=proc.time()[3]	
    out=testZ(y=y,X=X,Zt=Zs$Z,eigenZd=eigenG,SKAT=SKAT,Score=Score,LR=LR)		
    }
    ptm2=proc.time()[3]
    if(SKAT==T){
      used.timeWindowtest=ptm2-ptm	
      cat(out$p.SKAT$p.value,"\t",file=saveAt,append=T)
      cat(i,ncol(Zs$Z),"used.timeWindowtest:", used.timeWindowtest,"\n")
    }
    if(Score==T){
      cat(out$p.Score,"\t",file=saveAt,append=T)
    }
    if(LR==T){
      cat(out$p.LR,"\t",file=saveAt,append=T)
    }
    cat("done window test\n")
    if(singleSNPtest==T){
    ##if marker is non-polymorphic, fixed effect will not work!!	
  	for(j in 1:win.count){
  	cat(j,", ")	
  	Zsj=Zs$Z[,j,drop=F]
  	if(GxE==T){
  	col.Zxj=((ncol(Xe)*(j-1)+1):(ncol(Xe)*j))
  	Zxj=Zx$Z[,col.Zxj,drop=F]
  	nZxj=ncol(Zxj)	
  	test=c((ncol(X)+ncol(Zsj))+(1:nZxj))
  	p.value=singleSNP.P3D(y,cbind(X,Zsj,Zxj),Var=tSNP.fit0,eigenG=eigenG,test=test)$p.value
  	}else{
  	test=c(ncol(X)+c(1:ncol(Zsj)))	
  	p.value=singleSNP.P3D(y,cbind(X,Zsj),Var=tSNP.fit0,eigenG=eigenG,test=test)$p.value	
  	}
  	cat(p.value,"\t",file=saveAt,append=T)
  }
  cat("\n")
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

