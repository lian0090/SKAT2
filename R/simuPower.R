#Get columwise multiplication
#if ncol(Z2)=p, there are p chunks of Z1*Z2[,j] 
colmult=function(Z1,Z2){
  #return Z3, Z3=(Z1*Z2[,1],Z1*Z2[,2],..Z1*Z2[,ncol(Z2)])
  if(!is.matrix(Z1) | !is.matrix(Z2)){stop("columult arguments must be matrix ")}
  p2=ncol(Z1)*ncol(Z2)
  Z=matrix(nrow=nrow(Z1),ncol=p2)
  colID.Z1=rep(0,ncol(Z))
  colID.Z2=rep(0,ncol(Z))
  for(j in 1:ncol(Z2)){  
    for(i in c(1:ncol(Z1))){
      col.ij=	(j-1)*ncol(Z1)+i
      Z[,col.ij]=Z1[,i]*Z2[,j]	
      colID.Z1[col.ij]=i
      colID.Z2[col.ij]=j
    }
  }
  
  return(list(Z=Z,colID.Z1=colID.Z1,colID.Z2=colID.Z2))
}	

getSetsSNPID<-
  function(snp.id,setsSNPnames){
    ##SetsSNPnames: a list of SNPnames for sets
    setsSNPID=list()
    nsets=length(setsSNPnames)
    for(i in 1:nsets){
      setsSNPID[[i]]=which(snp.id%in%setsSNPnames[[i]])
    }
    return(setsSNPID)	    
  }


##return the index for snps in each sets
getSetsSNPID.StartEnd<-
  function(winsize,SNPstart=NULL,SNPend=NULL,chr=NULL,testchr=NULL,nsets=NULL){
    ##get sets Chr information
    if(is.null(chr) | is.null(testchr)){
      if(is.null(SNPstart)| is.null(SNPend)){stop("Must specify SNPstart and SNPend when chr and testchr not specified")}
      nchr=1	
      totalSNP=SNPend-SNPstart+1
      totalSets=ceiling(totalSNP/winsize)
      SNPChrStart=1
      SNPChrEnd=totalSNP
      SetsChrStart=1
    }else{	
      nchr=length(testchr)
      SNPChrStart=SNPChrEnd=SetsChrStart=rep(0,nchr)
      testchr=sort(testchr)
      SetsChrStart[1]=1
      totalSets=0
      for(i in 1:nchr){
        chr.i=testchr[i]
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
    if(is.null(nsets)){
      setsIDX=c(1:totalSets)
    }else{
      if(nsets>totalSets){
        setsIDX=c(1:totalSets)
      }else{	
        setsIDX=sample(1:totalSets,nsets,replace=F)
      }
    }
    
    nsets=length(setsIDX)
    setsSNPID=list()
    for(i in 1:nsets){
      chr.set=max(which(SetsChrStart<=setsIDX[i]))
      SetsSNPStarti=SNPChrStart[chr.set]+winsize*(setsIDX[i]-SetsChrStart[chr.set])
      SetsSNPEndi=min(SetsSNPStarti+winsize-1,SNPChrEnd[chr.set])
      setsSNPID[[i]]=SetsSNPStarti:SetsSNPEndi
    }
    return(setsSNPID)	    
  }



simu_ug.QTL=function(SNPstart,SNPend,geno,nQTL,kg=1,Sign){
  
  bcQTLs=sample(SNPstart:SNPend,nQTL)    
  Zg=geno[,bcQTLs]
  ug=simuBeta(Z=Zg,k=kg,Type="Normal",Sign=Sign)$u  	 
  return(ug)
}
simu_ug.eigenG=function(eigenG,kg=1){
  ug=eigenG$U1%*%rnorm(length(eigenG$d1),mean=0,sd=sqrt(kg^2/mean(eigenG$d1)))
  return(ug)
}

Get_BetaMAF<-function(Type, MAF, MaxValue=1.6){
  if(!(Type %in% c("LogMAF","FixedMAF"))){
    stop("Get_BetaMAF: Type must be in  LogMAF, FixedMAF")
  }
  n<-length(MAF)
  re<-rep(0,n)
  IDX<-which(MAF > 0)
  if(Type == "LogMAF"){
    re[IDX]<-abs(log10(MAF[IDX]))*MaxValue/4
  } else if (Type == "FixedMAF") {
    re[IDX]<-MaxValue
  }
  #Lian added this line
  re[which(re>MaxValue)]=MaxValue	
  return(re)
  
}

Get_CausalSNPs<-
  function(MAF, p,Type=c("Normal","LogMAF","FixedMAF","Equal")[2],Causal.Ratio, Causal.MAF.Cutoff){
    if(Type=="Normal" | Type =="Equal"){
      IDX=c(1:p)
    }else{
      IDX<-which(MAF < Causal.MAF.Cutoff)
    }
    if(length(IDX) == 0){
      msg<-sprintf("No SNPs with MAF < %f",Causal.MAF.Cutoff)
      return(NULL)
    }
    N.causal<-round(Causal.Ratio * length(IDX))
    if(N.causal < 1){
      N.causal = 1
    }
    #print(N.causal)
    #print(Causal.Ratio)
    #print(length(IDX))
    if(length(IDX)>1){
      ##if length(IDX)=1, sample(IDX) is acutally sampleing 1:IDX
      re<-sort(sample(IDX,N.causal))
    }else{
      re=IDX
    }
    return(re)
  }

Get_testSNPsMAF<-function(MAF,openLowerTestMAF=0,openUpperTestMAF=NULL){
  if(!is.null(openUpperTestMAF)){
    IDXu<-which(MAF < openUpperTestMAF)
  }else{
    IDXu=1:length(MAF)
  }
  if(!is.null(openLowerTestMAF)){
    IDXl<-which(MAF > openLowerTestMAF)
  }else{
    IDXl=1:length(MAF)
  }
  IDX=sort(intersect(IDXl,IDXu))
  if(length(IDX)==0){return(NULL)}
  return(IDX)
}



##MinMAF: minimum MAF frequency allowed for test. If MAF<MinMAF,this marker is not included in association test
simuBeta<-function(Z,k=NULL, Type="Normal", MAF=NULL,causalID=NULL,Causal.Ratio=1,Causal.MAF.Cutoff=0.03,Sign=0,MaxValue=1.6,scaleZ=T)
{
  if(!(Type %in% c("Normal","LogMAF","FixedMAF","Equal"))){
    stop("simuBeta: Type must be in Normal, LogMAF, FixedMAF, Equal")
  }
  p=ncol(Z)
  if(is.null(causalID)){
    causalID=Get_CausalSNPs(MAF=MAF,p=p,Type=Type,Causal.Ratio=Causal.Ratio,Causal.MAF.Cutoff=Causal.MAF.Cutoff)
  }
  n.causal=length(causalID)
  beta=rep(0,ncol(Z))
  Z=meanImpute(Z)
  Z=scale(Z,T,scale=scaleZ)
  if(n.causal>0){
    Z_causal=Z[,causalID,drop=F]
    if(Type=="Normal"){
      ##if you change sumvar to meanvar, the relative performance of singleSNPtest and Score test might change	
      #sumvar=sum(apply(Z_causal,2,var))
      #beta_causal=abs(rnorm(ncol(Z_causal),0,sqrt(k/sumvar)))
      beta_causal=abs(rnorm(ncol(Z_causal),0,k))
    }
    if (Type=="Equal"){
      #sumvar=sum(apply(Z_causal,2,var))
      #beta_causal=rep(sqrt(k/sumvar),ncol(Z_causal))
      beta_causal=rep(k,ncol(Z_causal))
    }
    if (Type=="LogMAF" | Type=="FixedMAF"){  				
      beta_causal=Get_BetaMAF(Type=Type,MAF=MAF[causalID],MaxValue=MaxValue)		
    }
    if(Sign > 0){
      temp.n<-floor(n.causal * Sign)
      if(temp.n > 0){
        temp.idx<-sample(1:n.causal, temp.n)
        beta_causal[temp.idx]<--beta_causal[temp.idx]
      }
    } 
    
    u=Z_causal%*%beta_causal
    beta[causalID]=beta_causal  
  }else{
    u=rep(0,nrow(Z))
  }
  return(list(Z=Z,beta=beta,u=u,causalID=causalID))
}


##simulate power and size for GxE
##plink returns G matrix as the mena variance of marker genotype. 
##if a subset of individual is selected, will eigenG work for a smaller number of individuals?
##only test the autosome SNPs
simuPower=function(geno,SNPstart=NULL,SNPend=NULL,chr=NULL,testchr=NULL,nsets=NULL,setsSNPID=NULL,eigenG=NULL,kg=1,ks=0.2,kx=0.1,winsize=20,seed=1,nQTL=0,Xf,Xe,windowtest=c("SKAT","Score",NULL)[1],saveAt=NULL,singleSNPtest=c("LR","t")[1],GxE=c(NA,"Normal","LogMAF","FixedMAF","Multiply","Equal")[4],betaType=c("Normal","LogMAF","FixedMAF","Equal")[1],MAF=NULL,Causal.MAF.Cutoff=0.03,openLowerTestMAF=NULL,openUpperTestMAF=NULL,Sign=0,removeZtFromG=T){
  #geno:  matrix, or gds.class object, snps in columns and individual in rows.
  #SNPstart: start position of snp to be tested
  #SNPend: end position of snp to be tested
  #nsets: number of sets for simulation
  #setsSNPID: A list of integer index for SNPs in each sets.
  #eigenG: eigen decoposition for G matrix
  #winsize: windowsize 
  #nQTL: number of major QTLs in the genetic background, defaul is 0
  #Xf: fixed effect (not included in GxE)
  #Xe: fixed effect (included for GxE) 
  #windowtest
  #singleSNPtest: LR, t-test, or NULL
  ##begin subsetting populations
  set.seed(seed)
  if(is.null(saveAt)){
    saveAt=paste("ks",ks,"kx",kx,sep="")
  }
  saveAt.Windowtest=file.path(saveAt,"Windowtest.dat")
  saveAt.SingleSNPtest=file.path(saveAt,paste("SingleSNP_",singleSNPtest,".dat",sep=""))
  saveAt.betaZs=file.path(saveAt,"betaZs.dat")
  saveAt.betaZx=file.path(saveAt,"betaZx.dat")
  if(!file.exists(saveAt))dir.create(saveAt,recursive=T)
  
  n.windowtest=length(windowtest)
  n.singleSNPtest=length(singleSNPtest)
  
  if(betaType=="LogMAF"| betaType=="FixedMAF"){if(is.null(MAF)) stop("must specify MAF if betaType is LogMAF or FixedMAF")}
  
  
  cat(
    "#kg=",kg,"\n",
    "#ks=",ks, "\n",
    "#kx=",kx, "\n",
    "#nsets=",nsets,"\n",
    "#winsize=", winsize,"\n",
    "#nQTL=",nQTL, "\n",
    "#n.windowtest=",n.windowtest,"\n",
    "#saveAt is ",saveAt,"\n")
  	file.create(saveAt.Windowtest,F)
  	file.create(saveAt.betaZs,F)
  	file.create(saveAt.betaZx,F)
  	for(k in 1:n.singleSNPtest){
    	file.create(saveAt.SingleSNPtest[k],F)
  	}
  N=nrow(Xf)
  
  if(is.null(setsSNPID)){
  	#get the index of SNPs in each sets
    setsSNPID=getSetsSNPID.StartEnd(winsize,SNPstart=SNPstart,SNPend=SNPend,chr=chr,testchr=testchr,nsets=nsets)
  }
  nsets=length(setsSNPID)
  nSNPs=sum(sapply(setsSNPID,length))
  #store results:should keep the initial values NA, so that for the markers not tested, this is still NA
  pvalue.window=matrix(NA,n.windowtest,nSNPs)
  pvalue.SingleSNP=matrix(NA,n.singleSNPtest,nSNPs)
 
 ##fixed effects
  X=cbind(Xf,Xe)
  beta_x=rep(0.5/ncol(X),ncol(X)) 
  beta_xe=beta_x[(ncol(Xf)+1):ncol(X)] 
 ##residuals 
  e=rnorm(N,0,1)
  
  ##simulate genetic background
  if(!is.null(eigenG)){
  	if(nQTL>0){
    	ug=simu_ug.QTL(SNPstart,SNPend,geno,nQTL,kg=kg,Sign=Sign)
  	}else{
    	ug=simu_ug.eigenG(eigenG,kg=kg)	
  	}
    y0=X%*%beta_x+ug+e  
  }else{
  	y0=X%*%beta_x+e
  }
  
  #commented out on June 25,2015. P3D.NULL should fit the true phenotypes instead of y0
  # if(!is.null(singleSNPtest)){
  # 	##fit P3D.NULL
  # 	tSNP.fit0=P3D.NULL(y0,X,eigenG)
  # }  
  ##begin simulating each window
  setStarti=1
  for(i in 1:nsets) {
    set.seed(i)
    cat("set:", i,"\n")
    setsSNPIDi=setsSNPID[[i]]
    Zs=geno[,setsSNPIDi]
    MAFi=MAF[setsSNPIDi]
    Zs=simuBeta(Z=Zs,k=ks,Type=betaType,MAF=MAFi,Causal.MAF.Cutoff=Causal.MAF.Cutoff,MaxValue=1.6,Sign=Sign) 
    if(betaType %in% c("LogMAF","FixedMAF")){
      testID.Zs=Get_testSNPsMAF(MAF=MAFi,openLowerTestMAF=openLowerTestMAF,openUpperTestMAF=openUpperTestMAF)
    }else{
      testID.Zs=1:ncol(Zs$Z)	
    }
    p.Zs=length(setsSNPIDi)
    p.testZs=length(testID.Zs)
    cat(Zs$beta,"\n",file=saveAt.betaZs,append=T)
    y=y0+Zs$u
    
    ###simulate GxE
    if (!is.null (GxE)){
      #GxE
      Zx=colmult(Xe,Zs$Z)
      Zx.colID.Xe=Zx$colID.Z1
      Zx.colID.Zs=Zx$colID.Z2
      Zx=Zx$Z
      p.Zx=ncol(Zx)
      causalID.Zx=which(Zx.colID.Zs %in% Zs$causalID)
      
      if(GxE %in% c("Normal","LogMAF","FixedMAF","Equal")){
        Zx=simuBeta(Z=Zx,k=kx,Type=GxE,MAF=MAFi[Zx.colID.Zs],causalID=causalID.Zx,MaxValue=0.8,Sign=Sign)	
        
      } else if(GxE == "Multiply"){
        Zx=list(Z=scale(Zx,T,F))
        Zx$beta=beta_xe[Zx.colID.Xe] * Zs$beta[Zx.colID.Zs]
        Zx$u=Zx$Z%*%Zx$beta	
        Zx$causalID=causalID.Zx
      }else{
        stop("GxE Type is not correct")
      }
      cat(Zx$beta,"\n",file=saveAt.betaZx,append=T)
      y=y+Zx$u
      if(p.testZs>0){
        testID.Zx=which(Zx.colID.Zs %in% testID.Zs)
        p.testZx=length(testID.Zx)
      }
    setCount=p.Zx #number of single tests in a window
    }else{
    	setCount=p.Zs
    }
    	    
    #end simulating GxE
    
    ##do test for each window
    #if p.testZs=0, all the tests become NA
    
    if(p.testZs==0){
      #this should be changed, if Zs is not tested, this does not mean it is not causal.  
      cat(rep(NA,n.windowtest),"\n",file=saveAt.Windowtest,append=T)	
      if(!is.null(GxE)){
        for(k in 1:n.singleSNPtest) {cat(rep(NA,p.Zs*ncol(Xe)),"\n",file=saveAt.SingleSNPtest[k],append=T)}
      }
      else{
        for(k in 1:n.singleSNPtest){cat(rep(NA,p.Zs),"\n",file=saveAt.SingleSNPtest[k],append=T)}
      }
    }
    
    
    if(p.testZs>0){
      #####start windowtest (always output both Score and SKAT test)
      if(!is.null(windowtest)){
        if(is.null(GxE)){
          ptm=proc.time()[3]	
          out=testWindow(y,X=X,Zt=Zs$Z[,testID.Zs,drop=F],eigenG=eigenG,removeZtFromG=removeZtFromG)
               }else{
          ptm=proc.time()[3]
          ##For GxE, you should not removeZtFromG, because when you remove the interaction matrix from G, you also removed the incidence
          #matrix for the fixed effects. The p-value in this case is 0. But I do not know why the p-value should be 0.
          out=testWindow(y,X=X,W=list(Zs$Z[,testID.Zs,drop=F]),Zt=Zx$Z[,testID.Zx,drop=F],eigenG=eigenG,removeZtFromG=F)
        }
        ptm2=proc.time()[3]
        
        used.timeWindowtest=ptm2-ptm	
        cat("used.timeWindowtest:", used.timeWindowtest,"\n")
        for(wi in 1:n.windowtest){
        	if(windowtest[wi]=="SKAT"){
        			pvalue.window.wi=out$p.SKAT$p.value
        	}
        	if(windowtest[wi]=="Score"){
        			pvalue.window.wi=out$p.Score
        	}
          cat(pvalue.window.wi," ",file=saveAt.Windowtest,append=T)  
          pvalue.window[wi,setStarti:(setStarti+setCount-1)]=pvalue.window.wi
        }
          cat("done window test\n")
        
      }
      #####end windowtest
      
      #####start single SNP test
      if(!is.null(singleSNPtest)){
        
        if(!is.null(eigenG)){
        	tSNP.fit0Var=P3D.NULL(y=y,X0=X,eigenG)$Var
        }else{
        	tSNP.fit0Var=NULL
        }
        
        ##if marker is non-polymorphic, fixed effect will not work!!	
        pmt.Me1=proc.time()[3]
        if(is.null(GxE)){
          Me=Meff(Zs$Z[,testID.Zs,drop=F])
        }else{
          Me=Meff(Zx$Z[,testID.Zx,drop=F])	
        }
        pmt.Me2=proc.time()[3]
        cat("used.time calculating Me:",pmt.Me2-pmt.Me1,"\n")
        
        for(j in 1:p.Zs){
          if (j %in% testID.Zs){
            #cat(setsSNPIDi[j],", ")	
            Zsj=Zs$Z[,j,drop=F]
            if(!is.null(GxE)){
              col.Zxj=which(Zx.colID.Zs==j)
              Zxj=Zx$Z[,col.Zxj,drop=F]
              nZxj=ncol(Zxj)	
              test=c((ncol(X)+ncol(Zsj))+(1:nZxj))
              p.value=singleSNP(y,X0=cbind(X,Zsj),Xt=Zxj,Var=tSNP.fit0Var,eigenG=eigenG,method=singleSNPtest)$p.value *Me
            }else{
              test=c(ncol(X)+c(1:ncol(Zsj)))	
              p.value=singleSNP(y,X0=X,Xt=Zsj,Var=tSNP.fit0Var,eigenG=eigenG,method=singleSNPtest)$p.value *Me	
            }
          }else p.value=rep(NA,n.singleSNPtest)
          
          for(k in 1:n.singleSNPtest){
            cat(p.value[k]," ",file=saveAt.SingleSNPtest[k],append=T)
            pvalue.SingleSNP[k,setStarti+j-1]=p.value[k]
            
          }
        }
        cat("Me",Me,"\n")
        cat("\n")
        ptm3=proc.time()[3]
        used.timeSingleSNP=ptm3-ptm2
        cat("used.timeSingleSNP:",used.timeSingleSNP,"\n")
      }
      
      for(k in 1:n.singleSNPtest){
        cat("\n",file=saveAt.SingleSNPtest[k],append=T)
      }
    }
    ##update setStarti
    if(is.null(GxE)){
     	setStarti=setStarti+p.Zs
     	}else{
     setStarti=setStarti+p.Zx
     }
  } 
  #results=list(pvalue.SingleSNP=pvalue.SingleSNP,pvalue.window=pvalue.window)
  #saveRDS(results,file=file.path(saveAt,"pvalues.rds"))  
  out=get_results(saveAt=saveAt,windowtest=c("SKAT","Score"),singleSNPtest=c("LR"))
  saveRDS(out,file=file.path(saveAt,"results.rds"))
  if(!is.null (GxE)){
  	poweri=getPower(p.window=out$p.window,p.singleSNP=out$p.singleSNP,beta=out$beta.Zx,alpha=0.05)
  	}else{
  	poweri=getPower(p.window=out$p.window,p.singleSNP=out$p.singleSNP,beta=out$beta.Zs,alpha=0.05)	
  		
  	}
  	
  saveRDS(poweri,file=file.path(saveAt,"poweri.rds"))
  return(poweri)
}

##do not put function argument as na.strings=na.strings in function definition
readLinesListf=function(con,split="\\s+",na.strings){
    dat=lapply(strsplit(readLines(con),split=split),function(a){
      whichNa=which(a %in% na.strings)
      if(length(whichNa)>0)a[whichNa]=NA;
      a=as.numeric(a);
      return(a)
    })
    return(dat)
  }

get_results=function(saveAt,windowtest,singleSNPtest,na.strings="NA"){
  out=list()
  out$p.singleSNP=vector("list",length(singleSNPtest))
  names(out$p.singleSNP)=singleSNPtest
  n.windowtest=length(windowtest)	
  n.singleSNPtest=length(singleSNPtest)
  
  saveAt.Windowtest=file.path(saveAt,"Windowtest.dat")
  saveAt.SingleSNPtest=file.path(saveAt,paste("SingleSNP_",singleSNPtest,".dat",sep=""))
  saveAt.betaZs=file.path(saveAt,"betaZs.dat")
  saveAt.betaZx=file.path(saveAt,"betaZx.dat")
  
  p.window=scan(saveAt.Windowtest,comment="#",na.strings=na.strings)
  out$p.window=matrix(p.window,ncol=n.windowtest,byrow=T)
  colnames(out$p.window)=windowtest
  for(i in 1:n.singleSNPtest){
    out$p.singleSNP[[i]]=readLinesListf(saveAt.SingleSNPtest[i],na.strings=na.strings)
  }
  out$beta.Zs=readLinesListf(saveAt.betaZs,na.strings=na.strings)
  out$beta.Zx=readLinesListf(saveAt.betaZx,na.strings=na.strings)
  return(out)
}


getPower.singleSNP=function(p.singleSNP,whNon0,wh0,alpha){
  n.0=length(which(wh0))
  n.Non0=length(which(whNon0))
  out=rep(0,2)
  whnotNA.singleSNP=sapply(p.singleSNP,function(a)!all(is.na(a)))
  whNon0.singleSNP=which(whNon0 & whnotNA.singleSNP)
  wh0.singleSNP=which(wh0 & whnotNA.singleSNP)
  if(length(whNon0.singleSNP)>0){
    p.singleSNP.power=sapply(p.singleSNP[whNon0.singleSNP],function(a) min(a,na.rm=T))	
    power.singleSNP=length(which(na.omit(unlist(p.singleSNP.power))<alpha))/n.Non0
    # cat("n.Non0",n.Non0,"\n")
    out[1]=power.singleSNP 
  }	  
  if(length(wh0.singleSNP)>0){
    #bonferroni correction
    p.singleSNP.size=sapply(p.singleSNP[wh0.singleSNP],function(a) min(a,na.rm=T))
    size.singleSNP=length(which(na.omit(unlist(p.singleSNP.size))<alpha))/n.0
    out[2]=size.singleSNP 
  }
  names(out)=c("power","size")
  return(out)
}

getPower.window=function(p.window,wh0,whNon0,alpha){
  n.window=ncol(p.window)
  n.0=length(which(wh0))
  n.Non0=length(which(whNon0))
  whnotNA.window=apply(p.window,1,function(a)!all(is.na(a)))
  whNon0.window=which(whNon0 & whnotNA.window)
  wh0.window=which(wh0 & whnotNA.window)
  out=rep(NA,n.window*2)
  namesout=NULL
  if(length(whNon0.window)>0){
    power.window=apply(p.window[whNon0.window,,drop=F],2,function(a)length(which(a<alpha))/n.Non0) 	
    out[1:n.window]=power.window
  }
  namesout=c(namesout,paste("power_",colnames(p.window),sep=""))
  #bonferroni correction
  if(length(wh0.window)>0){
    size.window=apply(p.window[wh0.window,,drop=F],2,function(a)length(which(a<alpha))/n.0) 
    out[(n.window+1):(2*n.window)]=size.window
  }
  namesout=c(namesout,paste("size_",colnames(p.window),sep=""))
  names(out)=namesout
  cat(out,"\n")
  return(out)
}

getPower=function(p.window,p.singleSNP,beta,alpha){
  n.singleSNP=length(p.singleSNP)
  nobs.window=nrow(p.window)
  nobs.singleSNP=min(sapply(p.singleSNP,length))
  nobs.beta=length(beta)
  nobs.min=min(c(nobs.window,nobs.singleSNP,nobs.beta))
  if(nobs.min < nobs.window){
    p.window=p.window[1:nobs.min,]
  }
  if(nobs.min < nobs.singleSNP){
    p.singleSNP=lapply(p.singleSNP,function(a) return(a[1:nobs.min]))
  }
  if(nobs.min < nobs.beta){
    beta=beta[1:nobs.min]
  }
  wh0=sapply(beta,function(a)all(a==0))
  whNon0=!wh0
  out=NULL
  namesout=NULL
  out1=getPower.window(p.window,wh0,whNon0,alpha)
  out=c(out,out1)
  namesout=c(namesout,names(out1))
  for(k in 1:n.singleSNP){
    outk=getPower.singleSNP(p.singleSNP[[k]],whNon0,wh0,alpha)
    namesoutk=paste(names(outk),"_singleSNP_",names(p.singleSNP)[k],sep="")
    out=c(out,outk)
    namesout=c(namesout,namesoutk)  
  }
  names(out)=namesout
  return(out)
}



