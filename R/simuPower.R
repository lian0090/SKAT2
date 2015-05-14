##simulate power and size for GxE
simuPower.GxE=function(geno,snp.id=NULL,N,eigenG,mu,var_e,kg=0,ks=0.2,kx=0,nsets=100,sets=NULL,winsize=30,seed=1,nQTL=100,Xf,Xe=NULL,SKAT=T,Score=T,LR=F,alpha=0.001,saveAt=NULL,singleSNPtest=F,rewrite=F,GxE=T){
  
  #geno:  matrix, or gds.class object
  #snp.id: names of SNPs, must be specified for gds object
  #N size of population, 
  #eigenG: eigen decoposition for G matrix
  #winsize: windowsize 
  #nQTL: number of QTLs in the genetic background
  #Xf: fixed effect (not included in GxE)
  #Xe: fixed effect (included for GxE) 
  #alpha: test size
  #singleSNPtest: whether to get p-value by single SNP test
  set.seed(seed)
  
  if(is.null(saveAt)){
    saveAt=paste("ks",ks,"kx",kx,".dat",sep="")
  }
  
  
  if(rewrite==F){
  if(file.exists(saveAt)){
    stop(saveAt, " exists in disk, please specify a new name for save file")
  }
  }
  cat("# mu =",mu,"\n",
      "#var_e=",var_e,"\n",
      "#kg=",kg,"\n",
      "#ks=",ks, "\n",
      "#kx=",kx, "\n",
      "#nsets=",nsets,"\n",
      "#winsize=", winsize,"\n",
      "#nQTL=",nQTL, "\n",
      "#alpha=",alpha,"\n",file=saveAt,append=F)
  
  if("gds.class" %in% class(geno)){
    if(! "gdsfmt" %in% rownames(installed.packages())) stop("must install gdsfmt package to use gds.class file")
    if(is.null(snp.id)){stop("must specify snp.id for gds object")}
    p=length(snp.id)
    bcQTLs=sample(1:p,100)
    Zgsnps=snp.id[bcQTLs]
    Zg=snpgdsGetGeno(genofile,sample.id=NULL,snp.id=Zgsnps,snpfirstdim=F)
    
  }
  if("matrix" %in% class(geno)){
    p=ncol(geno)
    bcQTLs=sample(1:p,100)
    Zg=geno[,bcQTLs]
  }
  X=cbind(Xf,Xe)
  beta_x=rep(0.5/ncol(X),ncol(X))
  beta_g=kg*var_e
  Zg_Z_u=simuBeta(Zg,kg,var_e)
  rm("Zg")
  
  tnsets=ceiling(p/winsize)
  if(is.null(sets)){
  if(is.null(nsets)){
    nsets=tnsets
    sets=c(1:nsets)
  }else{
    sets=sample(1:tnsets,nsets,replace=F)
  }
  }else{nsets=length(sets)}
  
  e=rnorm(N,0,sqrt(var_e))
  y0=X%*%beta_x+e
 
  if(singleSNPtest==T){
  	##fit P3D.NULL
  	tSNP.fit0=P3D.NULL(y0,X,eigenG)
  }
  
  
  for(i in 1:nsets) {
  	    seti=sets[i]
  	    cat(i,"set:", seti,"\n")
  	    ptm=proc.time()[3]
    win.start=(seti-1)*winsize+1
    win.count=min(winsize,p-winsize*(seti-1))
    win.end=win.start+win.count-1
    if("gds.class" %in% class(geno)){
      Zs=read.gdsn(index.gdsn(genofile, "genotype"), start=c(1,win.start), count=c(-1,win.count))	
    }	
    if("matrix" %in% class(geno)){
      Zs=geno[,seq((seti-1)*winsize+1,length.out=win.count)]
    }
    Zs_Z_u=simuBeta(Zs,ks,var_e) 
    rm("Zs")
    y=y0+Zg_Z_u$u+Zs_Z_u$u
    if(GxE==T){
    #GxE
    Zx=colmult(Xe,Zs_Z_u$Z)
    Zx_Z_u=simuBeta(Zx,kx,var_e)
    rm("Zx")	
    y=y+Zx_Z_u$u
    out=testZ(y=y,X=X,W=Zs_Z_u$Z,kw=c(ncol(Zs_Z_u$Z)),Zt=Zx_Z_u$Z,eigenZd=eigenG,SKAT=SKAT,Score=Score,LR=LR)
    }else{
    out=testZ(y=y,X=X,Zt=Zx_Z_u$Z,eigenZd=eigenG,SKAT=SKAT,Score=Score,LR=LR)		
    }
    ptm2=proc.time()[3]
    if(SKAT==T){
      cat(out$p.SKAT$p.value,"\t",file=saveAt,append=T)
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
  	Zsj=Zs_Z_u$Z[,j,drop=F]
  	if(GxE==T){
  	col.Zxj=((ncol(Xe)*(j-1)+1):(ncol(Xe)*j))
  	Zxj=Zx_Z_u$Z[,col.Zxj,drop=F]
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
  }

    
    cat("\n",file=saveAt,append=T)
    used.timeSKAT=ptm2-ptm
    used.timeSingleSNP=ptm3-ptm2
    cat(i,ncol(Zx_Z_u$Z),used.timeSKAT,used.timeSingleSNP,"\n")
    }
  
  p.SKAT=matrix(scan(saveAt,comment="#"),nrow=nsets,byrow=T)
  power=apply(p.SKAT,2,function(a)mean(a>alpha)) 
  return(power)
}

