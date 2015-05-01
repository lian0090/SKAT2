##functions
#Impute marker genotypes by column mean
meanImpute=function(X){
	X=apply(X,2,function(a){if(any(is.na(a))){a[which(is.na(a))]=mean(a,na.rm=T)};return(a)})
	return(X)
	}
	
#Get columwise multiplication 
colmult=function(Z1,Z2){
    if(!is.matrix(Z1) | !is.matrix(Z2)){stop("columult arguments must be matrix ")}
		p2=ncol(Z1)*ncol(Z2)
	Z3=matrix(nrow=nrow(Z1),ncol=p2)

	for(j in 1:ncol(Z2)){
	for(i in c(1:ncol(Z1))){
	Z3[,(j-1)*ncol(Z1)+i]=Z1[,i]*Z2[,j]	
		}
	}
	return(Z3)
	}	

#calculate matrix trace
tr=function(X){
	out=sum(diag(X))
	return(out)
}	

#simulate beta for random effects
simuBeta=function(Z,k,var_e=1){
	Z=meanImpute(Z)
	Z=scale(Z,T,F)
    sumvar=sum(apply(Z,2,var))
    beta=rnorm(ncol(Z),0,sqrt(k*var_e/sumvar))
    u=Z%*%beta
    return(list(Z=Z,beta=beta,u=u))
}

#get loglikelihood for Var
getLoglik=function(Var,y,X,W,kw,eigenZd,logVar,tauRel){
if(is.null(names(Var))){stop("Var must have names")}	
U1=eigenZd$U1
d1=eigenZd$d1
tU1X=crossprod(U1,X)
tU1y=crossprod(U1,y)
tXX=crossprod(X)
tXy=crossprod(X,y)
tyy=sum(y^2)
  if(!is.null(W)){
    if(is.null(kw))stop("kw must be speficied for W")
    nw=length(kw)
    if(sum(kw)!=ncol(W))stop("sum of kw should be equal to the number of columns in W")
    tauw=rep(0,nw)
    tU1W=crossprod(U1,W)
    tXW=crossprod(X,W)
    tWW=crossprod(W,W)
    tWy=crossprod(W,y)
  }else {
    nw=0
    tauw=NULL
    tU1W=NULL
    tXW=NULL
    tWW=NULL
    tWy=NULL
  }

out=neg2Log(Var=Var,tU1y=tU1y,tU1X=tU1X,tXX=tXX,tXy=tXy,tyy=tyy,d1=d1,n=n,tU1W=tU1W,tXW=tXW,tWW=tWW,tWy=tWy,kw=kw,tauRel=tauRel,logVar=logVar)
out=-1/2*out
return(out)
}

##simulate power and size for GxE
simuPower=function(geno,snp.id=NULL,N,eigenG,mu,var_e,kg=0,ks=0.2,kx=0,nsets=100,winsize=30,seed=1,nQTL=100,Xf,Xe=NULL,SKAT=T,Score=T,LR=F,alpha=0.001,saveAt=NULL){
  
  #geno:  matrix, or gds.class object
  #snp.id: names of SNPs, must be specified for gds object
  #N size of population, 
  #eigenG: eigen decoposition for G matrix
  #winsize: windowsize 
  #nQTL: number of QTLs in the genetic background
  #Xf: fixed effect (not included in GxE)
  #Xe: fixed effect (included for GxE) 
  #alpha: test size
  set.seed(seed)
  
  if(is.null(saveAt)){
    saveAt=paste("ks",ks,".dat",sep="")
  }
  
  if(file.exists(saveAt)){
    stop(saveAt, "exists in disk, please specify a new name for save file")
  }
  cat("#","mu=",mu, "var_e=",var_e, "kg=",kg,"ks=",ks, "kx=",kx, "nsets=",nsets, "winsize=", winsize,"nQTL=",nQTL, "alpha=",alpha,"\n",file=saveAt,append=F)
  
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
  
  if(is.null(nsets)){
    nsets=tnsets
    sets=c(1:nsets)
  }else{
    sets=sample(1:nsets,nsets,replace=F)
  }
  for(i in 1:nsets) {
    seti=sets[i]
    if("gds.class" %in% class(geno)){
      Zs=read.gdsn(index.gdsn(genofile, "genotype"), start=c(1,(seti-1)*winsize+1), count=c(-1,min(winsize,p-winsize*seti)))	
    }	
    if("matrix" %in% class(geno)){
      Zs=geno[,seq((seti-1)*winsize+1,length.out=min(winsize,p-winsize*seti))]
    }
    Zs_Z_u=simuBeta(Zs,ks,var_e) 
    rm("Zs")
    #GxE
    Zx=colmult(Xe,Zs_Z_u$Z)
    Zx_Z_u=simuBeta(Zx,kx,var_e)
    rm("Zx")	
    e=rnorm(N,0,sqrt(var_e))
    y=X%*%beta_x+Zg_Z_u$u+Zs_Z_u$u+Zx_Z_u$u+e
    ptm=proc.time()[3]
    out=testZ(y=y,X=X,W=cbind(Zs_Z_u$Z),kw=c(ncol(Zs_Z_u$Z)),Zt=Zx_Z_u$Z,eigenZd=eigenG,SKAT=T,Score=T,LR=F)
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
    cat("\n",file=saveAt,append=T)
    used.time=ptm2-ptm
    cat(i,used.time,"\n")
  }
  p.SKAT=matrix(scan(saveAt,comment="#"),nrow=nsets,byrow=T)
  power=apply(p.SKAT,2,mean) 
  return(power)
}
