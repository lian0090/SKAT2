 print.GWAS<-function(GWASobj){
  cat("$p.value\n")
  print(head(GWASobj$p.value))
  cat("\n")
  cat("$fit0\n")
  print(GWASobj$fit0)
}


print.lmm<-function(fit0){
  cat("Call:\n")
  print(substitute(fitNULL(foo),list(foo=fit0$FaST$formula)))
  cat("Variance components for the NULL model:\n")
  print(fit0$Var) 
}

# # #GWAS on a marker windows .
GWAS= function(formula, GxE.formula, data=NULL, setsize=20, sets=NULL, methods = NULL) {
  #if mehtods is supplied, it will supersede the formula, and test the last term as specified by methods. 
  
  #methods=c("Score", "SKAT","SSNP.P3D.LR","SSNP.P3D.t")[1]
  #the test method for the last term is based on the 
  #get information for GxE terms
  
  containM=as.logical(attr(terms(subbars(GxE.formula)),"factors")[1,])
  term.labels=attr(terms(subbars(GxE.formula)),"term.labels")
  nterms=length(term.labels)
  parsed.terms=findterms(GxE.formula)
  deparsed.terms=sapply(parsed.terms,deparse)

  movingTerm=findTrueTerms(parsed.terms[[1]])[[1]]
  ##if GxE is the only term, the first element of GxE term is used for moving window. 
  if(length(movingTerm)>1){
  if(movingTerm[[1]]==as.name(":")){
    movingTerm=movingTerm[[2]]
  }
  }
  ##if formula is too long, deparse formula will break it into two elements
  combinedform=as.formula(paste( paste0(deparse(formula),collapse=""),paste0(deparse(RHSForm(GxE.formula)),collapse=""),sep="+"))
  combinedform=replaceTerm(combinedform,movingTerm,1)
  last.term=findTrueTerms(parsed.terms[[nterms]])[[1]]
  if(is.null(methods)){
    if(!is.null(findbars(parsed.terms[[nterms]]))) methods="Score" else methods="SSNP.P3D.LR"
  }
  ##update methods 
  methodsZ=intersect(methods,c("Score","SKAT"))
  n.methodsZ=length(methodsZ)
  methodsX=intersect(methods,c("SSNP.P3D.LR","SSNP.P3D.t","SSNP.LR","SSNP.t"))
  n.methodsX=length(methodsX)
  nmethods = length(methods)
 
  
  #update sets   
  if(is.null(sets))eval(substitute(sets<-rep(c(1:ceiling(ncol(M)/setsize)), each = setsize)[1:ncol(M)],list(M=movingTerm)))
  nSets = length(unique(sets))
  #container for output         
  out = matrix(nrow = nSets, ncol = nmethods)
  colnames(out)=c(methodsZ,methodsX)
  
  ##start fitting NULL model or get FaST         
  if(nterms==0){
    stop("must supply test terms\n")
  }
  whichM=which(containM)
  whichNoM=which(!containM)
  formM=formNoM=NULL
  NULL.NoM=setdiff(whichNoM,nterms)
  NULL.M=setdiff(whichM,nterms)
  
  ##NULL model for the part not in the moving window.   
  
  if(length(NULL.NoM)>0){     
    NULLformNoM=as.formula(paste0(deparse(formula),"+",paste0(deparsed.terms[NULL.NoM],collapse="+")))
  }else{
    NULLformNoM=formula
  }		
  envG=new.env()
 
  ##if no terms in the moving window is in the NULL model, we only need to fit a single NULL model.
  if(length(NULL.M)==0){
    fit0=fitNULL(NULLformNoM,data=data)
    whichNa=fit0$FaST$whichNa
    n0=fit0$FaST$n0
    fr=fit0$FaST$fr
  }else{
    #every set needs its only NULL model.
    FaST0NoM=getFaST(NULLformNoM,data=data)
    whichNa=FaST0NoM$whichNa
    n0=FaST0NoM$n0
    fr=fit0$fr
   # NULLformMi=as.formula(paste0("~.+",paste0(gsub(term.labels[1],"Mi",deparsed.terms[NULL.M],fixed=T),collapse="+")))
    NULLformMi=replaceTerm(as.formula(paste0("~.+",paste0(deparsed.terms[NULL.M],collapse="+"))),movingTerm,quote(Mi))
  }
  
 #The last term is one to be tested. if it is not in the moving window. for example M+E, E is not the moving window, therefore, the test matrix only need to be obtained once.
 if((!containM[nterms])){
   form.testMatrix=as.formula(paste0("~-1+",term.labels[nterms]))
   testMatrix=model.matrix(form.testMatrix,data=model.frame(form.testMatrix,data=data,na.action=na.pass))
   testMatrix=checkVar(testMatrix)
   testMatrix=checkNA(testMatrix,n0=n0,whichNa=whichNa)
 }
 
 
  
  for (i in 1:nSets) {
    cat("set",i,"\n")
    if(setsize>1){whichSeti=which(sets==i)}else whichSeti=i
    Mi=as.matrix(eval(substitute(M[,whichSeti],list(M=movingTerm))))
    Mi=meanImpute(Mi)
    data$Mi=Mi
    ##should add methods to check the rank of Mi
    if(length(NULL.M)>0){
      ##replace the label of moving window by Mi.
      ##this can be optimized by direclty add, substract,or replace on FaST0NoM to save memory
      FaST0M=FaST.add(FaST0NoM,newformula=NULLformMi,data=data)
      fit0=fitNULL.FaST(FaST0M)
    }
    
    #get test matrix for every moving window if the test term contains the moving window. 
    if(containM[nterms]){
      #testMatrix=model.matrix(as.formula(paste0("~-1+",gsub(term.labels[1],"Mi",term.labels[nterms],fixed=T))))
      form.testMatrix=as.formula(paste0("~-1+",deparse(replaceTerm(last.term,movingTerm, quote(Mi)))))
      testMatrix=model.matrix(form.testMatrix,data=model.frame(form.testMatrix,data=data,na.action=na.pass))
      if(length(whichNa)>0)testMatrix=testMatrix[-whichNa,]
      testMatrix=checkVar(testMatrix)
    }
    
    
    ##extract sets from the first term
    if(n.methodsZ>0){
      out[i,1:n.methodsZ]=testZ(fit0,Zt=testMatrix,methods=methodsZ)$p.value[methodsZ]
    }
   
    if(n.methodsX>0){
      Me=Meff(testMatrix)
      pvalue.X=matrix(nrow=n.methodsX,ncol=ncol(testMatrix))
      for(j in 1:ncol(testMatrix)) pvalue.X[,j]=testX(fit0,Xt=testMatrix[,j],methods=methodsX)$p.value[methodsX]
      out[i,(n.methodsZ+1):(n.methodsZ+n.methodsX)]=apply(pvalue.X,1,min)*Me
    }        
  }
  
  #return(out)
  outGWAS=list(p.value=out,fit0=fit0)
  class(outGWAS)="GWAS"
  return(outGWAS)
}


testZ = function(fit0, Zt, methods = "Score") {
  
  out = vector("list",length=4)
  names(out)=c("Score","sdScore","p.value","Q")
  out$p.value=rep(NA,length=2)
  names(out$p.value)=c("Score","SKAT")
  
  updateFaST.Zt(fit0$FaST, Zt)
  if(is.null(fit0$FaST$Zt)){
    return(out)
  }
  
  if (class(fit0) ==c("lm")) {
    X=fit0$FaST$X
    resid = residuals(fit0)
    s2 = summary(fit0)$sigma^2
    Q = sum(crossprod(resid, Zt)^2)/s2/2
    tXZt = crossprod(X, Zt)
    Zw.1 = (crossprod(Zt) - t(tXZt) %*% solve(crossprod(X)) %*% tXZt)/2
    lambda = eigen(Zw.1, symmetric = TRUE, only.values = TRUE)$values
    lambda1 = lambda
    IDX1 <- which(lambda >= 0)
    # eigenvalue bigger than mean(lambda1[IDX1])/100000 
    IDX2 <- which(lambda1 > mean(lambda1[IDX1])/1e+05)
    lambda <- lambda1[IDX2]
    
    if ("SKAT" %in% methods) {
      out$p.value['SKAT'] <- Get_PValue.Lambda(lambda, Q)$p.value
      out$Q = Q
    }
    if ("Score" %in% methods) {
      varS = (sum(Zw.1^2)/(s2^2))/2
      S = (Q - (sum(diag(Zw.1))))/s2
      out$Score = S
      out$sdScore = sqrt(varS)
      out$p.value['Score'] = pnorm(S, mean = 0, sd = sqrt(varS), lower.tail = F)
    }
    return(out) 
  } else if(class(fit0)== "lmm"){
    #logVar, paramterize variance components with log when using REML to restrict variance component to be larger than 0.
    ##SKAT test or LR test
    getQ = ("SKAT" %in% methods)
    getS = ("Score" %in% methods)
    Qdis = getDL(fit0$Var, fit0$FaST, getQ = getQ, getS = getS, getNeg2Log = F, REML = T)
    
    if ("SKAT" %in% methods) {
      Q = Qdis$Q
      lambda = Qdis$lambda
      lambda1 = lambda
      IDX1 <- which(lambda >= 0)
      #eigenvalue bigger than mean(lambda1[IDX1])/100000 
      IDX2 <- which(lambda1 > mean(lambda1[IDX1])/1e+05)
      lambda <- lambda1[IDX2]
      p.value <- Get_PValue.Lambda(lambda, Q)$p.value
      out$p.value['SKAT'] = p.value
      out$Q = Q
    }
    #Score test
    if ("Score" %in% methods) {
      S = Qdis$S
      sdS = Qdis$sdS
      out$Score = S
      out$sdScore = sdS
      out$p.value['Score'] = pnorm(S/sdS, lower.tail = F)
    }
  }
    
  return(out)
}






