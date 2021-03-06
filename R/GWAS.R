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


##try to make the GWAS function easier to use. In the form of NULL.form and test.form. windowby.NULL, windowby.test
# GWAS2=function(Null.form,test.form,data=NULL,setsize=20,sets=NULL,methods=NULL){

# #if mehtods is supplied, it will supersede the formula, and test the last term as specified by methods. 
  
  # #methods=c("Score", "SKAT","SSNP.P3D.LR","SSNP.P3D.t")[1]
  # #the test method for the last term is based on the 
  # #get information for GxE terms
  
  # #make sure input are formula
  # if(class(NULL.formula)!="formula"){stop("NULL.formula must be a formula")}
  # if(class(test.formula)!="formula"){stop("test.formula must be a formula")}
 	
# }

# # #GWAS on a marker windows .
GWAS= function(formula, GxE.formula, data=NULL, setsize=20, sets=NULL, methods = NULL) {
  #if mehtods is supplied, it will supersede the formula, and test the last term as specified by methods. 
  
  #methods=c("Score", "SKAT","SSNP.P3D.LR","SSNP.P3D.t")[1]
  #the test method for the last term is based on the 
  #get information for GxE terms
  
  #make sure input are formula
  if(class(formula)!="formula"){stop("formula must be a formula")}
  if(class(GxE.formula)!="formula"){stop("GxE.formula must be a formula")}
  
  containM=as.logical(attr(terms(subflags(GxE.formula)),"factors")[1,])
  term.labels=attr(terms(GxE.formula),"term.labels")
  nterms=length(term.labels)
  parsed.terms=sapply(term.labels,asCall) #parsed.terms=findterms(GxE.formula)
  deparsed.terms=term.labels #deparsed.terms=sapply(parsed.terms,deparse)
  
  movingTerm=subflags(parsed.terms[[1]]) #movingTerm=findTrueTerms(parsed.terms[[1]])[[1]]
  ##if GxE is the only term, the first element of GxE term is used for moving window. 
  if(length(movingTerm)>1){
    if(movingTerm[[1]]==as.name(":")){
      movingTerm=movingTerm[[2]]
    }
  }
  ##if formula is too long, deparse formula will break it into two elements
  #the largest width.cutoff is 500
  combinedform=as.formula(paste( paste0(deparse(formula),collapse=""),paste0(deparse(RHSForm(GxE.formula)),collapse=""),sep="+"))
  combinedform=replaceTerm(combinedform,movingTerm,1)
  last.term=parsed.terms[[nterms]]#last.term=findTrueTerms(parsed.terms[[nterms]])[[1]]
  if(is.null(methods)){
    if(length(last.term)>1 && last.term[[1]]==as.name(".R"))methods="Score" else methods="SSNP.P3D.LR"    
    #if(!is.null(findbars(parsed.terms[[nterms]]))) methods="Score" else methods="SSNP.P3D.LR"
  }
  ##update methods 
  methodsZ=intersect(methods,c("Score","SKAT"))
  n.methodsZ=length(methodsZ)
  methodsX=intersect(methods,c("SSNP.P3D.LR","SSNP.P3D.t","SSNP.LR","SSNP.t"))
  n.methodsX=length(methodsX)
  nmethods = length(methods)
  P3D.methods=intersect(methodsX,c("SSNP.P3D.LR","SSNP.P3D.t"))
  NP.methods=setdiff(methods,P3D.methods)
  NPX.methods=setdiff(methodsX,P3D.methods)
  n.P3D.methods=length(P3D.methods)
  
  
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
  #which term is in NULL model and does not contain the moving window
  NULL.NoM=setdiff(whichNoM,nterms)
  ##which term is in NULL model and contains the moving window
  NULL.M=setdiff(whichM,nterms)
  ##NULL.MR: whether the NULL model contains any random effect moving window
  if(length(NULL.M)>0){
    NULL.MR=any(sapply(parsed.terms[NULL.M],function(a)if(length(a)>0 && a[[1]]==as.name(".R"))TRUE else FALSE))
  }else NULL.MR=F
  
  ##NULL model for the part not in the moving window.   
  
  if(length(NULL.NoM)>0){     
    NULLformNoM=as.formula(paste0(deparse(formula),"+",paste0(deparsed.terms[NULL.NoM],collapse="+")))
  }else{
    NULLformNoM=formula
  }		
  
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
    n0=FaST0NoM$fr
    fr=FaST0NoM$fr
    # NULLformMi=as.formula(paste0("~.+",paste0(gsub(term.labels[1],"Mi",deparsed.terms[NULL.M],fixed=T),collapse="+")))
    NULLformMi=replaceTerm(as.formula(paste0("~.+",paste0(deparsed.terms[NULL.M],collapse="+"))),movingTerm,quote(Mi))
    ##If moving window are all fixed effects, P3D will take population parameters from without the moving windows (markers)
    if(length(P3D.methods)>0 && (!NULL.MR)){
      fit0NoM=fitNULL.FaST(FaST0NoM)
    }else fit0NoM=NULL
  }
  
  #The last term is one to be tested. if it is not in the moving window. for example M+E, E is not the moving window, therefore, the test matrix only need to be obtained once.
  if((!containM[nterms])){
    form.testMatrix=as.formula(paste0("~-1+",deparse(subflags(last.term))))
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
      if(length(NP.methods)>0|(n.P3D.methods>0 && (NULL.MR))){fit0=fitNULL.FaST(FaST0M)}
      if(n.P3D.methods>0 && (!NULL.MR) ){fit0P3D=updatefitNULLFaST(fit0NoM,FaST0M)} 
    }
    
    #get test matrix for every moving window if the test term contains the moving window. 
    if(containM[nterms]){
      #testMatrix=model.matrix(as.formula(paste0("~-1+",gsub(term.labels[1],"Mi",term.labels[nterms],fixed=T))))
      #form.testMatrix=as.formula(paste0("~-1+",deparse(replaceTerm(last.term,movingTerm, quote(Mi)))))
      form.testMatrix=as.formula(paste0("~-1+",deparse(replaceTerm(subflags(last.term),movingTerm, quote(Mi)))))
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
      if(n.P3D.methods>0 ){
        ##only when NULL model contains Moving term, and all moving term are fixed in the NULL model, fit0P3D is different from fit0
        if(length(NULL.M)>0 && (!NULL.MR)){
          for(j in 1:ncol(testMatrix)) {
            pvalue.X[1:n.P3D.methods,j]=testX(fit0P3D,Xt=testMatrix[,j],methods=P3D.methods)$p.value[P3D.methods]
          }
        }else{
          for(j in 1:ncol(testMatrix)) {
            pvalue.X[1:n.P3D.methods,j]=testX(fit0,Xt=testMatrix[,j],methods=P3D.methods)$p.value[P3D.methods]
          }
        }
      }
      if(length(NPX.methods)>0) {
        for(j in 1:ncol(testMatrix)) {
          pvalue.X[(n.P3D.methods+1):n.methodsX,j]=testX(fit0,Xt=testMatrix[,j],methods=NPX.methods)$p.value[NPX.methods]
        }
      }
      out[i,(n.methodsZ+1):(n.methodsZ+n.methodsX)]=apply(pvalue.X,1,min)*Me
    }        
  }
  
  #return(out)
  outGWAS=list(p.value=out,fit0=fit0)
  class(outGWAS)="GWAS"
  return(outGWAS)
}




