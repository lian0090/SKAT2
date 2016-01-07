##return environments without dollar signs for variable names
getFaST=function(formula,data=NULL,getG=T,Ztest=NULL,...){
  ##... further parameters passed to model.frame such as na.action
  noGform=formula
  RHS.noGform=noCallPattern(RHSForm(formula),call.pattern=c(as.name(".eigenG"),as.name(".G")))
  if(is.null(RHS.noGform))RHSForm(noGform)<-1 else RHSForm(noGform)<-RHS.noGform
  environment(noGform)=environment(formula)
  if (is.null(data)) {
    fr=model.frame(subflags(noGform),...)
  } else {
    fr=model.frame(subflags(noGform),data=data,...)
  }
  
  ##get a whole model.frame before you get X, y, Z to avoid environment contamination.
  #if(is.null(environment(formula)$n0)){
  whichNa=as.integer(attr(fr,"na.action"))
  n=nrow(fr)
  n0=n+length(whichNa)
  y <- model.response(fr)
  
  ##get fixed term
  fixedform <- noGform
  nb <- noCallPattern(RHSForm(fixedform),call.pattern=as.name(".R"))
  
  if(is.null(nb)){
    RHSForm(fixedform)=1
  }else {
    RHSForm(fixedform)=nb;  
  }
  X=model.matrix(fixedform,fr)
  if(ncol(X)==0)X=NULL
  
  
  #get random term
  method="lm"
  out=list()
  out$Z=NULL
  Z=NULL
  eigenG=NULL
  listZw=NULL
  listG=NULL
  whichEigenG=NULL
  n.randomTerm=0
  namesRandomTerm=NULL
  if(!getG){
    bars=findCallPattern(RHSForm(formula),".R")
    if(!is.null(bars)){
      Z=lapply(bars,FUN=function(x) model.matrix(as.formula(substitute(~-1+foo,list(foo=x[[2]]))),fr))
      namesRandomTerm= unlist(lapply(bars, function(x) deparse(x[[2]])))
      names(Z) <-namesRandomTerm
      n.randomTerm=length(Z)
      }
    out=list(Z=Z,X=X,y=y,whichNa=whichNa,n0=n0,eigenG=eigenG,listZw=listZw,fr=fr,listG=listG,method=method,whichFaST=whichEigenG,n.randomTerm=n.randomTerm,namesRandomTerm=namesRandomTerm)
  }else{
    #default eigenG, listZw are NULL, default method is "lm"
    randomTerm=findCallPattern(formula,c(".R",".G",".eigenG"))
    if(is.null(randomTerm)){
      out=list(Z=Z,X=X,y=y,whichNa=whichNa,n0=n0,eigenG=eigenG,listZw=listZw,fr=fr,listG=listG,method=method,whichFaST=whichEigenG,n.randomTerm=n.randomTerm,namesRandomTerm=namesRandomTerm)
    }else{
      n.randomTerm=length(randomTerm)
      listRandomTerm=vector(mode="list",length=n.randomTerm)
      termFlag=rep(NA,n.randomTerm)
      namesRandomTerm=rep(NA,n.randomTerm)
      kw=rep(0,n.randomTerm)
      for(i in 1:n.randomTerm){
        listRandomTerm[[i]]=list(Z=NULL,G=NULL,eigenG=NULL)
        termFlag[i]=deparse(randomTerm[[i]][[1]])
        trueTerm=randomTerm[[i]][[2]]
        namesRandomTerm[i]=deparse(trueTerm)
        if(termFlag[i]==".G"){
          G<-eval(trueTerm)
          if(length(whichNa)>0)G<-G[-whichNa,-whichNa]
          eigenG=getEigenG(G=G)
          kw[i]=length(eigenG$d1)
          listRandomTerm[[i]]$G=G
          listRandomTerm[[i]]$eigenG=eigenG
        }
        if(termFlag[i]==".eigenG"){
          eigenG<-eval(trueTerm)
          if(length(whichNa)>0)eigenG<-getEigenG(Zg=sweep_prod(eigenG$U1[-whichNa,],sqrt(eigenG$d1),F,2))
          kw[i]=length(eigenG$d1)
          listRandomTerm[[i]]$eigenG=eigenG
        }
        if(termFlag[i]==".R"){
          Z=model.matrix(as.formula(substitute(~-1+trueTerm,list(trueTerm=trueTerm))),fr)
          eigenG<-getEigenG(Zg=Z)
          kw[i]=length(eigenG$d1)
          listRandomTerm[[i]]$Z=Z
          listRandomTerm[[i]]$eigenG=eigenG
        }
      }
      ##re-set eigenG to NULL  
      eigenG=NULL
      whichEigenG=which.max(kw)
      kwT=sum(kw[-whichEigenG])
      if(kwT>=n){
        method="brute"
        listG=vector(mode="list",length=n.randomTerm)
        names(listG)=namesRandomTerm
        for(i in 1:n.randomTerm){
          if(is.null(listRandomTerm[[i]]$G)){
            listG[[i]]=with(listRandomTerm[[i]]$eigenG,tcrossprod(sweep_prod(U1,sqrt(d1),F,2)))
          }else{
            listG[[i]]=listRandomTerm[[i]]$G
          }
        }
        out=list(Z=Z,X=X,y=y,whichNa=whichNa,n0=n0,eigenG=eigenG,listZw=listZw,fr=fr,listG=listG,method=method,whichFaST=whichEigenG,n.randomTerm=n.randomTerm,namesRandomTerm=namesRandomTerm)
        
        
      
        }else{
        method="FaST"
        eigenG=listRandomTerm[[whichEigenG]]$eigenG
        attr(eigenG,"termName")=namesRandomTerm[whichEigenG]
        listRandomTerm=listRandomTerm[-whichEigenG]
        nw=n.randomTerm-1
        if(nw>0){
          listZw=vector(mode="list",nw)
          for(i in 1:nw){
            listZw[[i]]=with(listRandomTerm[[i]]$eigenG,sweep_prod(U1,sqrt(d1),F,2))
          }
          names(listZw)=namesRandomTerm[-whichEigenG]
        }
        out=list(Z=Z,X=X,y=y,whichNa=whichNa,n0=n0,eigenG=eigenG,listZw=listZw,fr=fr,listG=listG,method=method,whichFaST=whichEigenG,n.randomTerm=n.randomTerm,namesRandomTerm=namesRandomTerm)
        out$U1 = eigenG$U1
        out$d1 = eigenG$d1
        out$nd = length(out$d1)
        out$n=length(out$y)
        out$tU1y = crossprod(out$U1, out$y)
        if (out$nd < out$n) {
          out$tyy = sum(out$y^2)
        }
        #must be outside of if!(is.null(Z)), otherwise, X and Ztest will not be updated into FaST
        updateFaST.X(out, X)
        updateFaST.Zw(out, listZw)
        updateFaST.Ztest(out, Ztest)
      }      
    }         
  } 
  
  out$formula=formula
  class(out) = c("FaST")
  return(out)
  
}

asCall=function(a){
  a=as.formula(paste("~",a))
  a=RHSForm(a)
  return(a)
}


RHSForm <- function(formula) {
  formula[[length(formula)]]
}

`RHSForm<-` <- function(formula,value) {
  formula[[length(formula)]] <- value
  formula
}




subflags=function(term){
  if(is.name(term)||!is.language(term))return(term)
  if(length(term)==2){
    if(is.call(term)&&term[[1]] == as.name(".G"))
      term=term[[2]]
    if(is.call(term)&&term[[1]]==as.name(".eigenG"))
      term=term[[2]]
    if(is.call(term)&&term[[1]]==as.name(".R"))
      term=term[[2]]
  }
  if(length(term)>=2){
    for(j in 2:length(term)) term[[j]]=subflags(term[[j]])
  }
  return(term)
  
}



###findterms as orignal form ()
findterms<- function(term) {
  if (term[[1]] == as.name("~")) 
    term = term[[length(term)]]
  
  fft = function(term) {
    #give a right hand side term
    if(length(term)>1){
      if (term[[1]] == as.name("+")) {
        return(c(fft(term[[2]]), fft(term[[3]])))
      } else {
        
        return(term)
      }
    }else{
      if(is.atomic(term))return(NULL) else return(term)
    }
  }
  
  
  term = fft(term)
  if (!is.list(term)) {
    term = list(term)
  }
  return(term)
}




replaceTerm=function(term, pattern=NULL, replace=quote(Mi)){
  
  ft<-function(term){
    if(length(term)>1){
      if(length(term)==length(pattern)){if(term==pattern) return(replace)}
      for(i in 1:length(term)){
        if(length(term[[i]])!=length(pattern)) {term[[i]]=ft(term[[i]])} else{
          if(term[[i]]==pattern){
            term[[i]]=replace
          }else term[[i]]=ft(term[[i]])
        }
      }
      return(term)
    }else{
      if(length(term)!=length(pattern)){return(term)}else{
        if(term==pattern){
          return(replace)
        }else return(term)
      }
    }
  }
  return(ft(term))
  
}


findTrueTerms=function(term){
  if(typeof(term)!="language" & typeof(term)!="symbol")(stop("term must be language"))
  if(length(term)>1){
    if(term[[1]]==as.name("~"))term=term[[length(term)]]
  }
  
  ft<-function(term){
    if(length(term)>1){
      if(term[[1]]==as.name("+")) return(c(ft(term[[2]]),ft(term[[3]])))
      if(term[[1]]==as.name("(")) return(ft(term[[2]]))
      if(term[[1]]==as.name("|")) return(ft(term[[3]]))
      if(term[[1]]==as.name(".eigenG"))return(ft(term[[2]]))
      if(term[[1]]==as.name(".G"))return(ft(term[[2]]))
      if(term[[1]]==as.name(".R"))return(ft(term[[2]]))
      return(term)
    }
    if(is.atomic(term))return(NULL)
    if(is.name(term))return(term)
  }
  term = ft(term)
  if (!is.list(term)) {
    term = list(term)
  }
  return(term)
  
}

#remove call pattern from formula,if one argument of a Call argument goes to NULL, this Call also diappears
#~.G(a) becomes NULL
#y~.G(a) becomes y
#y~.G(a)+b becomes y~.G(a)
noCallPattern=function(term,call.pattern=c(".eigenG",".G")){
  call.pattern=sapply(call.pattern,as.name)
  if (!any(sapply(call.pattern,function(a)deparse(a)%in%all.names(term)))) 
    return(term)
  if (is.call(term) && any(sapply(call.pattern,function(a)term[[1]]==a))) return(NULL) 
  if (length(term) == 2) {
    nb <- noCallPattern(term[[2]],call.pattern)
    if (is.null(nb)) 
      return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- noCallPattern(term[[2]],call.pattern)
  nb3 <- noCallPattern(term[[3]],call.pattern)
  if (is.null(nb2)) 
    return(nb3)
  if (is.null(nb3)) 
    return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}

##return a list of items that has call.pattern
findCallPattern=function(term,call.pattern){
  call.pattern=sapply(call.pattern,as.name)
  #call.pattern=as.name(".eigenG")
  fc=function(term,call.pattern){
    if (is.name(term) || !is.language(term)) 
      return(NULL)
    if (any(sapply(call.pattern,function(a)term[[1]]==a)))
      return(term)
    if (length(term) == 2) 
      return(fc(term[[2]],call.pattern))
    return(c(fc(term[[2]],call.pattern), fc(term[[3]],call.pattern)))
  }
  term=fc(term,call.pattern)
  if(!is.list(term) & !is.null(term))term=list(term)
  return(term)
}


LHSForm <- function(formula) {
  if (length(formula)==2) NULL else formula[[2]]
}








#currently only add more terms to Fast
# should be able to replace and subtract in the future
##do not update Ztest.

FaST.add<-function(FaST,newformula,data=NULL){
  #... further arguments passed to data.frame.
  oldformula=FaST$formula
  newformula=update.formula(oldformula,newformula)
  oldterms=findterms(oldformula)
  newterms=findterms(newformula)
  addterms=setdiff(newterms,oldterms)
  dropterms=setdiff(oldterms,newterms)
  
  if(length(dropterms)>0){
    stop("currently do not support drop terms from the oldformula")
  }
  if(length(addterms)==0){
    return(FaST)
  }
  if(class(data)=="model.frame"){
    #warning("data is model.frame, only pure names in the formula can be properly dealt\n")
  }
  ##if original FaST does not contain U1, get FaST from scratch.
  if(is.null(FaST$U1)){
    n0=FaST$n0
    whichNa=FaST$whichNa
    FaST=getFaST(newformula,data=data)
    FaST$n0=n0
    FaST$whichNa=whichNa
  }else{
    add.formula=as.formula(paste0("~-1+",paste0(addterms,collapse="+")))
    environment(add.formula)=environment(newformula)
    addMatrix=checkData(add.formula,data=data,getG=F,na.action=na.pass)
    addMatrix.Z=addMatrix$Z
    addMatrix.X=addMatrix$X
    rm("addMatrix")
    ##update Zw for Z
    if(!is.null(addMatrix.Z)){
      updateFaST.Zw(FaST,c(FaST[["listZw"]],addMatrix.Z))
    }
    if(!is.null(addMatrix.X)){
      updateFaST.X(FaST, cbind(FaST[["X"]],addMatrix.X))
    }	
  }
  
  FaST[["formula"]]<-newformula
  return(FaST)
}


updateFaST.Zw = function(FaST, listZw) {
  
  eval.parent(substitute({
    FaST[["Zw"]] = NULL
    FaST[["kw"]] = NULL
    FaST[["nw"]] = 0
    FaST[["tU1W"]] = NULL
    FaST[["tXW"]] = NULL
    FaST[["tWW"]] = NULL
    FaST[["tWy"]] = NULL
    FaST[["tWZt"]] = NULL
  }))
  
  if (!is.null(listZw)) {
    
    if(nrow(listZw[[1]])==FaST$n0){
      if (length(FaST$whichNa) > 0) {
        listZw = lapply(listZw, function(a) {
          a[-whichNa, , drop = F]
        })
      }
    }else if (nrow(listZw[[1]]) != FaST$n) {
      stop("Zw must have the same number of rows as the original y or the current y")
    }
    kw = sapply(listZw, ncol)
    nw = length(kw)
    Zw = do.call(cbind, listZw)
    eval.parent(substitute({
      FaST[["listZw"]] = listZw
      FaST[["Zw"]] = Zw
      FaST[["kw"]] = kw
      FaST[["nw"]] = nw
      if (!is.null(FaST$U1)) FaST[["tU1W"]] = crossprod(FaST$U1, Zw)
      if (FaST$nd < FaST$n) {
        FaST[["tXW"]] = crossprod(FaST$X, Zw)
        FaST[["tWW"]] = crossprod(Zw)
        FaST[["tWy"]] = crossprod(Zw, FaST$y)
        if (!is.null(FaST$Ztest)) FaST[["tWZt"]] = crossprod(Zw, FaST$Ztest)
      }
    }))
  }else{
    FaST[["listZw"]]=NULL
    #change FaST["listZw"] last to avoid changes to listZw if listZw is an element of FaST.
  }
}


updateFaST.X = function(FaST, X) {
  if (is.null(X)) {
    X = matrix(rep(1, FaST$n))
  } else {
    X=as.matrix(X)
    if(nrow(X)==FaST$n0){
      if (length(FaST$whichNa) > 0) {
        X = as.matrix(X[-whichNa, , drop = F])
      }
    }else {
      if (nrow(X) != FaST$n) {
        stop("X must have the same number of rows as the original y or the current y")
      }
    }
  }
  eval.parent(substitute({
    FaST[["X"]] = X
    FaST[["tU1X"]] = NULL
    FaST[["tXX"]] = NULL
    FaST[["tXy"]] = NULL
    FaST[["tXW"]] = NULL
    FaST[["tXZt"]] = NULL
  }))
  if (FaST$method=="FaST"){   
  eval.parent(substitute(
   {
      FaST[["tU1X"]] = crossprod(FaST$U1, X)
      if (FaST$nd < FaST$n) {
        FaST[["tXX"]] = crossprod(X)
        FaST[["tXy"]] = crossprod(X, FaST$y)
        if (!is.null(FaST$Zw)) FaST[["tXW"]] = crossprod(X, FaST$Zw)
        if (!is.null(FaST$Ztest)) FaST[["tXZt"]] = crossprod(X, FaST$Ztest)
      }
    }
  ))
  }
}

updateFaST.Ztest = function(FaST, Ztest) {
  eval.parent(substitute({
    FaST[["tU1Zt"]] = NULL
    FaST[["tyZt"]] = NULL
    FaST[["tXZt"]] = NULL
    FaST[["tZtZt"]] = NULL
    FaST[["tWZt"]] = NULL
  }))
  if (!is.null(Ztest)) {
    
    if(nrow(Ztest)==FaST$n0){if(length(FaST$whichNa)>0)Ztest=Ztest[-FaST$whichNa,]}else{
      if (nrow(Ztest) != FaST$n) {
        stop("Ztest must have the same number of rows as the original y or the current y")
      }
    }
    varXt = apply(Ztest, 2, function(a) var(a, na.rm = T))
    whichVar = which(varXt > 0)
    if (length(whichVar) > 0) {
      Ztest = Ztest[, whichVar, drop = F]
    } else {
      warning("No variants in Ztest")
      Ztest = NULL
    }
    eval.parent(substitute({
      FaST[["Ztest"]] = Ztest
      if (!is.null(Ztest)) {
        if (!is.null(FaST$U1)) {
          FaST$tU1Zt = crossprod(FaST$U1, Ztest)
          if (FaST$nd < FaST$n) {
            FaST[["tyZt"]] = crossprod(FaST$y, Ztest)
            FaST[["tXZt"]] = crossprod(FaST$X, Ztest)
            FaST[["tZtZt"]] = crossprod(Ztest)
            if (!is.null(FaST$Zw)) {
              FaST[["tWZt"]] = crossprod(FaST$Zw, Ztest)
            }
          }
        }
      }
    }))
    
  }else{
    FaST[["Ztest"]] = NULL  ##do this lage, if Ztest is from FaST$Ztest, then setting FaST[["Zt"]] to NULL also will set Ztest to NULL
  }
}


checkNA=function(X,n0,whichNa){
  X=as.matrix(X)
  if(nrow(X)!=n0)eval(substitute(stop("in checkNA number of rows of X must be equal to n0")))
  if(length(whichNa)>0){
    if (length(whichNa) > 0) 
      X = X[-whichNa,,drop=F ]
  }  
  X=meanImpute(X)
  return(X)
}
checkVar=function(X){
  X=as.matrix(X)
  varX = apply(X, 2, var)
  which.noVar = which(varX == 0)
  n.noVar=length(which.noVar)
  if(n.noVar==ncol(X)){
    eval(substitute(warning(paste("No variants in",X))))
    return(NULL)
  }
  if (n.noVar > 0) {
    X = X[, -which.noVar, drop = F]
  } 
  return(X)	
}








