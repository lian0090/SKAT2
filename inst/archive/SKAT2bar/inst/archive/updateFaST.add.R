updateFaST.add<-function(FaST,newformula,updateFormula=T){
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
    #do nothing 
  }
  
  ##if original FaST does not contain U1, get FaST from scratch.
  if(is.null(FaST$U1)){
    eval.parent(substitute({FaST=getFaST(newformula)}))
    if(!updateFormula){
      eval.parent(subsitute(FaST[["formula"]]<-oldformula))
    }
  }else{
    addMatrix=getXyZ(as.formula(paste0("~-1+",paste0(addterms,collapse="+"))))
    addMatrix.Z=addMatrix$Z
    addMatrix.X=addMatrix$X
    rm("addMatrix")
    ##update Zw for Z
    if(!is.null(addMatrix.Z)){
      eval.parent(substitute({updateFaST.Zw(FaST,c(FaST[["listZw"]],addMatrix.Z))}))
    }
    if(!is.null(addMatrix.X)){
      eval.parent(substitute({updateFaST.X(FaST, cbind(FaST[["X"]],addMatrix.X))}))
    }	
  }
  if(updateFormula){
    eval.parent(subsitute(FaST[["formula"]]<-newformula))
  }
  
}