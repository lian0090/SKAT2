testZ = function(fit0, Ztest, methods = "Score") {
  
  out = vector("list",length=4)
  names(out)=c("Score","sdScore","p.value","Q")
  out$p.value=rep(NA,length=2)
  names(out$p.value)=c("Score","SKAT")
  
  updateFaST.Ztest(fit0$FaST, Ztest)
  if(is.null(fit0$FaST$Ztest)){
    return(out)
  }
  
  if (class(fit0) ==c("lm")) {
    X=fit0$FaST$X
    Ztest=fit0$FaST$Ztest
    resid = residuals(fit0)
    s2 = summary(fit0)$sigma^2
    Q = sum(crossprod(resid, Ztest)^2)/s2/2
    tXZt = crossprod(X, Ztest)
    tZtPZts2 = (crossprod(Ztest) - t(tXZt) %*% solve(crossprod(X)) %*% tXZt)
    lambda = eigen(tZtPZts2 , symmetric = TRUE, only.values = TRUE)$values/2    
    if ("SKAT" %in% methods) {
      out$p.value['SKAT'] <- Get_PValue.Lambda(lambda, Q)$p.value
      out$Q = Q
    }
    if ("Score" %in% methods) {
      varS = (sum(tZtPZts2 ^2)/(s2^2))/2
      S = (Q - (sum(diag(tZtPZts2)))/2)/s2
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




