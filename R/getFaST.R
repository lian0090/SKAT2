##return environments without dollar signs for variable names
checkData=function(formula,data,getG=F,...){
  ##... further parameters passed to model.frame such as na.action
  noGform=noCallPattern(subbars(formula),call.pattern=c(as.name(".eigenG"),as.name(".G")))
  noGform=subbars(noGform)  
  environment(noGform)=environment(formula)
  
  if (is.null(data)) {
    fr=model.frame(noGform,...)
  } else {
    fr=model.frame(noGform,data=data,...)
  }
  #names(fr)=gsub(".*\\$(.*)$","\\1",names(fr))
  #denv=list2env(fr)  
  # denv=parent.env(environment())
  # formula=as.formula(formula,env=denv)
 
  ##get a whole model.frame before you get X, y, Z to avoid environment contamination.
  if(is.null(environment(formula)$n0)){
    whichNa=as.integer(attr(fr,"na.action"))
    n0=nrow(fr)+length(whichNa)
  }else{
    whichNa=environment(formula)$whichNa
    n0=environment(formula)$n0
  }
  
  ##get fixed term
  fixedform <- noGform
  nb <- nobars(RHSForm(fixedform))
  
  if(is.null(nb)){
    RHSForm(fixedform)=1
  }else {
    RHSForm(fixedform)=nb;  
  }
  X=model.matrix(fixedform,fr)
  if(ncol(X)==0)X=NULL
  
  #get random term
  bars=findbars(RHSForm(formula))
  if(is.null(bars))Z=NULL else{
    Z=lapply(bars,FUN=function(x) model.matrix(as.formula(substitute(~-1+foo,list(foo=x[[3]]))),fr))
    names(Z) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
  }
  
  y <- model.response(fr)
  

  ##G and eigenG must be dealt with separately 
  if(getG){
    parsedterms0=findterms(formula)
    whichEigenG=which(sapply(parsedterms0,function(a)ifelse(is.call(a),a[[1]]==quote(.eigenG),a==quote(.eigenG))))
    whichG=which(sapply(parsedterms0,function(a)ifelse(is.call(a),a[[1]]==quote(.G),a==quote(.G))))
    if(length(whichEigenG)>0){
     ##this recomputes eigenG. There are more efficient ways of updating eigenG. Do it in a later time 
      eigenG<-eval(parsedterms0[[whichEigenG]][[2]])
      if(class(eigenG)!="eigenG")stop(".eigenG() must be used on an eigenG object ")
      if(length(whichNa)>0)eigenG<-getEigenG(Zg=sweep_prod(eigenG$U1[-whichNa,],sqrt(eigenG$d1),F,2))
      listZw=Z
    }else if(length(whichG)>0){
      G<-eval(parsedterms0[[whichG]][[2]])
      if(length(whichNa)>0)G<-G[-whichNa,-whichNa]
      eigenG=getEigenG(G=G)
      listZw=Z
    }else if(!is.null(Z)){
        ncolZ=sapply(Z,function(a){ncola=ncol(a);ifelse(is.null(ncola),1,ncola)})
        whichMaxcol=which.max(ncolZ)
        whichEigenG=ifelse(whichMaxcol,whichMaxcol,1)  
        eigenG=getEigenG(Zg=Z[[whichEigenG]])
        if (length(Z) > 1) {
          listZw = Z[-whichEigenG]
        }else listZw=NULL
    }else{
        eigenG=NULL
        listZw=NULL
      }
  }else{
    eigenG=NULL
    listZw=NULL
  }


   return(list(Z=Z,X=X,y=y,whichNa=whichNa,n0=n0,eigenG=eigenG,listZw=listZw,fr=fr))


}


asCall=function(a){
  a=as.formula(paste("~",a))
  a=RHSForm(a)
  return(a)
}
expandDoubleVerts <- function(term)
{
    expandDoubleVert <- function(term) {
	frml <- formula(paste0("~", deparse(term[[2]])))
	## need term.labels not all.vars to capture interactions too:
	newtrms <- paste0("0+", attr(terms(frml), "term.labels"))
	if(attr(terms(frml), "intercept")!=0)
	    newtrms <- c("1", newtrms)
	as.formula(paste("~(",
			 paste(vapply(newtrms, function(trm)
				      paste0(trm, "|", deparse(term[[3]])), ""),
			       collapse=")+("), ")"))[[2]]
    }

    if (!is.name(term) && is.language(term)) {
	if (term[[1]] == as.name("(")) {
	    term[[2]] <- expandDoubleVerts(term[[2]])
	}
	stopifnot(is.call(term))
	if (term[[1]] == as.name('||'))
	    return( expandDoubleVert(term) )
	## else :
	term[[2]] <- expandDoubleVerts(term[[2]])
	if (length(term) != 2) {
	    if(length(term) == 3)
		term[[3]] <- expandDoubleVerts(term[[3]])
	}
    }
    term
}

RHSForm <- function(formula) {
    formula[[length(formula)]]
}

`RHSForm<-` <- function(formula,value) {
    formula[[length(formula)]] <- value
    formula
}
 findbars<-
function (term) 
{
    fb <- function(term) {
        if (is.name(term) || !is.language(term)) 
            return(NULL)
        if (term[[1]] == as.name("(")) 
            return(fb(term[[2]]))
        stopifnot(is.call(term))
        if (term[[1]] == as.name("|")) 
            return(term)
        if (length(term) == 2) 
            return(fb(term[[2]]))
        c(fb(term[[2]]), fb(term[[3]]))
    }
    expandSlash <- function(bb) {
        makeInteraction <- function(x) {
            if (length(x) < 2) 
                return(x)
            trm1 <- makeInteraction(x[[1]])
            trm11 <- if (is.list(trm1)) 
                trm1[[1]]
            else trm1
            list(substitute(foo:bar, list(foo = x[[2]], bar = trm11)), 
                trm1)
        }
        slashTerms <- function(x) {
            if (!("/" %in% all.names(x))) 
                return(x)
            if (x[[1]] != as.name("/")) 
                stop("unparseable formula for grouping factor")
            list(slashTerms(x[[2]]), slashTerms(x[[3]]))
        }
        if (!is.list(bb)) 
            expandSlash(list(bb))
        else unlist(lapply(bb, function(x) {
            if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]]))) 
                lapply(unlist(makeInteraction(trms)), function(trm) substitute(foo | 
                  bar, list(foo = x[[2]], bar = trm)))
            else x
        }))
    }
    modterm <- expandDoubleVerts(if (is(term, "formula")) 
        term[[length(term)]]
    else term)
    expandSlash(fb(modterm))
}


 nobars<-
function (term) 
{
    if (!any(c("|", "||") %in% all.names(term))) 
        return(term)
    if (is.call(term) && term[[1]] == as.name("|")) 
        return(NULL)
    if (is.call(term) && term[[1]] == as.name("||")) 
        return(NULL)
    if (length(term) == 2) {
        nb <- nobars(term[[2]])
        if (is.null(nb)) 
            return(NULL)
        term[[2]] <- nb
        return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) 
        return(nb3)
    if (is.null(nb3)) 
        return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}
subbars <- function(term)
{
    if (is.name(term) || !is.language(term)) return(term)
    if (length(term) == 2) {
        term[[2]] <- subbars(term[[2]])
        return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name('|'))
        term[[1]] <- as.name('+')
    if (is.call(term) && term[[1]] == as.name('||'))
        term[[1]] <- as.name('+')
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
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
removeDollar=function(term){
  ft<-function(term){
    if(length(term)>1){
      if(term[[1]]==as.name("$")) term=ft(term[[3]]) else {
        if(length(term)==2)term[[2]]=ft(term[[2]])
        if(length(term)==3){
          term[[2]]=ft(term[[2]])
          term[[3]]=ft(term[[3]])
        }
      }
      
    }
      return(term)
   
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

#find .eigenG 
noCallPattern=function(term,call.pattern=c(as.name(".eigenG"),as.name(".G"))){
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


LHSForm <- function(formula) {
  if (length(formula)==2) NULL else formula[[2]]
}







getFaST = function(formula,data=NULL,Zt = NULL) {
  out=checkData(formula,data,getG=T)
  eigenG<-out$eigenG
  if(!is.null(eigenG)){
    ##need more careful check to deal with Na
		out$U1 = eigenG$U1
		out$d1 = eigenG$d1
		out$nd = length(out$d1)
    out$n=length(out$y)
		out$tU1y = crossprod(out$U1, out$y)
		if (out$nd < out$n) {
			out$tyy = sum(out$y^2)
		}
	}
	#must be outside of if!(is.null(Z)), otherwise, X and Zt will not be updated i nto FaST
	updateFaST.X(out, out$X)
	updateFaST.Zw(out, out$listZw)
	updateFaST.Zt(out, out$Zt)
	out$formula=formula
	class(out) = c("FaST")
	return(out)
}

#currently only add more terms to Fast
# should be able to replace and subtract in the future
##do not update Zt.

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
		FaST[["listZw"]] = NULL
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
				if (!is.null(FaST$Zt)) FaST[["tWZt"]] = crossprod(Zw, FaST$Zt)
			}
		}))
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

	eval.parent(substitute({

		if (!is.null(FaST$U1)) {
			FaST[["tU1X"]] = crossprod(FaST$U1, X)

			if (FaST$nd < FaST$n) {
				FaST[["tXX"]] = crossprod(X)
				FaST[["tXy"]] = crossprod(X, FaST$y)
				if (!is.null(FaST$Zw)) FaST[["tXW"]] = crossprod(X, FaST$Zw)
				if (!is.null(FaST$Zt)) FaST[["tXZt"]] = crossprod(X, FaST$Zt)
			}
		}
	}))
}

updateFaST.Zt = function(FaST, Zt) {
	eval.parent(substitute({
		FaST[["Zt"]] = NULL
		FaST[["tU1Zt"]] = NULL
		FaST[["tyZt"]] = NULL
		FaST[["tXZt"]] = NULL
		FaST[["tZtZt"]] = NULL
		FaST[["tWZt"]] = NULL
	}))
	if (!is.null(Zt)) {
		
    if(nrow(Zt)==FaST$n0){if(length(FaST$whichNa)>0)Zt=Zt[-FaST$whichNa,]}else{
      if (nrow(Zt) != FaST$n) {
        stop("Zt must have the same number of rows as the original y or the current y")
      }
    }
		varXt = apply(Zt, 2, function(a) var(a, na.rm = T))
		whichVar = which(varXt > 0)
		if (length(whichVar) > 0) {
			Zt = Zt[, whichVar, drop = F]
		} else {
			warning("No variants in Zt")
			Zt = NULL
		}
		eval.parent(substitute({
			FaST[["Zt"]] = Zt
			if (!is.null(Zt)) {
				if (!is.null(FaST$U1)) {
					FaST$tU1Zt = crossprod(FaST$U1, Zt)
					if (FaST$nd < FaST$n) {
						FaST[["tyZt"]] = crossprod(FaST$y, Zt)
						FaST[["tXZt"]] = crossprod(FaST$X, Zt)
						FaST[["tZtZt"]] = crossprod(Zt)
						if (!is.null(FaST$Zw)) {
							FaST[["tWZt"]] = crossprod(FaST$Zw, Zt)
						}
					}
				}
			}
		}))
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








