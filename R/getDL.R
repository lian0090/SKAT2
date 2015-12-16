#Z=sweep(op(X),v,MARGINopX,"*")
sweep_prod=function(X,v,transX,MARGINopx){
 return(.Call("Rsweep_prod",X,v,as.integer(transX),MARGINopx) )
}

##reduce the initial values of var
reduceInitVar = function(nw, VarRel) {
		namesPar = c("VarE", "VarG")
      if (nw > 0) {
	namesPar = c("VarE", "VarG", paste("VarW", c(1:nw), sep = ""))
		}

	##return the independent var. for example, VarW1=VarW2, then VarW1 will be removed from the output
	if (!is.null(VarRel)) {
		splitvar = strsplit(VarRel, split = "=")
		unlist.splitvar = unlist(splitvar)
		#only the independent Var need to be positive if using logVar
		#exp transformation will only be applied to par, not the dependent var

		n.Rel = length(splitvar)
		if (length(unlist.splitvar) > 2 * length(splitvar)) 
			stop("Every equation must contain only one relationship")
		#strip off the digits, decimal point (.) and +,-,*,/
		if (any(!(gsub("([[:digit:][:punct:]e]+)", "", unlist.splitvar) %in% c("VarG", "VarW", "")))) {
			stop("in VarRel, must only specify VarG or VarWD, where D is any integer less or equal to the number of Zw matrix")
		}
		left.var = gsub(".*(Var[wd]{1}\\d*).*", "\\1", sapply(splitvar, function(a) a[1]))
		right.var = gsub(".*(Var[wd]{1}\\d*).*", "\\1", sapply(splitvar, function(a) a[2]))
		if (any(duplicated(left.var))) {
			stop("duplicated terms on the lest side of equations")
		}
		if (length(intersect(right.var, left.var)) > 0) {
			stop("one variable can only appear on one side of equation")
		}

		for (i in 1:length(splitvar)) {
			splitvar.i = gsub(".*(Var[dw]{1}\\d*).*", "\\1", splitvar[[i]])

			n.split = length(splitvar.i)
			##only keep the last var from relationship equation
			namesPar = setdiff(namesPar, splitvar.i[-n.split])
		}
	}
     VarReduced = rep(0.5, length(namesPar))
     names(VarReduced) = namesPar

	return(VarReduced)
}

##expand the independent var based on VarRel. For example VarW1=-VarW2, VarW2 will be input as Var, and the value of VarW1 will be filled based on VarW1=-VarW2
expandVar = function(Var, logVar, VarRel) {
	#if input Var is in log scale, this function will put it back to the original scale. 
	if (logVar == T) {
		Var = exp(Var)
	}


	if (!is.null(VarRel)) {
		for (i in 1:length(Var)) {
			assign(names(Var)[i], Var[i])
		}
		namesPar = names(Var)
		if (any(!gsub("\\d*", "", namesPar) %in% c("VarE", "VarG", "VarW"))) {
			stop("Var names must be VarE, VarG, or VarWD")
		}

		for (i in 1:length(VarRel)) {
			eval(parse(text = VarRel[[i]]))
		}
		split.var = strsplit(VarRel, split = "=")
		names.VarW = grep("VarW", c(unlist(split.var), namesPar), value = T)
		names.VarW = gsub(".*(VarW\\d+).*", "\\1", names.VarW)
		names.VarW = unique(names.VarW)

		if (length(names.VarW) > 0) {
			index.VarW = as.numeric(gsub("VarW(\\d+)", "\\1", names.VarW))
			names.VarW = names.VarW[order(index.VarW)]

			VarW = vector()
			for (i in 1:length(names.VarW)) {
				VarW[i] = get(names.VarW[i])
			}
			Var = c(VarE, VarG, VarW)
			names(Var) = c("VarE", "VarG", names.VarW)
		} else {
			Var = c(VarE, VarG)
			names(Var) = c("VarE", "VarG")
		}
	}

	return(Var)
}

##does not depend on the real names of Var. 
getDL = function(Var, FaST, getQ = F, getS = F, getNeg2Log = T, REML = T, getAlphaHat=F) {
	if (FaST$nw == 0) {
		VarW = NULL
	} else {
		VarW = Var[-c(1:2)]
	}
	##.Call will not be able to take NULL values. 
	FaSTnames=c("d1","n","tU1y","tU1X","tXX","tXy","tyy","kw","nw","tU1W","tXW","tWW","tWy","tZtZt","tU1Zt","tXZt","tyZt","tWZt")
	for(i in 1:length(FaSTnames)){
		namei=FaSTnames[i]
		if(is.null(FaST[[namei]])){
			FaST[[namei]]=NA
		}
	}
	
	out <- .Call("C_getDL", Var[1], Var[2], FaST$d1, FaST$n, FaST$tU1y, FaST$tU1X, FaST$tXX, FaST$tXy, FaST$tyy, VarW, FaST$kw, FaST$nw, FaST$tU1W, FaST$tXW, FaST$tWW, FaST$tWy, FaST$tZtZt, FaST$tU1Zt, FaST$tXZt, FaST$tyZt, 
		FaST$tWZt, as.integer(getQ), as.integer(getS), as.integer(getNeg2Log), as.integer(REML),as.integer(getAlphaHat))
	return(out)

}






