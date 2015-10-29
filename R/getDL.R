
get_tau = function(Var, logVar, tauRel) {

	if (logVar == T) {
		Var = exp(Var)
	}

	for (i in 1:length(Var)) {
		assign(names(Var)[i], Var[i])
	}
	namesPar = names(Var)
	if (any(!gsub("\\d*", "", namesPar) %in% c("var_e", "taud", "tauw"))) {
		stop("Var names must be var_e, taud, or tauwD")
	}

	if (is.null(tauRel)) {
		names.tauw = grep("tauw", namesPar, value = T)
	} else {
		for (i in 1:length(tauRel)) {
			eval(parse(text = tauRel[[i]]))

		}
		split.tau = strsplit(tauRel, split = "=")
		names.tauw = grep("tauw", c(unlist(split.tau), namesPar), value = T)
		names.tauw = gsub(".*(tauw\\d+).*", "\\1", names.tauw)
		names.tauw = unique(names.tauw)
	}
	if (length(names.tauw) > 0) {
		index.tauw = as.numeric(gsub("tauw(\\d+)", "\\1", names.tauw))
		names.tauw = names.tauw[order(index.tauw)]

		tauw = vector()
		for (i in 1:length(names.tauw)) {
			tauw[i] = get(names.tauw[i])
		}
		names(tauw) = names(tauw)
	} else tauw = NULL
	Var = list(var_e = var_e, taud = taud, tauw = tauw)
	return(Var)
}
getDL.XYZ = function(var_e, taud, tauw = NULL, eigenG, X, y, Zw = NULL, kw = NULL, Zt = NULL) {
	n = length(y)
	U1 = eigenG$U1
	d1 = eigenG$d1
	tXX = crossprod(X)
	tU1y = crossprod(U1, y)
	tU1X = crossprod(U1, X)
	tXy = crossprod(X, y)
	tyy = sum(y^2)
	if (!is.null(Zw)) {
		tU1W = crossprod(U1, Zw)
		tXW = crossprod(X, Zw)
		tWW = crossprod(Zw)
		tWy = crossprod(Zw, y)
	} else {
		tU1W = tXW = tWW = tWy = NULL

	}
	if (!is.null(Zt)) {
		tZtZt = crossprod(Zt)
		tU1Zt = crossprod(U1, Zt)
		tXZt = crossprod(X, Zt)
		tyZt = crossprod(y, Zt)
		if (!is.null(Zw)) {
			tWZt = crossprod(Zw, Zt)
		} else {
			tWZt = NULL
		}
	} else {
		tZtZt = tU1Zt = tXZt = tyZt = NULL
	}
	out = getDL(var_e, taud, d1 = d1, n = n, tU1y = tU1y, tU1X = tU1X, tXX = tXX, tXy = tXy, tyy = tyy, tauw = tauw, kw = kw, tU1W = tU1W, tXW = tXW, 
		tWW = tWW, tWy = tWy, tZtZt = tZtZt, tU1Zt = tU1Zt, tXZt = tXZt, tyZt = tyZt, tWZt = tWZt, getQ = F, getS = F, getNeg2Log = T, REML = T)

	return(out)
}



getDL = function(var_e, taud, d1, n, tU1y, tU1X, tXX, tXy, tyy, tauw = NULL, kw = NULL, tU1W = NULL, tXW = NULL, tWW = NULL, tWy = NULL, tZtZt = NULL, 
	tU1Zt = NULL, tXZt = NULL, tyZt = NULL, tWZt = NULL, getQ = F, getS = F, getNeg2Log = T, REML = T) {
	if (is.null(tauw)) 
		tauw = NA
	if (is.null(kw)) 
		kw = NA
	if (is.null(tU1W)) 
		tU1W = NA
	if (is.null(tXW)) 
		tXW = NA
	if (is.null(tWW)) 
		tWW = NA
	if (is.null(tWy)) 
		tWy = NA
	if (is.null(tZtZt)) 
		tZtZt = NA
	if (is.null(tU1Zt)) 
		tU1Zt = NA
	if (is.null(tXZt)) 
		tXZt = NA
	if (is.null(tyZt)) 
		tyZt = NA
	if (is.null(tWZt)) 
		tWZt = NA
	getQ = as.integer(getQ)
	getS = as.integer(getS)
	getNeg2Log = as.integer(getNeg2Log)
	REML = as.integer(REML)
	out <- .Call("C_getDL", var_e, taud, d1, n, tU1y, tU1X, tXX, tXy, tyy, tauw, kw, tU1W, tXW, tWW, tWy, tZtZt, tU1Zt, tXZt, tyZt, tWZt, getQ, 
		getS, getNeg2Log, REML)
	return(out)

}






