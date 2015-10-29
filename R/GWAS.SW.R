
#GWAS on a single window.
GWAS.SW = function(y, X0 = NULL, Z0 = list(), Zt, methods = c("Score", "SKAT"), removeZtFromZ0 = F, tU1X = NULL, tU1y = NULL, tXX = NULL, tXy = NULL, tyy = NULL) {

	if (is.null(X0)) {
		X0 = matrix(rep(1, length(y)))
	}
	X0 = as.matrix(X0)
	Zt = as.matrix(Zt)
	whNA = which(is.na(y))
	if (length(whNA) > 0) {
		y = y[-whNA]
		X0 = X0[-whNA, , drop = F]
	}
	#remove non-variants
	varXt = apply(Zt, 2, var)
	whichVar = which(varXt > 0)
	if (length(whichVar) > 0) {
		Zt = Zt[, whichVar, drop = F]
	} else {
		warning("No variants in Zt")
		return(NA)
	}
	if (length(Z0) < 1) {
		eigenG = NULL
	} else {
		if (!("eigenG" %in% class(Z0[[1]]))) {
			cat("It is better to supply the first element of Z as an eigenG object to save computation\n")
			eigenG = getEigenG(Zg = Z0[[1]])
		} else {
			eigenG = Z[[1]]
		}
		if (length(Z0) == 1) {
			Zw = NULL
		} else {
			Zw = Z0[-1]
		}
		if (length(whNA) > 0) {
			eigenG$U1 = eigenG$U1[-whNA, ]
		}

		U1 = eigenG$U1
		d1 = eigenG$d1
		nd = length(d1)
		n = length(y)
		if (is.null(tU1X)) 
			tU1X = crossprod(U1, X0)
		if (is.null(tU1y)) 
			tU1y = crossprod(U1, y)
		if (nd < n) {
			if (is.null(tXy)) 
				tXy = crossprod(X0, y)
			if (is.null(tXX)) 
				tXX = crossprod(X0, X0)
			if (is.null(tyy)) 
				tyy = sum(y^2)
		}
		rm(list = c("U1", "d1"))
	}



	if (("Score" %in% methods) | ("SKAT" %in% methods)) {
		windowtest = intersect(c("Score", "SKAT"), methods)
		if (removeZtFromZ0 == T) {
			Zw = c(Zw, list(Zt))
			nw = length(Zw)
			out = testZ(y, X = X0, eigenG = eigenG, Zt = Zt, Zw = Zw, tauRel = paste("tauw", nw, "=-taud", sep = ""), windowtest = windowtest, tU1X = tU1X, tU1y = tU1y, 
				tXX = tXX, tXy = tXy, tyy = tyy)
		} else {
			out = testZ(y, X = X0, eigenG = eigenG, Zt = Zt, Zw = Zw, tauRel = NULL, windowtest = windowtest, tU1X = tU1X, tU1y = tU1y, tXX = tXX, tXy = tXy, tyy = tyy)
		}

	}


}


testZ = function(y, X, Zw = NULL, tauRel = NULL, Zt, eigenG, windowtest, tU1X = NULL, tU1y = NULL, tXX = NULL, tXy = NULL, tyy = NULL, logVar = T, optimizer = "bobyqa") {

	#check input
	X = as.matrix(X)
	Zt = as.matrix(Zt)
	if (any(is.na(y))) {
		#optim function will report not being able to evalue function at intial values when there is NA
		if ((!is.null(tU1X)) | (!is.null(tU1y)) | !is.null(tXX) | (!is.null(tXy)) | (!is.null(tyy))) {
			stop("there should be no missing values when tU1X, tU1y, tXX, tXy or tyy is supplied")
		}
		whNAy = which(is.na(y))
		y = y[-whNAy]
		X = X[-whNAy, , drop = F]
		Zt = Zt[-whNAy, , drop = F]
		eigenG$U1 = eigenG$U1[-whNAy, ]
	} else {
		whNAy = NULL
	}
	out = list()
	#Null model with no random effects
	#classical SKAT test
if (!is.null(windowtest)) {
		if (is.null(eigenG)) {
			mod = lm(y ~ -1 + X)
			resid = residuals(mod)
			s2 = summary(mod)$sigma^2
			Q = sum(crossprod(resid, Zt)^2)/s2/2
			tXZt = crossprod(X, Zt)
			Zw.1 = (crossprod(Zt) - t(tXZt) %*% solve(crossprod(X)) %*% tXZt)/2
			lambda = eigen(Zw.1, symmetric = TRUE, only.values = TRUE)$values
			lambda1 = lambda
			IDX1 <- which(lambda >= 0)
			# eigenvalue bigger than mean(lambda1[IDX1])/100000 
			IDX2 <- which(lambda1 > mean(lambda1[IDX1])/1e+05)
			lambda <- lambda1[IDX2]
			varS = (sum(Zw.1^2)/(s2^2))/2
			S = (Q - (sum(diag(Zw.1))))/s2
			out$p.SKAT <- Get_PValue.Lambda(lambda, Q)
			out$Q = Q
			out$S = S
			out$p.Score = pnorm(S, mean = 0, sd = sqrt(varS), lower.tail = F)
			return(out)
		}
	}

	#logVar, paramterize variance components with log when using REML to restrict variance component to be larger than 0.
	
	d1 = eigenG$d1
	U1 = eigenG$U1

	if (is.null(tU1X)) 
		tU1X = crossprod(U1, X)
	if (is.null(tU1y)) 
		tU1y = crossprod(U1, y)

	n = length(y)
	nd = length(d1)


	if (nd < n) {
		if (is.null(tXy)) 
			tXy = crossprod(X, y)
		if (is.null(tXX)) 
			tXX = crossprod(X, X)
		if (is.null(tyy)) 
			tyy = sum(y^2)
	}

	if (!is.null(Zw)) {

		if (!is.list(Zw)) {
			stop("Zw must be a list of all other random effect incidence matrix")
		}
		kw = sapply(Zw, ncol)
		nw = length(kw)
		Zw = do.call(cbind, Zw)
		if (sum(kw) != ncol(Zw)) 
			stop("sum of kw should be equal to the number of columns in Zw")
		#remove NA values
		if (length(whNAy) > 0) {
			Zw = Zw[-whNAy, , drop = F]
		}
		tauw = rep(0, nw)
		tU1W = crossprod(U1, Zw)
		tXW = crossprod(X, Zw)
		tWW = crossprod(Zw, Zw)
		tWy = crossprod(Zw, y)
		if (!is.null(windowtest)) {
			tWZt = crossprod(Zw, Zt)
		} else {
			tWZt = NULL
		}
	} else {
		nw = 0
		kw = NULL
		tauw = NULL
		tU1W = NULL
		tXW = NULL
		tWW = NULL
		tWy = NULL
		tWZt = NULL
	}

	if (!is.null(windowtest)) {
		tU1Zt = crossprod(U1, Zt)
		tyZt = crossprod(y, Zt)
		tXZt = crossprod(X, Zt)
		tZtZt = crossprod(Zt)
	}

	##test with low rank Zh	
	namesPar = c("var_e", "taud")
	if (nw > 0) {
		namesPar = c("var_e", "taud", paste("tauw", c(1:nw), sep = ""))
		if (!is.null(tauRel)) {
			splittau = strsplit(tauRel, split = "=")
			unlist.splittau = unlist(splittau)
			#only the independent Var need to be positive if using logVar
			#exp transformation will only be applied to par, not the dependent var
#if(any(grepl("-[[:punct:][:digit:]]*tau",unlist.splittau))){if(logVar==T)stop("logVar must be set to False when the tau values are of different signs")}
n.Rel = length(splittau)
			if (length(unlist.splittau) > 2 * length(splittau)) 
				stop("Every equation must contain only one relationship")
			#strip off the digits, decimal point (.) and +,-,*,/
			if (any(!(gsub("([[:digit:][:punct:]e]+)", "", unlist.splittau) %in% c("taud", "tauw", "")))) {
				stop("in tauRel, must only specify taud or tauwD, where D is any integer less or equal to the number of Zw matrix")
			}
			left.tau = gsub(".*(tau[wd]{1}\\d*).*", "\\1", sapply(splittau, function(a) a[1]))
			right.tau = gsub(".*(tau[wd]{1}\\d*).*", "\\1", sapply(splittau, function(a) a[2]))
			if (any(duplicated(left.tau))) {
				stop("duplicated terms on the lest side of equations")
			}
			if (length(intersect(right.tau, left.tau)) > 0) {
				stop("one variable can only appear on one side of equation")
			}

			for (i in 1:length(splittau)) {
				splittau.i = gsub(".*(tau[dw]{1}\\d*).*", "\\1", splittau[[i]])

				n.split = length(splittau.i)
				##only keep the last tau from relationship equation
				namesPar = setdiff(namesPar, splittau.i[-n.split])
			}
		}
	}
	par = rep(0.5, length(namesPar))
	names(par) = namesPar

	fit0 = fit.optim(par = par, fn = neg2Log, logVar = logVar, d1 = d1, n = n, tU1y = tU1y, tU1X = tU1X, tXX = tXX, tXy = tXy, tyy = tyy, tU1W = tU1W, tXW = tXW, tWW = tWW, 
		tWy = tWy, kw = kw, tauRel = tauRel, optimizer = optimizer)
	out$fit0 = fit0

	##SKAT test or LR test
	
	if ("SKAT" %in% windowtest | "Score" %in% windowtest) {
		var_e = fit0$outVar$var_e
		taud = fit0$outVar$taud
		if (nw > 0) 
			tauw = fit0$outVar$tauw
		getQ = ("SKAT" %in% windowtest)
		getS = ("Score" %in% windowtest)
		Qdis = getDL(var_e, taud = taud, d1 = d1, n = n, tU1y = tU1y, tU1X = tU1X, tXX = tXX, tXy = tXy, tyy = tyy, tauw = tauw, kw = kw, tU1W = tU1W, tXW = tXW, tWW = tWW, 
			tWy = tWy, tZtZt = tZtZt, tU1Zt = tU1Zt, tXZt = tXZt, tyZt = tyZt, tWZt = tWZt, getQ = getQ, getS = getS, getNeg2Log = F, REML = T)

		if ("SKAT" %in% windowtest) {
			Q = Qdis$Q
			lambda = Qdis$lambda
			lambda1 = lambda
			IDX1 <- which(lambda >= 0)
			#eigenvalue bigger than mean(lambda1[IDX1])/100000 
			IDX2 <- which(lambda1 > mean(lambda1[IDX1])/1e+05)
			lambda <- lambda1[IDX2]
			p.value <- Get_PValue.Lambda(lambda, Q)
			out$p.SKAT = p.value
			out$Q = Q
		}
		#Score test
		if ("Score" %in% windowtest) {
			S = Qdis$S
			sdS = Qdis$sdS
			out$Score = S
			out$sdScore = sdS
			out$p.Score = pnorm(S/sdS, lower.tail = F)
		}

	}
	return(out)
}

