##testing
#library(BGLR)
# data(mice)
# fast1=getFaST(y=mice.pheno$Obesity.BMI,Z=list(Z1=mice.X[,1:20]))
# updateFaST.Zt(fast1,Zt=mice.X[,20:22])
# updateFaST.Zw(fast1,list(mice.X[,1:2]))

getFaST = function(y, X = NULL, Z = NULL, Zt = NULL) {
	out = list()
	out$n0 = length(y)
	whichNa = which(is.na(y))
	if (length(whichNa) > 0) {
		y = y[-whichNa]
	}
	out$y = y
	out$whichNa = whichNa
	out$n = length(y)
	listZw = NULL
	if (!is.null(Z)) {
		if (!("eigenG" %in% class(Z[[1]]))) {
			cat("It is better to supply the first element of Z as an eigenG object to save computation\n")
			eigenG = getEigenG(Zg = Z[[1]])
		} else {
			eigenG = Z[[1]]
		}
		if (length(whichNa) > 0) {
			eigenG$U1 = eigenG$U1[-whichNa, ]
		}
		if (length(Z) > 1) {
			listZw = Z[-1]
		}
		out$U1 = eigenG$U1
		out$d1 = eigenG$d1
		out$nd = length(out$d1)
		out$tU1y = crossprod(out$U1, out$y)
		if (out$nd < out$n) {
			out$tyy = sum(out$y^2)
		}
	}	
	#must be outside of if!(is.null(Z)), otherwise, X and Zt will not be updated i nto FaST
	updateFaST.X(out, X)
	updateFaST.Zw(out, listZw)
	updateFaST.Zt(out, Zt)
	class(out) = c("FaST")
	return(out)
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
		if (nrow(listZw[[1]]) != FaST$n0) {
			stop("Zw must have the same number of rows as the original y")
		}
		if (length(FaST$whichNa) > 0) {
			listZw = lapply(listZw, function(a) {
				a[-whichNa, , drop = F]
			})
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
		if (nrow(X) != FaST$n0) {
			stop("X must have the same number of rows as the original y")
		}
		if (length(FaST$whichNa) > 0) {
			X = as.matrix(X[-whichNa, , drop = F])
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
		if (nrow(Zt) != FaST$n0) {
			stop("Zt must have the same number of rows as the original y")
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



