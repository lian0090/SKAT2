#t-test on individual markers
#It is not exactly the same as that from emma, due to the small difference in estimating variance components. 
#If I use the same variance component from emma to put into getDL, I will get exactly the same pvalue 
#Population structure previously determined. 
##perform association mapping for provided markers while correcting for multiple test.


GWAS.P3D = function(y, Xt, X0 = NA, Z0 = NA, fit0 = NA, multipleCorrection = T, method = "LR") {
	#P3D0 allows previously defined P3D0, this might be useful if you are constantly testing you code for small number of markers  
	FaST0=getFaST(y,X=X0,Z=Z0)
		
	if (is.na(fit0)) {
		fit0 = fitlmm.FaST(FaST0)
	}
	whichNa =FaST0$whichNa
	Xt = as.matrix(Xt)
	##remove non-variants
	if (length(whichNa) > 0) 
		Xt = Xt[-whichNa, ]
	    varXt = apply(Xt, 2, var)
	    whichVar = which(varXt > 0)
	if (length(whichVar) > 0) {
		Xt = Xt[, whichVar, drop = F]
	} else {
		warning("No variants in Xt")
		return(NA)
	}

	if (multipleCorrection == T) {
		Me = Meff(Xt)
		cat("effective number of test is ", Me, "\n")
	} else {
		Me = 1
	}
	p.value = rep(NA, ncol(Xt))

	p.value = apply(Xt, 2, function(a) {
		singleXt.FaST(FaST0, fit0=fit0, Xt = a, method = method, P3D=T)$p.value * Me
	})
	return(list(Me = Me, p.value = p.value, H0var = Var))
}


testX = function(FaST0, fit0 = NA, Xt, method = "LR",P3D=T) {
	##Var is the variance for the NULL model, without fitting the test SNPs 
	#P3D=T, NULL model and alternative model share the same variance components.	
out$p.value = rep(0, length(method))
	names(out$p.value) = method

	out = list()
	Xt = as.matrix(Xt)
	ntest = ncol(Xt)
	test = c((ncol(FaST0$X) + 1):(ncol(FaST0$X) + ntest))

	if (ntest > 1) {
		if ("t" %in% method) {
			stop("use LR test when there is more than one fix effect to be tested")
		}
	}

	##tests for NULL Z0
	if ("LR" %in% method) {
		if (is.na(fit0)) {
			cat("it is better to supply fit0 when using LR test ")
			fit0 = fitlmm.FaST(FaST0)
		}
	}
	FaST1 = updateFaST.X(FaST0, X = cbind(X0, Xt))


	if (class(fit0) == "lm") {
		fit1 = fitlmm.FaST(updateFaST.X(FaST0, X = cbind(X0, Xt)))

		for (i in 1:length(method)) {
			if (method[i] == "LR") {
				Q = -2 * (logLik(fit0) - logLik(fit1))
				out$p.value[i] = pchisq(Q, df = ntest, lower.tail = F)
			}
			if (method[i] == "t") {
				out$p.value[i] = summary(fit1)$coefficients[n.beta, 4]
			}

		}
		##end methods for NULL Z0
		} else {

		Var = fit0$Var

		if (P3D == T) {
			Var1 = Var
		} else {
			fit1 = fitlmm.FaST(FaST1)
			Var1 = fit1$Var
		}
		outDL = getDL(Var1, FaST1, getNeg2Log = T, REML = F)
		for (i in 1:length(method)) {
			if ("LR" == method[i]) {
				trypvalue = try({
					ln0 = fit0$loglik
					ln1 = -1/2 * outDL$neg2logLik
					Q = -2 * (ln0 - ln1)
					p.value = pchisq(Q, df = ntest, lower.tail = F)
					out$ML1 = ln1
					out$ML0 = ln0
					out$LR = Q
					out$Var = Var
					out$Var1 = Var1
				})
			} else if ("t" == method[i]) {
				trypvalue = try({
					beta = outDL$hat_alpha
					vbeta = outDL$invtXVinvX
					tscore = abs(beta[test])/sqrt(vbeta[test, test])
					##note: the df for t-distribution is not corrected by Satterthwaite's method. Likelihood ratio test should be better.	
					p.value = 2 * pt(tscore, df = n - n.beta, lower.tail = F)
				})
			}
			if (inherits(trypvalue, "try-error")) {
				p.value = NA
			}
			out$p.value[i] = p.value
		}
	}
	return(out)
}
