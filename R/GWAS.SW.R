
#GWAS on a single window.
GWAS.SW = function(y, X0 = NULL, Xt, G = NULL, W = NULL, methods = c("Score", "SKAT"), removeXtFromG = F, tU1X=NULL, tU1y=NULL, tXX=NULL, tXy=NULL, tyy=NULL) {

	if (is.null(X0)) {
		X0 = matrix(rep(1, length(y)))
	}
	X0 = as.matrix(X0)
	Xt = as.matrix(Xt)
	whNA = which(is.na(y))
	if (length(whNA) > 0) {
		y = y[-whNA]
		X0 = X0[-whNA, , drop = F]
	}
	#remove non-variants
	varXt=apply(Xt,2,var)
	whichVar=which(varXt > 0)
	if(length(whichVar)>0){
		Xt=Xt[,whichVar,drop=F]
	}else{
		warning("No variants in Xt")
		return(NA)
	}
	if (!is.null(G)) {
		if ("matrix" %in% class(G)) {
			cat("To save computation time, it is better to supply eigenG instead of G \n")
			eigenG = getEigenG(G)
		} else if (("eigenG" %in% class(G))) {
			eigenG = G
		} else {
			stop("G must be a matrix or eigenG")
		}
		if (length(whNA) > 0) {
			eigenG$U1 = eigenG$U1[-whNA, ]
		}

		U1 = eigenG$U1
		d1 = eigenG$d1
		nd=length(d1)
		n=length(y)
		if (is.null(tU1X)) 
			tU1X = crossprod(U1, X)
		if (is.null(tU1y)) 
			tU1y = crossprod(U1, y)
		if (nd < n) {
			if (is.null(tXy)) 
				tXy = crossprod(X, y)
			if (is.null(tXX)) 
				tXX = crossprod(X, X)
			if (is.null(tyy)) 
				tyy = sum(y^2)
		}
		rm(list = c("G", "U1", "d1"))
	} else {
		eigenG = NULL
	}



	if (("Score" %in% methods) | ("SKAT" %in% methods)) {
		windowtest = intersect(c("Score", "SKAT"),methods)
		if (removeXtFromG == T) {
			W = c(W, list(Zt))
			nw = length(W)
			out = testZ(y, X=X0, eigenZd = eigenG, Zt = Xt, W = W, tauRel = paste("tauw", nw, "=-taud", sep = ""), windowtest = windowtest, 
				tU1X=tU1X,tU1y=tU1y,tXX=tXX,tXy=tXy,tyy=tyy)
		} else {
			out = testZ(y, X=X0, eigenZd = eigenG, Zt = Xt, W = W, tauRel = NULL, windowtest = windowtest,tU1X=tU1X,tU1y=tU1y,tXX=tXX,tXy=tXy,tyy=tyy)
		}
		
	}

	
}
