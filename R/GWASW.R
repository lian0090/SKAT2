
# # #GWAS on a marker windows .
# GWASW= function(y, MxE.formula=~(1|M)+(1|E)+(1|M:E), X = NULL, Z = NULL,  sets = rep(c(1:ceiling(ncol(M)/10)), each = 10)[1:ncol(M)], methods = c("Score", "SKAT")) {
 # fit0 = getFaST(y, X = X, Z = Z, Zt = NULL)   
    # nSets = length(unique(sets))
    # nmethods = length(methods)
    # out = matrix(nrow = nSets, ncol = nmethods)
    	# bars=findbars(MxE.formula)
    	# sapply(bar.terms,function(a){bar.terms[[1]]%in%a})
    	# nobars=findnobars(MxE.formula)
    	# subbar.form=subbars(MxE.formula)
    	# subbar.terms=lapply(bars,FUN=all.vars)

    	# termlabels=attr(terms(subbar.form),"term.label")
    	    # nterms=length(termlabels)
    	    # if(nterms==0){
    	    	# stop("must supply test terms")
    	    # }
  # if(nterms>=1){
    	    	# if(is.null(bars)){
    	    		# stop("Must supply at least one random effect")
    	    	# }

    	# if(nterms==1){
    # for (i in 1:nSets) {
    	         	    	# ##extract sets from the first term
    	    	 # substitute(Mi<-M[,which(sets==i)],list(M=quote(termlabels[1])))
    	    	
    	   # updateFaST.Zt(fit0$FaST, Zt = Zti)
    	    	 # if(nterms>1){
    	    	 	
    	    	 	
    	    	 # }
        # }
        # if(nterms)
    	    	
    	    
    	# barsi[[1]][3]=quote(Mi)
    	# bars
    	# MxEi.formula=substitute(MxE.formula,)
        # MxE.formula[[3]]
        # Mi = M[, which(sets == i)]
        # updateFaST.Zt(fit0$FaST, Zt = Zti)
        # if()
        # out[i, ] = testZ(fit0, methods = methods)[paste("p.", methods, sep = "")]
    # }
    # return(out)
# }

# GWASW=function(,){
	
	
	
	
	
# }



# FormToMatrix=function(formula){

# #window will always be divided based on the first term
# #test will always be done on the last term.
# denv=parent.env(environment())
# formula=as.formula(formula,env=denv)
# fr.form <- subbars(formula) # substitute "|" by "+" 
# environment(fr.form) <- environment(formula)

# ##get a whole model.frame before you get X, y, Z to avoid environment contamination.
# fr=model.frame(fr.form)

# ##get fixed term
# fixedform <- MxE.formula
# nb <- nobars(RHSForm(fixedform))

# if(is.null(nb))X=NULL else {
  # RHSForm(fixedform)<-nb
  # X=model.matrix(fixedform,fr)
# }

# #get random term
# bars=findbars(RHSForm(formula))
# Z=lapply(bars,FUN=function(x) model.matrix(as.formula(substitute(~-1+foo,list(foo=x[[3]]))),fr))
# names(Z) <- unlist(lapply(bars, function(x) deparse(x[[3]])))

# y <- model.response(fr)
# return(list(Z=Z,X=X,y=y))

# }

# #an example for GxE analysis 
# # GWAS.GxE.MW = function(y, Xf, Xe, Z0, M, sets = rep(c(1:ceiling(ncol(M)/10)), each = 10)[1:ncol(M)],  methods = c("Score", "SKAT"),removeMsetFromZ0 = F) {
 # # if (is.na(Z0)) {
        # # if (removeZtFromZ0 == T) {
            # # nw = length(Z0) - 1
            # # Z = c(Z0, list(Zti))
            # # tauRel = paste("tauw", nw, "=-taud", sep = "")
        # # } else {
            # # Z = Z0
            # # tauRel = NA
        # # }
    # # } else {
        # # Z = NA
        # # tauRel = NA
    # # }


# # }

# testXZ=function(fit0,XZ=list(X=,Z=)){
	
	
# }


testZ = function(fit0, Zt, methods = "Score") {

	out = vector("list",length=5)
	names(out)=c("Score","sdScore","p.Score","p.SKAT","Q")
	
    
    	updateFaST.Zt(fit0$FaST, Zt)
    if(is.null(fit0$FaST$Zt)){
    	return(out)
    }

	if (class(fit0)=="lm") {
		out=with(fit0$FaST, {
			out = vector("list",length=5)
        	names(out)=c("Score","sdScore","p.Score","p.SKAT","Q")
			resid = residuals(fit0)
			s2 = summary(fit0)$sigma^2
			Q = sum(crossprod(resid, Zt)^2)/s2/2
			tXZt = crossprod(X, Zt)
			Zw.1 = (crossprod(Zt) - t(tXZt) %*% solve(crossprod(X)) %*% tXZt)/2
			lambda = eigen(Zw.1, symmetric = TRUE, only.values = TRUE)$values
			lambda1 = lambda
			IDX1 <- which(lambda >= 0)
			# eigenvalue bigger than mean(lambda1[IDX1])/100000 
			IDX2 <- which(lambda1 > mean(lambda1[IDX1])/1e+05)
			lambda <- lambda1[IDX2]

			if ("SKAT" %in% methods) {
				out$p.SKAT <- Get_PValue.Lambda(lambda, Q)$p.value
				out$Q = Q
			}
			if ("Score" %in% methods) {
				varS = (sum(Zw.1^2)/(s2^2))/2
				S = (Q - (sum(diag(Zw.1))))/s2
				out$Score = S
				out$sdScore = sqrt(varS)
				out$p.Score = pnorm(S, mean = 0, sd = sqrt(varS), lower.tail = F)
			}
			return(out)
		}
		)
		return(out)
		
	} else {
		#logVar, paramterize variance components with log when using REML to restrict variance component to be larger than 0.
		##SKAT test or LR test
		getQ = ("SKAT" %in% methods)
		getS = ("Score" %in% methods)
		Qdis = getDL(fit0$Var, fit0$FaST, getQ = getQ, getS = getS, getNeg2Log = F, REML = T)

		if ("SKAT" %in% methods) {
			Q = Qdis$Q
			lambda = Qdis$lambda
			lambda1 = lambda
			IDX1 <- which(lambda >= 0)
			#eigenvalue bigger than mean(lambda1[IDX1])/100000 
			IDX2 <- which(lambda1 > mean(lambda1[IDX1])/1e+05)
			lambda <- lambda1[IDX2]
			p.value <- Get_PValue.Lambda(lambda, Q)$p.value
			out$p.SKAT = p.value
			out$Q = Q
		}
		#Score test
		if ("Score" %in% methods) {
			S = Qdis$S
			sdS = Qdis$sdS
			out$Score = S
			out$sdScore = sdS
			out$p.Score = pnorm(S/sdS, lower.tail = F)
		}
	}


	return(out)
}

