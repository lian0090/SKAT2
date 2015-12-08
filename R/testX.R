#t-test on individual markers
#It is not exactly the same as that from emma, due to the small difference in estimating variance components. 
#If I use the same variance component from emma to put into getDL, I will get exactly the same pvalue 
#Population structure previously determined. 
##perform association mapping for provided markers while correcting for multiple test.

# GWAS.P3D = function(null.formula, Xt, multipleCorrection = T, method = "LR") {
# 	#P3D0 allows previously defined P3D0, this might be useful if you are constantly testing you code for small number of markers  
# 	
# 	data=formToMatrix(null.formula)
# 	y=data$y
# 	X0=data$X
# 	fit0=fitNULL(y,X=X0,Z=data$Z)		
# 	whichNa =FaST0$whichNa
# 	##remove non-variants
# 	Xt=checkX(Xt,FaST0$n0,FaST0$whichNa)
# 	 
# 	if (multipleCorrection == T) {
# 		Me = Meff(Xt)
# 		cat("effective number of test is ", Me, "\n")
# 	} else {
# 		Me = 1
# 	}
# 	p.value = rep(NA, ncol(Xt))
# 
# 	p.value = apply(Xt, 2, function(a) {
# 		testX(FaST0, fit0=fit0, Xt = a, method = method, P3D=T)$p.value * Me
# 	})
# 	return(list(Me = Me, p.value = p.value, H0var = Var))
# }


testX = function(fit0, Xt, methods = c("SSNP.P3D.LR","SSNP.P3D.t","SSNP.LR","SSNP.t")[1]) {
	##Var is the variance for the NULL model, without fitting the test SNPs 
	#P3D=T, NULL model and alternative model share the same variance components.	
	##if Xt is NULL return NA

	if(is.null(Xt)) return(NA)	
	if (grepl("\\.FaST",class(fit0))) {
		#for t-test when P3D is False. Otherwise, use lm or lmm for fit0.
		if(any(grepl("\\.P3D",methods)) | any(grepl("\\.LR",methods))) stop("must supply a fitted fit0 from fitNULL if P3D is true or if LR in methods")	
		Var0=NA
		ln0=NA
	}else{
		#fit0 as lmm or lm is only needed for LR or P3D. Not used for SSNP.t
	 if(class(fit0)=="lmm"){
	 	Var0 = fit0$Var
	 	ln0 = -1/2*getDL(Var0, fit0$FaST, getNeg2Log = T, REML = F)$neg2logLik
	 }
	 if(class(fit0)=="lm") {
	 	Var0=NA
	 	ln0=logLik(fit0,REML=F)
	 }
	}
	 class.fit0=class(fit0)
	 
	 
	Xt = as.matrix(Xt)
	ntest = ncol(Xt)
	test = c((ncol(fit0$FaST$X) + 1):(ncol(fit0$FaST$X) + ntest))
	

	#get Var and ln0 befrore update FaST0
	 
	 #this can be optimized (fit0 has its own ML likelihood) because we only need to get neg2logLik once. 
	 FaST1=fit0$FaST
	 rm("fit0")
	 updateFaST.X(FaST1, X = cbind(FaST1$X, Xt))
   n.beta=ncol(FaST1$X)
     
   namesout=c("LR","p.value","logML0","logML1","Var0","Var1")
	 out=vector("list",length=length(namesout))	 
	 names(out)=namesout
	 out=lapply(out,function(a){a=rep(NA,length=length(methods));names(a)=methods;return(a)})
   out$Var0=as.list(out$Var0)
   out$Var1=as.list(out$Var1)
	 LR.methods=grep("\\.LR",methods,value=T)
	 t.methods=grep("\\.t",methods,value=T)
	 P3D.methods=grep("\\.P3D",methods,value=T)
	 NP.methods=setdiff(methods,P3D.methods)
	
	if (class.fit0 %in% c("lm","lm.FaST")) {
		##P3D is irrelavent for lm models
		fit1 = fitNULL.FaST(FaST1)
			if (length(LR.methods)>0) {
				LR = -2 * (ln0 - logLik(fit1))
				p.LR = pchisq(LR, df = ntest, lower.tail = F)
				out$p.value[LR.methods]=p.LR
			}
			if (length(t.methods)>0) {
				p.t = min(summary(fit1)$coefficients[test, 4],na.rm=T)
				out$p.value[t.methods]=p.t
				}
		##end methods for NULL Z0
		} else if (class.fit0 %in% c("lmm","lmm.FaST")){
			        
  get.test.lmm<-function(out,FaST1,ln0,Var0,Var1,methods,ntest){
       outDL = getDL(Var1, FaST1, getNeg2Log = T, REML = F, getAlphaHat=T)

    for (i in 1:length(methods)) {
			if (grepl("\\.LR",methods[i])) {
				trypvalue = try({
					ln1 = -1/2 * outDL$neg2logLik
					Q = -2 * (ln0 - ln1)
					p.value = pchisq(Q, df = ntest, lower.tail = F)
					out$p.value[methods[i]]=p.value
					out$logML1[methods[i]]=ln1
					out$logML0[methods[i]] = ln0
					out$LR[methods[i]]= Q
					out$Var0[[methods[i]]] = Var0
					out$Var1[[methods[i]]] = Var1
				})
			} else if (grepl("\\.t", methods[i])) {
 				trypvalue = try({
 					beta = outDL$hat_alpha
 					vbeta = outDL$invtXVinvX
 					tscore=sapply(test,function(j){abs(beta[j])/sqrt(vbeta[j, j])})
 					##note: the df for t-distribution is not corrected by Satterthwaite's method. Likelihood ratio test should be better.	
 					p.t= 2 * pt(tscore, df = n - n.beta, lower.tail = F)
 					out$p.value[methods[i]]=min(p.t,na.rm=T)
 				})
			  
 			}
		}	
       return(out)
     }

		if(length(P3D.methods)>0) {
     	       Var1 = Var0
            out=get.test.lmm(out,FaST1=FaST1,ln0=ln0, Var0,Var1,methods=P3D.methods,ntest)
          }

     	if(length(NP.methods)>0){
     		 fit1 = fitNULL.FaST(FaST1)
     	     Var1=fit1$Var
     	     out=get.test.lmm(out,FaST1=FaST1,ln0=ln0,Var0=Var0,Var1=Var1,methods=NP.methods,ntest)
     	}
     }
        	
	
	return(out)

}


