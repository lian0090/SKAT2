#  if (!is.null(Z)) {
#	classZ=sapply(Z,class)
#whichNumeric=which(classZ=="numeric")
#if(length(whichNumeric)>0){
#	for(i in whichNumeric){
#		Z[[i]]=as.matrix(Z[[i]])
#	}
#}
#whichEigenG=which(classZ=="eigenG")
#if(length(whichEigenG)>1){
#	stop("only one eigenG object is allowed in Z")
#}
#if (length(whichEigenG)==0) {
#	cat("It is better to supply the one element of Z as an eigenG object to save computation\n")
#	ncolZ=sapply(Z,function(a){ncola=ncol(a);ifelse(is.null(ncola),1,ncola)})
#	whichMaxcol=which.max(ncolZ)
#	whichEigenG=ifelse(whichMaxcol,whichMaxcol,1)   
#	eigenG=getEigenG(Zg=Z[[whichEigenG]])
#	} else {
#		eigenG = Z[[whichEigenG]]
#	}
#	if (length(whichNa) > 0) {
#		eigenG$U1 = eigenG$U1[-whichNa, ]
#	}
#	if (length(Z) > 1) {
#		listZw = Z[-whichEigenG]
#	}