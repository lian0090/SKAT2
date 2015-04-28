##functions
#Impute marker genotypes by column mean
meanImpute=function(X){
	X=apply(X,2,function(a){if(any(is.na(a))){a[which(is.na(a))]=mean(a,na.rm=T)};return(a)})
	return(X)
	}
	
#Get columwise multiplication 
colmult=function(Z1,Z2){
		p2=ncol(Z1)*ncol(Z2)
	Z3=matrix(nrow=nrow(Z1),ncol=p2)

	for(j in 1:ncol(Z2)){
	for(i in c(1:ncol(Z1))){
	Z3[,(j-1)*ncol(Z1)+i]=Z1[,i]*Z2[,j]	
		}
	}
	return(Z3)
	}	

#calculate matrix trace
tr=function(X){
	out=sum(diag(X))
	return(out)
}	

#simulate beta for random effects
simuBeta=function(Z,k,var_e=1){
	Z=meanImpute(Z)
	Z=scale(Z,T,F)
    sumvar=sum(apply(Z,2,var))
    beta=rnorm(ncol(Z),0,sqrt(k*var_e/sumvar))
    u=Z%*%beta
    return(list(Z=Z,beta=beta,u=u))
}
