##functions
#Impute marker genotypes by column mean
meanImpute=function(X){apply(X,2,function(a){if(any(is.na(a))){a[which(ia.na(a))]=mean(a)};return(a)})
	
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
	