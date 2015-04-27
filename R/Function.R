##functions
#Impute marker genotypes by column mean
meanImpute=function(X){apply(X,2,function(a){if(any(is.na(a))){a[which(ia.na(a))]=mean(a)};return(a)})