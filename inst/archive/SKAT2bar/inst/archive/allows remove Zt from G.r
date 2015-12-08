#GWAS on a marker windows .
GWAS= function(y, X0 = NA, Z0 = NA, M, sets = rep(c(1:ceiling(ncol(M)/10)), each = 10)[1:ncol(M)], methods = c("Score", "SKAT")) {
    FaST = getFaST(y, X = X0, Z = Z, Zt = NA, tauRel = tauRel)
    if(removeMsetFromZ0){
        if(!is.na(FaST$U1)){
        nw=FaST$nw+1
        Zw = c(Z0, list(Zt[,1:2]))
        
    }

    }

    if (is.na(Z0)) {
        if (removeZtFromZ0 == T) {
            nw = length(Z0) - 1
            Z = c(Z0, list(Zti))
            tauRel = paste("tauw", nw, "=-taud", sep = "")
        } else {
            Z = Z0
            tauRel = NA
        }
    } else {
        Z = NA
        tauRel = NA
    }

    
    nSets = length(unique(sets))
    nmethods = length(methods)
    out = matrix(nrow = nSets, ncol = nmethods)
    for (i in 1:nSets) {
        Zti = M[, which(sets == i)]
        if(!is.na(FaST$U1) & removeMsetFromZ0){
            FaST=updateFaST.Zw(FaST,FaST$listZw)
            Z = c(Z0, list(Zti))
            tauRel = paste("tauw", nw, "=-taud", sep = "")



        }
        FaST = updateFaST.Zt(FaST, Zt = Zti)
        out[i, ] = testZ.FaST(FaST, windowtest = methods)[paste("p.", methods, sep = "")]
    }
    return(out)
}
