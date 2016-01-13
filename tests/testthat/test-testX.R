data(mouse)


test_that("an eigenG object also named eigenG, this should work well",{
  n=500
  p=20
  pheno=mouse.pheno[1:n,]
  Z1=mouse.X[1:n,1:p]
  X=cbind(model.matrix(~pheno$GENDER)[,-1],pheno$CageDensity)
  y=pheno$Obesity.BMI
  G=tcrossprod(scale(mouse.X[1:n,],T,F))
  eigenG=getEigenG(G=G)
  fit0=fitNULL(y~X+.eigenG(eigenG),data=pheno)
  out=testX(fit0,Z1[,1])
  a= 0.8228676
  names(a)="SSNP.P3D.LR"
  expect_equal(out$p.value, a,tolerance = 1e-5)
  
    })

