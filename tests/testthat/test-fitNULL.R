data(mouse)


test_that("Var and loglik",{
  Z1=mouse.X[,1:20]
  Z2=mouse.X[,21:40]
  fit0=fitNULL(Obesity.BMI~GENDER+.eigenG(mouse.eigenG)+.R(Z1),data=mouse.pheno)
  a= c(2.284944e-03,1.221591e-07)
  names(a)=c("mouse.G","Z1")
  expect_equal(fit0$Var,a,tolerance = 1e-5)
    })

test_that("fit more than two random effects",{
  Z1=mouse.X[,1:20]
  Z2=mouse.X[,21:40]
  Z3=mouse.X[,41:60]
  fit0=fitNULL(Obesity.BMI~GENDER+.R(Z1)+.R(Z2)+.R(Z3),data=mouse.pheno)
  a= c(0.002735162, 0.799769622, 0.630553858)
  names(a)=c("Z1","Z2","Z3")
  expect_equal(fit0$Var,a,tolerance = 1e-5)
})