data(mouse)


test_that("Var and loglik",{
  Z1=mouse.X[,1:20]
  Z2=mouse.X[,21:40]
  fit0=fitNULL(Obesity.BMI~GENDER+.eigenG(mouse.eigenG)+.R(Z1),data=mouse.pheno)
  a= c(2.284943e-03, 1.221593e-07, 1.373424e-10)
  names(a)=c("var_e","mouse.G","Z1")
  expect_equal(fit0$Var,a,tolerance = 1e-5)
    })

test_that("fit more than two random effects",{
  Z1=mouse.X[,1:20]
  Z2=mouse.X[,21:40]
  Z3=mouse.X[,41:60]
  fit0=fitNULL(Obesity.BMI~GENDER+.R(Z1)+.R(Z2)+.R(Z3),data=mouse.pheno)
  a= c(2.722769e-03, 1.912356e-11, 1.444994e-12, 3.691089e-07)
  names(a)=c("var_e","Z1","Z2","Z3")
  expect_equal(fit0$Var,a,tolerance = 1e-5)
})