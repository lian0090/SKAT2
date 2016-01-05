data(mouse)


test_that("No random effects in NULL model",{
  Z1=mouse.X[,1:20]
  Z2=mouse.X[,21:40]
  fit0=fitNULL(Obesity.BMI ~GENDER,data=mouse.pheno)
  out=testZ(fit0,Ztest=Z1,method="SKAT")
  a=0.3618101
  names(a)=c("SKAT")
  expect_equal(out$p.value[2],a,tolerance=1e-5)
  a=6806.817
  expect_equal(out$Q,a,tolerance=1e-5)
    })

test_that("fit more than two random effects",{
  Z1=mouse.X[,1:20]
  Z2=mouse.X[,21:40]
  fit0=fitNULL(Obesity.BMI~GENDER+.eigenG(mouse.eigenG)+.R(Z1),data=mouse.pheno)
  out=testZ(fit0,Z2,methods=c("Score","SKAT"))
  p.value=c(0.4116452,0.2620584)
  names(p.value)=c("Score","SKAT")
  a= list(Score=324760.8,sdScore=1454273,p.value=p.value,Q=1676405)
  expect_equal(out,a,tolerance = 1e-5)
})