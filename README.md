

## Functions implemented in the package
- `GWAS`: GWAS= function(null.formula, MxE.formula,  setsize=1, sets=NULL, methods = NULL)
  * null.formula: formula for NULL model 
  * MxE.formula: formula for marker by environment interaction
  * setsize: size of the marker window/set to be tested. 
  * sets: an integer vector the same size as the number of markers. 
- `fitNULL`
- `testZ`
- `testX`

## Basic call of the GWAS function
 ```R
 library(SKAT2)
 data(mouse)
pheno=mouse.pheno
y=pheno$Obesity.BMI

GWAS(y~pheno$GENDER,~mouse.X[,1:100])
 ```
It currently does not support NA values well. especially when eigenG is supplied and there is NA value. Cgamma is NA. (tXX,tXZt,... does not exist when kd>n, check the code for Cgamma)

 
## Examples

### data preparation
```R
library(SKAT2)
data(mouse)
#G=tcrossprod(scale(mouse.X,T,F))
#eigenG=getEigenG(G)
```

### Example 1: GWAS by SKAT/Score: no random component in the NULL model.
##compare results with the original SKAT package
###code from SKAT2
```R
GWAS(formula=Obesity.BMI~GENDER,
GxE.formula=~.R(mouse.X[,1:100]),
data=mouse.pheno,methods="SKAT")
           SKAT
[1,] 0.36181008
[2,] 0.48028832
[3,] 0.27704551
[4,] 0.01444503
[5,] 0.29171852
```
###code for SKAT package
```R
obj<-SKAT_Null_Model(Obesity.BMI~GENDER, 
out_type="C", data=mouse.pheno)
p=matrix(0,5,1)
for(i in 1:5){
Zi=mouse.X[,(i-1)*20+1:20]
p[i,]=SKAT(Zi, obj, weights=rep(1,ncol(Zi)),
 is_check_genotype=F)$p.value   }
> p
           [,1]
[1,] 0.36181008
[2,] 0.48028832
[3,] 0.27704551
[4,] 0.01444503
[5,] 0.29171852
```

```R
GWAS(formula=Obesity.BMI~GENDER,GxE.formula=~.R(mouse.X[,1:100]),data=mouse.pheno,methods=c("Score","SKAT"))
            Score       SKAT
[1,] 5.332589e-01 0.36181008
[2,] 7.375804e-01 0.48028832
[3,] 3.119485e-01 0.27704551
[4,] 7.072944e-12 0.01444503
[5,] 3.163844e-01 0.29171852
```

###Example 2: GWAS by Score: fit cage as random effect in the in the NULL model.
```R
mouse.pheno$cage=as.factor(mouse.pheno$cage)
GWAS(formula=Obesity.BMI~GENDER+.R(cage),GxE.formula=~.R(mouse.X[,1:100]),data=mouse.pheno)
           Score
[1,] 0.598052536
[2,] 0.408630191
[3,] 0.424617959
[4,] 0.000140982
[5,] 0.474230337
```

###Example 3: GWAS by Score: fit Genomic background in the NULL model
.G(G) tells the model that G is a Genomic matrix
```R
##cage is a 525 level radnom variable. If we fit both .R(cage) and eigenG, it will be extremely slow. Now, I will only fit eigenG, no cage.
pheno$cage=as.factor(pheno$cage)
#GWAS(formula=Obesity.BMI~GENDER+.G(mouse.G),GxE.formula=~(1|mouse.X[,1:100]),data=pheno)
GWAS(formula=Obesity.BMI~GENDER+.eigenG(mouse.eigenG),GxE.formula=~.R(mouse.X[,1:100]),data=mouse.pheno)
         Score
[1,] 0.6469787
[2,] 0.4116452
[3,] 0.8181726
[4,] 0.6267310
[5,] 0.6903958
```
###Example 4: GWAS by SKAT/Score: GxE term
```R
> GWAS(formula=Obesity.BMI~GENDER+.eigenG(mouse.eigenG),GxE.formula=~.R(mouse.X[,1:100]:GENDER),data=mouse.pheno)

           Score
[1,] 0.337835569
[2,] 0.003583466
[3,] 0.299666718
[4,] 0.797022229
[5,] 0.348932006

GWAS(formula=Obesity.BMI~GENDER+.eigenG(mouse.eigenG),GxE.formula=~.R(mouse.X[,1:100])+.R(mouse.X[,1:100]:GENDER),data=mouse.pheno)
           Score
[1,] 0.337832101
[2,] 0.003583462
[3,] 0.299666219
[4,] 0.797022323
[5,] 0.348931831
```

###Example 5: GWAS single markers
```
pvalue=GWAS(formula=Obesity.BMI~GENDER+.eigenG(mouse.eigenG),GxE.formula=~mouse.X[,1:100],data=mouse.pheno,setsize=1)
> head(pvalue$p.value)
     SSNP.P3D.LR
[1,]   0.5685092
[2,]   0.5962195
[3,]   0.9435306
[4,]   0.7538891
[5,]   0.1313761
[6,]   0.5962195
```



##Work on single marker windows. More general cases. 
###test the window of SNPs (Z1)
```R
Z1=mouse.X[,1:20]
Z2=mouse.X[,21:40]
fit0=fitNULL(Obesity.BMI ~GENDER,data=mouse.pheno)
out=testZ(fit0,Zt=Z1,method="SKAT")
> out
$Score
NULL

$sdScore
NULL

$p.value
    Score      SKAT 
       NA 0.3618101 

$Q
[1] 6806.817
library(SKAT2)
data(mouse)
Z1=mouse.X[,1:2]
Z2=mouse.X[,21:40]
fit0=fitNULL(Obesity.BMI~GENDER+.eigenG(mouse.eigenG)+(1|Z1),data=mouse.pheno)
testZ(fit0,Z2,methods=c("Score","SKAT"))
$Score
[1] 324618.7

$sdScore
[1] 1454048

$p.value
    Score      SKAT 
0.4116698 0.2620787 

$Q
[1] 1676088
library(SKAT2)
data(mouse)
Z1=mouse.X[,1:2]
Z2=mouse.X[,21:40]
testMat=model.matrix(~-1+Z1:GENDER,data=mouse.pheno)
####As long as I do not fit Litter together with another random effect, everything works fine
fit0=fitNULL(Obesity.BMI~GENDER+rnorm(nrow(Z1))+.eigenG(mouse.eigenG)+(1|Z1),data=mouse.pheno)
testZ(fit0,testMat,methods=c("Score","SKAT"))
###But whenever I fit Litter together with another random effects, R stops: memory not mapped.
fit0=fitNULL(Obesity.BMI~Litter+.eigenG(mouse.eigenG)+(1|Z1),data=mouse.pheno)
testZ(fit0,testMat,methods=c("Score","SKAT"))

```


## Test markers as fixed effects with P3D method (population structure previously determined)



### Step 1,  fit the NULL model to get the population parameters
```R
fit0=fitNULL(y~X1+.eigenG(eigenG))
testX(fit0,Xt=Z1[,1])
```
Output
```R
$p.value
[1] 0.8228676

> fit0=fitNULL(y~X1+.eigenG(eigenG)+(1|Z1))
> testX(fit0,Z2[,1])
$LR
 SSNP.P3D.LR 
7.527054e-05 

$p.value
SSNP.P3D.LR 
  0.9930778 

$logML0
SSNP.P3D.LR 
     628.26 

$logML1
SSNP.P3D.LR 
     628.26 

$Var0
$Var0$SSNP.P3D.LR
        VarE         VarG        VarW1 
2.273050e-03 7.077773e-07 2.251884e+02 


$Var1
$Var1$SSNP.P3D.LR
        VarE         VarG        VarW1 
2.273050e-03 7.077773e-07 2.251884e+02 
```


## Installation

SKAT2 is not available on [CRAN](http://cran.r-project.org/) yet. However, it can be installed directly from GitHub using the [devtools](https://github.com/hadley/devtools) package.

1. Install `devtools` package: `install.packages('devtools')`
2. Load `devtools` package: `library(devtools)`
3. Install `BGData` package from GitHub: `install_github('lian0090/SKAT2')`
4. Load `SKAT2` package: `library(SKAT2)`

## Possible Problems during installation
### OS X Mavericks
error ld: library not found for -lgfortran 
Solutions:
open terminal and type:
```
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```
Explanations can be found here 
http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/
