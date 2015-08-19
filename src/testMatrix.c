/*test in R
R CMD SHLIB matrix.c
set.seed(1);
A=matrix(rnorm(10),nrow=5,ncol=2);
B=matrix(rnorm(20),nrow=2,ncol=10);
tA=t(A);
tB=t(B);


setwd("/Users/lianlian/Dropbox/github/SKAT2/inst/Ccode")
dyn.load("matrix.so")
A=matrix(rnorm(10),nrow=10,ncol=10);
determinant(A)
.Call("Rdeterminant",A,1);
.Call("RGetEigenValue",A)
.Call("Rsweep_prod",A,c(1:5),0,1);
sweep(A,1,c(1:5),"*")
.Call("Rsweep_prod",A,c(1:2),1,1);
sweep(t(A),1,c(1:2),"*")
.Call("Rsweep_prod",A,c(1:5),1,2);

.Call("R_matprod",A,B,0,0);
.Call("R_matprod",tA,tB,1,1);
A*5
matrix(.Call("R_matprod2",as.vector(A),10,1,5,1,1,0,0),nrow=nrow(A),ncol=ncol(A));
matrix(.Call("R_matprod2",as.vector(A),10,1,5,1,1,0,0),nrow=nrow(A),ncol=ncol(A));

Matrix inverse
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
h8 <- hilbert(8); h8
sh8 <- solve(h8)
.Call("R_matinv",h8);
*/
//testing in R
SEXP Rdeterminant(SEXP Ain, SEXP R_useLog){
	PROTECT(Ain=AS_NUMERIC(Ain));
	int n=INTEGER(getAttrib(Ain,R_DimSymbol))[0];
	SEXP modulus;
	PROTECT(modulus=allocVector(REALSXP,1));
	int useLog=INTEGER_VALUE(R_useLog);
	REAL(modulus)[0]=determinant(REAL(Ain), n, useLog);
	UNPROTECT(2);
	return(modulus);
}


SEXP  RsEigenValue(SEXP R_X)
{
	PROTECT(R_X=AS_NUMERIC(R_X));
	int n=INTEGER(getAttrib(R_X,R_DimSymbol))[0];
	SEXP R_values;
	PROTECT(R_values=allocVector(REALSXP,n));
	sEigenValue(REAL(R_X), n, REAL(R_values));
	UNPROTECT(2);
	return(R_values);
}

SEXP Rsweep_prod(SEXP R_X, SEXP R_v, SEXP R_transX, SEXP R_MARGINopX){

PROTECT(R_X=AS_NUMERIC(R_X));
PROTECT(R_v=AS_NUMERIC(R_v));
SEXP Rdim=getAttrib(R_X,R_DimSymbol);
int nrX=INTEGER(Rdim)[0];
int ncX=INTEGER(Rdim)[1];
int transX=INTEGER_VALUE(R_transX);
int MARGINopX=INTEGER_VALUE(R_MARGINopX);
SEXP R_Z;
if(transX){
PROTECT(R_Z=allocMatrix(REALSXP,ncX,nrX));
}else{
PROTECT(R_Z=allocMatrix(REALSXP,nrX,ncX));
}

double *X=NUMERIC_POINTER(R_X);
double *v=NUMERIC_POINTER(R_v);
double *Z=NUMERIC_POINTER(R_Z);

char trans[3]="NT";

sweep_prod(X,v,Z,nrX,ncX,trans[transX],MARGINopX);


UNPROTECT(3);
return R_Z;
}
//this shows that dgemm can be directly used for matrix times scaler as long as you put the dimensions right.Because after all, every thing is just a scaler.
SEXP R_matprod2(SEXP R_x,  SEXP R_nrx, SEXP R_ncx, SEXP R_y, SEXP R_nry, SEXP R_ncy, SEXP R_transa, SEXP R_transb){
PROTECT(R_x=AS_NUMERIC(R_x));
PROTECT(R_y=AS_NUMERIC(R_y));
int dimx[2]={INTEGER_VALUE(R_nrx), INTEGER_VALUE(R_ncx)};//INTEGER(R_ncx)[0] will have to work with R_ncx=AS_INTEGER(R_ncx);
int dimy[2]={INTEGER_VALUE(R_nry), INTEGER_VALUE(R_ncy)};
int transa=INTEGER_VALUE(R_transa);
int transb=INTEGER_VALUE(R_transb);
//double alpha =NUMERIC_VALUE(R_alpha);
//double beta=NUMERIC_VALUE(R_beta);
double *x;
x=NUMERIC_POINTER(R_x);
double *y;
y=NUMERIC_POINTER(R_y);
SEXP R_z;
PROTECT(R_z=allocMatrix(REALSXP,dimx[transa],dimy[!transb]));

char trans[3]="NT";

double *z;
z=NUMERIC_POINTER(R_z);


matprod(x, dimx[0], dimx[1], y, dimy[0],  dimy[1],   z,  trans[transa],  trans[transb],1,0);

UNPROTECT(3);
return(R_z);


}

SEXP R_matprod(SEXP R_x, SEXP R_y,SEXP R_transa, SEXP R_transb)
{
PROTECT(R_x=AS_NUMERIC(R_x));
PROTECT(R_y=AS_NUMERIC(R_y));
SEXP Rdimx=getAttrib(R_x,R_DimSymbol);
SEXP Rdimy=getAttrib(R_y,R_DimSymbol);
int dimx[2]={INTEGER(Rdimx)[0], INTEGER(Rdimx)[1]};
int dimy[2]={INTEGER(Rdimy)[0], INTEGER(Rdimy)[1]};
int transa=INTEGER_VALUE(R_transa);
int transb=INTEGER_VALUE(R_transb);
double *x;
x=NUMERIC_POINTER(R_x);
double *y;
y=NUMERIC_POINTER(R_y);

SEXP R_z;


PROTECT(R_z=allocMatrix(REALSXP,dimx[transa],dimy[!transb]));

printf("dimension of Z is %d,%d\n",dimx[transa],dimy[!transb]);

char trans[3]="NT";

double *z;
z=NUMERIC_POINTER(R_z);

matprod(x, dimx[0], dimx[1], y, dimy[0],  dimy[1],   z,  trans[transa],  trans[transb],1,0);

UNPROTECT(3);
return(R_z);

}



SEXP R_matinv(SEXP R_A)
{
	PROTECT(R_A=AS_NUMERIC(R_A));
	SEXP Rdimx=getAttrib(R_A,R_DimSymbol);
	int M=INTEGER(Rdimx)[0];
	int N=INTEGER(Rdimx)[1];
	if(M!=N)error("must be square matrix\n");
	SEXP R_z;
	PROTECT(R_z=allocMatrix(REALSXP,N,N));
	double *z;
	z=NUMERIC_POINTER(R_z);
    memcpy(z, REAL(R_A), sizeof(double)*N * N);
	matinv(z,N);
	UNPROTECT(2);
	return(R_z);
}
	
	


