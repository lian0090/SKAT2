#ifndef MATRIX_H
#define MATRIX_H

//A will be written over with Ainv
 void matinv(double *A,int N);

//Z=sweep(op(X),v,MARGINopX,"*");
//transX='N' or 'T'
 void sweep_prod(double *X, double *v, double *Z, int nrX, int ncX, int transX, int MARGINopX);
//z :=      alpha*op( x)*op( y ) + beta*z,
//when beta  is supplied as zero then C need not be set on input
////transa,transb='N' or 'T'
 void matprod(double *x, int nrx, int ncx,
		    double *y, int nry, int ncy, double *z, char transa, char transb,double alpha, double beta);
 void sEigenValue(double *x, int n, double *values);
 double determinant(double *X, int n, int useLog);

#endif
