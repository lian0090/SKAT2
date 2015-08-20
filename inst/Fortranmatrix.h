#ifndef FORTRANMATRIX_H
#define FORTRANMATRIX_H

//uses dgetrf, dgetri, dgemm and dsyevr from FORTRAN

//blas functions
/* DGEMM - perform one of the matrix-matrix operations    */
/* C := alpha*op( A )*op( B ) + beta*C */
extern void
dgemm_(const char *transa, const char *transb, const int *m,
       const int *n, const int *k, const double *alpha,
       const double *a, const int *lda,
       const double *b, const int *ldb,
       const double *beta, double *c, const int *ldc);




/* DGETRF - compute an LU factorization of a general M-by-N */
/* matrix A using partial pivoting with row interchanges */
extern void
dgetrf_(const int* m, const int* n, double* a, const int* lda,
                  int* ipiv, int* info);



/* DGETRI - compute the inverse of a matrix using the LU */
/* factorization computed by DGETRF */
extern void
dgetri_(const int* n, double* a, const int* lda,
                    int* ipiv, double* work, const int* lwork, int* info);





/* DSYEVR - compute all eigenvalues and, optionally, eigenvectors   */
/* of a real symmetric matrix A					   */
extern void
dsyevr_(const char *jobz, const char *range, const char *uplo,
                 const int *n, double *a, const int *lda,
                 const double *vl, const double *vu,
                 const int *il, const int *iu,
                 const double *abstol, int *m, double *w,
                 double *z, const int *ldz, int *isuppz,
                 double *work, const int *lwork,
                 int *iwork, const int *liwork,
                 int *info);


#endif
