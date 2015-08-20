*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER K,LDA,LDB,LDC,M,N
*       CHARACTER TRANSA,TRANSB
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGEMM  performs one of the matrix-matrix operations
*>
*>    C := alpha*op( A )*op( B ) + beta*C,
*>
*> where  op( X ) is one of
*>
*>    op( X ) = X   or   op( X ) = X**T,
*>
*> alpha and beta are scalars, and A, B and C are matrices, with op( A )
*> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n',  op( A ) = A.
*>
*>              TRANSA = 'T' or 't',  op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c',  op( A ) = A**T.
*> \endverbatim
*>
*> \param[in] TRANSB
*> \verbatim
*>          TRANSB is CHARACTER*1
*>           On entry, TRANSB specifies the form of op( B ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSB = 'N' or 'n',  op( B ) = B.
*>
*>              TRANSB = 'T' or 't',  op( B ) = B**T.
*>
*>              TRANSB = 'C' or 'c',  op( B ) = B**T.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry,  M  specifies  the number  of rows  of the  matrix
*>           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry,  N  specifies the number  of columns of the matrix
*>           op( B ) and the number of columns of the matrix C. N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>           On entry,  K  specifies  the number of columns of the matrix
*>           op( A ) and the number of rows of the matrix op( B ). K must
*>           be at least  zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*>           part of the array  A  must contain the matrix  A,  otherwise
*>           the leading  k by m  part of the array  A  must contain  the
*>           matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*>           least  max( 1, k ).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*>           part of the array  B  must contain the matrix  B,  otherwise
*>           the leading  n by k  part of the array  B  must contain  the
*>           matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*>           least  max( 1, n ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION.
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*>           supplied as zero then C need not be set on input.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*>           Before entry, the leading  m by n  part of the array  C must
*>           contain the matrix  C,  except when  beta  is zero, in which
*>           case C need not be set on entry.
*>           On exit, the array  C  is overwritten by the  m by n  matrix
*>           ( alpha*op( A )*op( B ) + beta*C ).
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>           On entry, LDC specifies the first dimension of C as declared
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2011
*
*> \ingroup double_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
SUBROUTINE dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
DOUBLE PRECISION alpha,beta
INTEGER k,lda,ldb,ldc,m,n
CHARACTER transa,transb
*     ..
*     .. Array Arguments ..
DOUBLE PRECISION a(lda,*),b(ldb,*),c(ldc,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
LOGICAL lsame
EXTERNAL lsame
*     ..
*     .. External Subroutines ..
EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
INTRINSIC max
*     ..
*     .. Local Scalars ..
DOUBLE PRECISION temp
INTEGER i,info,j,l,ncola,nrowa,nrowb
LOGICAL nota,notb
*     ..
*     .. Parameters ..
DOUBLE PRECISION one,zero
parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
nota = lsame(transa,'N')
notb = lsame(transb,'N')
IF (nota) THEN
nrowa = m
ncola = k
ELSE
nrowa = k
ncola = m
END IF
IF (notb) THEN
nrowb = k
ELSE
nrowb = n
END IF
*
*     Test the input parameters.
*
info = 0
IF ((.NOT.nota) .AND. (.NOT.lsame(transa,'C')) .AND.
+    (.NOT.lsame(transa,'T'))) THEN
info = 1
ELSE IF ((.NOT.notb) .AND. (.NOT.lsame(transb,'C')) .AND.
+         (.NOT.lsame(transb,'T'))) THEN
info = 2
ELSE IF (m.LT.0) THEN
info = 3
ELSE IF (n.LT.0) THEN
info = 4
ELSE IF (k.LT.0) THEN
info = 5
ELSE IF (lda.LT.max(1,nrowa)) THEN
info = 8
ELSE IF (ldb.LT.max(1,nrowb)) THEN
info = 10
ELSE IF (ldc.LT.max(1,m)) THEN
info = 13
END IF
IF (info.NE.0) THEN
CALL xerbla('DGEMM ',info)
RETURN
END IF
*
*     Quick return if possible.
*
IF ((m.EQ.0) .OR. (n.EQ.0) .OR.
+    (((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one))) RETURN
*
*     And if  alpha.eq.zero.
*
IF (alpha.EQ.zero) THEN
IF (beta.EQ.zero) THEN
DO 20 j = 1,n
DO 10 i = 1,m
c(i,j) = zero
10             CONTINUE
20         CONTINUE
ELSE
DO 40 j = 1,n
DO 30 i = 1,m
c(i,j) = beta*c(i,j)
30             CONTINUE
40         CONTINUE
END IF
RETURN
END IF
*
*     Start the operations.
*
IF (notb) THEN
IF (nota) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
DO 90 j = 1,n
IF (beta.EQ.zero) THEN
DO 50 i = 1,m
c(i,j) = zero
50                 CONTINUE
ELSE IF (beta.NE.one) THEN
DO 60 i = 1,m
c(i,j) = beta*c(i,j)
60                 CONTINUE
END IF
DO 80 l = 1,k
IF (b(l,j).NE.zero) THEN
temp = alpha*b(l,j)
DO 70 i = 1,m
c(i,j) = c(i,j) + temp*a(i,l)
70                     CONTINUE
END IF
80             CONTINUE
90         CONTINUE
ELSE
*
*           Form  C := alpha*A**T*B + beta*C
*
DO 120 j = 1,n
DO 110 i = 1,m
temp = zero
DO 100 l = 1,k
temp = temp + a(l,i)*b(l,j)
100                 CONTINUE
IF (beta.EQ.zero) THEN
c(i,j) = alpha*temp
ELSE
c(i,j) = alpha*temp + beta*c(i,j)
END IF
110             CONTINUE
120         CONTINUE
END IF
ELSE
IF (nota) THEN
*
*           Form  C := alpha*A*B**T + beta*C
*
DO 170 j = 1,n
IF (beta.EQ.zero) THEN
DO 130 i = 1,m
c(i,j) = zero
130                 CONTINUE
ELSE IF (beta.NE.one) THEN
DO 140 i = 1,m
c(i,j) = beta*c(i,j)
140                 CONTINUE
END IF
DO 160 l = 1,k
IF (b(j,l).NE.zero) THEN
temp = alpha*b(j,l)
DO 150 i = 1,m
c(i,j) = c(i,j) + temp*a(i,l)
150                     CONTINUE
END IF
160             CONTINUE
170         CONTINUE
ELSE
*
*           Form  C := alpha*A**T*B**T + beta*C
*
DO 200 j = 1,n
DO 190 i = 1,m
temp = zero
DO 180 l = 1,k
temp = temp + a(l,i)*b(j,l)
180                 CONTINUE
IF (beta.EQ.zero) THEN
c(i,j) = alpha*temp
ELSE
c(i,j) = alpha*temp + beta*c(i,j)
END IF
190             CONTINUE
200         CONTINUE
END IF
END IF
*
RETURN
*
*     End of DGEMM .
*
END


SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
INTEGER            IPIV( * )
DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
DOUBLE PRECISION   ONE
PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. External Functions ..
INTEGER            ILAENV
EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
INFO = 0
IF( M.LT.0 ) THEN
INFO = -1
ELSE IF( N.LT.0 ) THEN
INFO = -2
ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
INFO = -4
END IF
IF( INFO.NE.0 ) THEN
CALL XERBLA( 'DGETRF', -INFO )
RETURN
END IF
*
*     Quick return if possible
*
IF( M.EQ.0 .OR. N.EQ.0 )
$   RETURN
*
*     Determine the block size for this environment.
*
NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
CALL DGETF2( M, N, A, LDA, IPIV, INFO )
ELSE
*
*        Use blocked code.
*
DO 20 J = 1, MIN( M, N ), NB
JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
IF( INFO.EQ.0 .AND. IINFO.GT.0 )
$         INFO = IINFO + J - 1
DO 10 I = J, MIN( M, J+JB-1 )
IPIV( I ) = J - 1 + IPIV( I )
10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
$                      IPIV, 1 )
*
*              Compute block row of U.
*
CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
$                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
$                     LDA )
IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
$                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
$                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
$                        LDA )
END IF
END IF
20    CONTINUE
END IF
RETURN
*
*     End of DGETRF
*
END




SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
INTEGER            IPIV( * )
DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRI computes the inverse of a matrix using the LU factorization
*  computed by DGETRF.
*
*  This method inverts U and then computes inv(A) by solving the system
*  inv(A)*L = inv(U) for inv(A).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the factors L and U from the factorization
*          A = P*L*U as computed by DGETRF.
*          On exit, if INFO = 0, the inverse of the original matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimal performance LWORK >= N*NB, where NB is
*          the optimal blocksize returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
*                singular and its inverse could not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
DOUBLE PRECISION   ZERO, ONE
PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
LOGICAL            LQUERY
INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB,
$                   NBMIN, NN
*     ..
*     .. External Functions ..
INTEGER            ILAENV
EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
EXTERNAL           DGEMM, DGEMV, DSWAP, DTRSM, DTRTRI, XERBLA
*     ..
*     .. Intrinsic Functions ..
INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
INFO = 0
NB = ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 )
LWKOPT = N*NB
WORK( 1 ) = LWKOPT
LQUERY = ( LWORK.EQ.-1 )
IF( N.LT.0 ) THEN
INFO = -1
ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
INFO = -3
ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
INFO = -6
END IF
IF( INFO.NE.0 ) THEN
CALL XERBLA( 'DGETRI', -INFO )
RETURN
ELSE IF( LQUERY ) THEN
RETURN
END IF
*
*     Quick return if possible
*
IF( N.EQ.0 )
$   RETURN
*
*     Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,
*     and the inverse is not computed.
*
CALL DTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
IF( INFO.GT.0 )
$   RETURN
*
NBMIN = 2
LDWORK = N
IF( NB.GT.1 .AND. NB.LT.N ) THEN
IWS = MAX( LDWORK*NB, 1 )
IF( LWORK.LT.IWS ) THEN
NB = LWORK / LDWORK
NBMIN = MAX( 2, ILAENV( 2, 'DGETRI', ' ', N, -1, -1, -1 ) )
END IF
ELSE
IWS = N
END IF
*
*     Solve the equation inv(A)*L = inv(U) for inv(A).
*
IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
DO 20 J = N, 1, -1
*
*           Copy current column of L to WORK and replace with zeros.
*
DO 10 I = J + 1, N
WORK( I ) = A( I, J )
A( I, J ) = ZERO
10       CONTINUE
*
*           Compute current column of inv(A).
*
IF( J.LT.N )
$         CALL DGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ),
$                     LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
20    CONTINUE
ELSE
*
*        Use blocked code.
*
NN = ( ( N-1 ) / NB )*NB + 1
DO 50 J = NN, 1, -NB
JB = MIN( NB, N-J+1 )
*
*           Copy current block column of L to WORK and replace with
*           zeros.
*
DO 40 JJ = J, J + JB - 1
DO 30 I = JJ + 1, N
WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
A( I, JJ ) = ZERO
30          CONTINUE
40       CONTINUE
*
*           Compute current block column of inv(A).
*
IF( J+JB.LE.N )
$         CALL DGEMM( 'No transpose', 'No transpose', N, JB,
$                     N-J-JB+1, -ONE, A( 1, J+JB ), LDA,
$                     WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
CALL DTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
$                  ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
50    CONTINUE
END IF
*
*     Apply column interchanges.
*
DO 60 J = N - 1, 1, -1
JP = IPIV( J )
IF( JP.NE.J )
$      CALL DSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
60 CONTINUE
*
WORK( 1 ) = IWS
RETURN
*
*     End of DGETRI
*
END



SUBROUTINE DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
$                   ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,
$                   IWORK, LIWORK, INFO )
*
*  -- LAPACK driver routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
CHARACTER          JOBZ, RANGE, UPLO
INTEGER            IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N
DOUBLE PRECISION   ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
INTEGER            ISUPPZ( * ), IWORK( * )
DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DSYEVR computes selected eigenvalues and, optionally, eigenvectors
*  of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
*  selected by specifying either a range of values or a range of
*  indices for the desired eigenvalues.
*
*  DSYEVR first reduces the matrix A to tridiagonal form T with a call
*  to DSYTRD.  Then, whenever possible, DSYEVR calls DSTEMR to compute
*  the eigenspectrum using Relatively Robust Representations.  DSTEMR
*  computes eigenvalues by the dqds algorithm, while orthogonal
*  eigenvectors are computed from various "good" L D L^T representations
*  (also known as Relatively Robust Representations). Gram-Schmidt
*  orthogonalization is avoided as far as possible. More specifically,
*  the various steps of the algorithm are as follows.
*
*  For each unreduced block (submatrix) of T,
*     (a) Compute T - sigma I  = L D L^T, so that L and D
*         define all the wanted eigenvalues to high relative accuracy.
*         This means that small relative changes in the entries of D and L
*         cause only small relative changes in the eigenvalues and
*         eigenvectors. The standard (unfactored) representation of the
*         tridiagonal matrix T does not have this property in general.
*     (b) Compute the eigenvalues to suitable accuracy.
*         If the eigenvectors are desired, the algorithm attains full
*         accuracy of the computed eigenvalues only right before
*         the corresponding vectors have to be computed, see steps c) and d).
*     (c) For each cluster of close eigenvalues, select a new
*         shift close to the cluster, find a new factorization, and refine
*         the shifted eigenvalues to suitable accuracy.
*     (d) For each eigenvalue with a large enough relative separation compute
*         the corresponding eigenvector by forming a rank revealing twisted
*         factorization. Go back to (c) for any clusters that remain.
*
*  The desired accuracy of the output can be specified by the input
*  parameter ABSTOL.
*
*  For more details, see DSTEMR's documentation and:
*  - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations
*    to compute orthogonal eigenvectors of symmetric tridiagonal matrices,"
*    Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.
*  - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and
*    Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25,
*    2004.  Also LAPACK Working Note 154.
*  - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric
*    tridiagonal eigenvalue/eigenvector problem",
*    Computer Science Division Technical Report No. UCB/CSD-97-971,
*    UC Berkeley, May 1997.
*
*
*  Note 1 : DSYEVR calls DSTEMR when the full spectrum is requested
*  on machines which conform to the ieee-754 floating point standard.
*  DSYEVR calls DSTEBZ and SSTEIN on non-ieee machines and
*  when partial spectrum requests are made.
*
*  Normal execution of DSTEMR may create NaNs and infinities and
*  hence may abort due to a floating point exception in environments
*  which do not handle NaNs and infinities in the ieee standard default
*  manner.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  RANGE   (input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the half-open interval (VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
********** For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and
********** DSTEIN are called
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, the lower triangle (if UPLO='L') or the upper
*          triangle (if UPLO='U') of A, including the diagonal, is
*          destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  VL      (input) DOUBLE PRECISION
*  VU      (input) DOUBLE PRECISION
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues. VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  ABSTOL  (input) DOUBLE PRECISION
*          The absolute error tolerance for the eigenvalues.
*          An approximate eigenvalue is accepted as converged
*          when it is determined to lie in an interval [a,b]
*          of width less than or equal to
*
*                  ABSTOL + EPS *   max( |a|,|b| ) ,
*
*          where EPS is the machine precision.  If ABSTOL is less than
*          or equal to zero, then  EPS*|T|  will be used in its place,
*          where |T| is the 1-norm of the tridiagonal matrix obtained
*          by reducing A to tridiagonal form.
*
*          See "Computing Small Singular Values of Bidiagonal Matrices
*          with Guaranteed High Relative Accuracy," by Demmel and
*          Kahan, LAPACK Working Note #3.
*
*          If high relative accuracy is important, set ABSTOL to
*          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that
*          eigenvalues are computed to high relative accuracy when
*          possible in future releases.  The current code does not
*          make any guarantees about high relative accuracy, but
*          future releases will. See J. Barlow and J. Demmel,
*          "Computing Accurate Eigensystems of Scaled Diagonally
*          Dominant Matrices", LAPACK Working Note #7, for a discussion
*          of which matrices define their eigenvalues to high relative
*          accuracy.
*
*  M       (output) INTEGER
*          The total number of eigenvalues found.  0 <= M <= N.
*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          The first M elements contain the selected eigenvalues in
*          ascending order.
*
*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))
*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
*          contain the orthonormal eigenvectors of the matrix A
*          corresponding to the selected eigenvalues, with the i-th
*          column of Z holding the eigenvector associated with W(i).
*          If JOBZ = 'N', then Z is not referenced.
*          Note: the user must ensure that at least max(1,M) columns are
*          supplied in the array Z; if RANGE = 'V', the exact value of M
*          is not known in advance and an upper bound must be used.
*          Supplying N columns is always safe.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  ISUPPZ  (output) INTEGER array, dimension ( 2*max(1,M) )
*          The support of the eigenvectors in Z, i.e., the indices
*          indicating the nonzero elements in Z. The i-th eigenvector
*          is nonzero only in elements ISUPPZ( 2*i-1 ) through
*          ISUPPZ( 2*i ).
********** Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,26*N).
*          For optimal efficiency, LWORK >= (NB+6)*N,
*          where NB is the max of the blocksize for DSYTRD and DORMTR
*          returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.  LIWORK >= max(1,10*N).
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  Internal error
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Inderjit Dhillon, IBM Almaden, USA
*     Osni Marques, LBNL/NERSC, USA
*     Ken Stanley, Computer Science Division, University of
*       California at Berkeley, USA
*     Jason Riedy, Computer Science Division, University of
*       California at Berkeley, USA
*
* =====================================================================
*
*     .. Parameters ..
DOUBLE PRECISION   ZERO, ONE, TWO
PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
LOGICAL            ALLEIG, INDEIG, LOWER, LQUERY, VALEIG, WANTZ,
$                   TRYRAC
CHARACTER          ORDER
INTEGER            I, IEEEOK, IINFO, IMAX, INDD, INDDD, INDE,
$                   INDEE, INDIBL, INDIFL, INDISP, INDIWO, INDTAU,
$                   INDWK, INDWKN, ISCALE, J, JJ, LIWMIN,
$                   LLWORK, LLWRKN, LWKOPT, LWMIN, NB, NSPLIT
DOUBLE PRECISION   ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN,
$                   SIGMA, SMLNUM, TMP1, VLL, VUU
*     ..
*     .. External Functions ..
LOGICAL            LSAME
INTEGER            ILAENV
DOUBLE PRECISION   DLAMCH, DLANSY
EXTERNAL           LSAME, ILAENV, DLAMCH, DLANSY
*     ..
*     .. External Subroutines ..
EXTERNAL           DCOPY, DORMTR, DSCAL, DSTEBZ, DSTEMR, DSTEIN,
$                   DSTERF, DSWAP, DSYTRD, XERBLA
*     ..
*     .. Intrinsic Functions ..
INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
IEEEOK = ILAENV( 10, 'DSYEVR', 'N', 1, 2, 3, 4 )
*
LOWER = LSAME( UPLO, 'L' )
WANTZ = LSAME( JOBZ, 'V' )
ALLEIG = LSAME( RANGE, 'A' )
VALEIG = LSAME( RANGE, 'V' )
INDEIG = LSAME( RANGE, 'I' )
*
LQUERY = ( ( LWORK.EQ.-1 ) .OR. ( LIWORK.EQ.-1 ) )
*
LWMIN = MAX( 1, 26*N )
LIWMIN = MAX( 1, 10*N )
*
INFO = 0
IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
INFO = -1
ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
INFO = -2
ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
INFO = -3
ELSE IF( N.LT.0 ) THEN
INFO = -4
ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
INFO = -6
ELSE
IF( VALEIG ) THEN
IF( N.GT.0 .AND. VU.LE.VL )
$         INFO = -8
ELSE IF( INDEIG ) THEN
IF( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) THEN
INFO = -9
ELSE IF( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) THEN
INFO = -10
END IF
END IF
END IF
IF( INFO.EQ.0 ) THEN
IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
INFO = -15
ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
INFO = -18
ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
INFO = -20
END IF
END IF
*
IF( INFO.EQ.0 ) THEN
NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
NB = MAX( NB, ILAENV( 1, 'DORMTR', UPLO, N, -1, -1, -1 ) )
LWKOPT = MAX( ( NB+1 )*N, LWMIN )
WORK( 1 ) = LWKOPT
IWORK( 1 ) = LIWMIN
END IF
*
IF( INFO.NE.0 ) THEN
CALL XERBLA( 'DSYEVR', -INFO )
RETURN
ELSE IF( LQUERY ) THEN
RETURN
END IF
*
*     Quick return if possible
*
M = 0
IF( N.EQ.0 ) THEN
WORK( 1 ) = 1
RETURN
END IF
*
IF( N.EQ.1 ) THEN
WORK( 1 ) = 7
IF( ALLEIG .OR. INDEIG ) THEN
M = 1
W( 1 ) = A( 1, 1 )
ELSE
IF( VL.LT.A( 1, 1 ) .AND. VU.GE.A( 1, 1 ) ) THEN
M = 1
W( 1 ) = A( 1, 1 )
END IF
END IF
IF( WANTZ )
$      Z( 1, 1 ) = ONE
RETURN
END IF
*
*     Get machine constants.
*
SAFMIN = DLAMCH( 'Safe minimum' )
EPS = DLAMCH( 'Precision' )
SMLNUM = SAFMIN / EPS
BIGNUM = ONE / SMLNUM
RMIN = SQRT( SMLNUM )
RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
*     Scale matrix to allowable range, if necessary.
*
ISCALE = 0
ABSTLL = ABSTOL
VLL = VL
VUU = VU
ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
ISCALE = 1
SIGMA = RMIN / ANRM
ELSE IF( ANRM.GT.RMAX ) THEN
ISCALE = 1
SIGMA = RMAX / ANRM
END IF
IF( ISCALE.EQ.1 ) THEN
IF( LOWER ) THEN
DO 10 J = 1, N
CALL DSCAL( N-J+1, SIGMA, A( J, J ), 1 )
10       CONTINUE
ELSE
DO 20 J = 1, N
CALL DSCAL( J, SIGMA, A( 1, J ), 1 )
20       CONTINUE
END IF
IF( ABSTOL.GT.0 )
$      ABSTLL = ABSTOL*SIGMA
IF( VALEIG ) THEN
VLL = VL*SIGMA
VUU = VU*SIGMA
END IF
END IF

*     Initialize indices into workspaces.  Note: The IWORK indices are
*     used only if DSTERF or DSTEMR fail.

*     WORK(INDTAU:INDTAU+N-1) stores the scalar factors of the
*     elementary reflectors used in DSYTRD.
INDTAU = 1
*     WORK(INDD:INDD+N-1) stores the tridiagonal's diagonal entries.
INDD = INDTAU + N
*     WORK(INDE:INDE+N-1) stores the off-diagonal entries of the
*     tridiagonal matrix from DSYTRD.
INDE = INDD + N
*     WORK(INDDD:INDDD+N-1) is a copy of the diagonal entries over
*     -written by DSTEMR (the DSTERF path copies the diagonal to W).
INDDD = INDE + N
*     WORK(INDEE:INDEE+N-1) is a copy of the off-diagonal entries over
*     -written while computing the eigenvalues in DSTERF and DSTEMR.
INDEE = INDDD + N
*     INDWK is the starting offset of the left-over workspace, and
*     LLWORK is the remaining workspace size.
INDWK = INDEE + N
LLWORK = LWORK - INDWK + 1

*     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in DSTEBZ and
*     stores the block indices of each of the M<=N eigenvalues.
INDIBL = 1
*     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in DSTEBZ and
*     stores the starting and finishing indices of each block.
INDISP = INDIBL + N
*     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
*     that corresponding to eigenvectors that fail to converge in
*     DSTEIN.  This information is discarded; if any fail, the driver
*     returns INFO > 0.
INDIFL = INDISP + N
*     INDIWO is the offset of the remaining integer workspace.
INDIWO = INDISP + N

*
*     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
*
CALL DSYTRD( UPLO, N, A, LDA, WORK( INDD ), WORK( INDE ),
$             WORK( INDTAU ), WORK( INDWK ), LLWORK, IINFO )
*
*     If all eigenvalues are desired
*     then call DSTERF or DSTEMR and DORMTR.
*
IF( ( ALLEIG .OR. ( INDEIG .AND. IL.EQ.1 .AND. IU.EQ.N ) ) .AND.
$    IEEEOK.EQ.1 ) THEN
IF( .NOT.WANTZ ) THEN
CALL DCOPY( N, WORK( INDD ), 1, W, 1 )
CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
CALL DSTERF( N, W, WORK( INDEE ), INFO )
ELSE
CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
CALL DCOPY( N, WORK( INDD ), 1, WORK( INDDD ), 1 )
*
IF (ABSTOL .LE. TWO*N*EPS) THEN
TRYRAC = .TRUE.
ELSE
TRYRAC = .FALSE.
END IF
CALL DSTEMR( JOBZ, 'A', N, WORK( INDDD ), WORK( INDEE ),
$                   VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ,
$                   TRYRAC, WORK( INDWK ), LWORK, IWORK, LIWORK,
$                   INFO )
*
*
*
*        Apply orthogonal matrix used in reduction to tridiagonal
*        form to eigenvectors returned by DSTEIN.
*
IF( WANTZ .AND. INFO.EQ.0 ) THEN
INDWKN = INDE
LLWRKN = LWORK - INDWKN + 1
CALL DORMTR( 'L', UPLO, 'N', N, M, A, LDA,
$                      WORK( INDTAU ), Z, LDZ, WORK( INDWKN ),
$                      LLWRKN, IINFO )
END IF
END IF
*
*
IF( INFO.EQ.0 ) THEN
*           Everything worked.  Skip DSTEBZ/DSTEIN.  IWORK(:) are
*           undefined.
M = N
GO TO 30
END IF
INFO = 0
END IF
*
*     Otherwise, call DSTEBZ and, if eigenvectors are desired, DSTEIN.
*     Also call DSTEBZ and DSTEIN if DSTEMR fails.
*
IF( WANTZ ) THEN
ORDER = 'B'
ELSE
ORDER = 'E'
END IF

CALL DSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL,
$             WORK( INDD ), WORK( INDE ), M, NSPLIT, W,
$             IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWK ),
$             IWORK( INDIWO ), INFO )
*
IF( WANTZ ) THEN
CALL DSTEIN( N, WORK( INDD ), WORK( INDE ), M, W,
$                IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ,
$                WORK( INDWK ), IWORK( INDIWO ), IWORK( INDIFL ),
$                INFO )
*
*        Apply orthogonal matrix used in reduction to tridiagonal
*        form to eigenvectors returned by DSTEIN.
*
INDWKN = INDE
LLWRKN = LWORK - INDWKN + 1
CALL DORMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z,
$                LDZ, WORK( INDWKN ), LLWRKN, IINFO )
END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
*  Jump here if DSTEMR/DSTEIN succeeded.
30 CONTINUE
IF( ISCALE.EQ.1 ) THEN
IF( INFO.EQ.0 ) THEN
IMAX = M
ELSE
IMAX = INFO - 1
END IF
CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
END IF
*
*     If eigenvalues are not in order, then sort them, along with
*     eigenvectors.  Note: We do not sort the IFAIL portion of IWORK.
*     It may not be initialized (if DSTEMR/DSTEIN succeeded), and we do
*     not return this detailed information to the user.
*
IF( WANTZ ) THEN
DO 50 J = 1, M - 1
I = 0
TMP1 = W( J )
DO 40 JJ = J + 1, M
IF( W( JJ ).LT.TMP1 ) THEN
I = JJ
TMP1 = W( JJ )
END IF
40       CONTINUE
*
IF( I.NE.0 ) THEN
W( I ) = W( J )
W( J ) = TMP1
CALL DSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
END IF
50    CONTINUE
END IF
*
*     Set WORK(1) to optimal workspace size.
*
WORK( 1 ) = LWKOPT
IWORK( 1 ) = LIWMIN
*
RETURN
*
*     End of DSYEVR
*
END

