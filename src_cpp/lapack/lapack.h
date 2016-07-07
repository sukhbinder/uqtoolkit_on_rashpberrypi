#ifndef LAPACK_H
#define LAPACK_H

#include "ftndefs.h"


// Computes eigenvalues and eigenvectors of a real symmetric matrix A
// See dsyevx.f for details on the arguments.
extern FTN_FUNC void FTN_NAME(dsyevx)(char*, char*, char*, int*,
				      double*, int*, double*, double*,
				      int*, int*, double*, int*, 
				      double*, double*, int*, double*, int*,
				      int*, int*, int*);




// Computes Cholesky decomposition of a real symmetric matrix A
// See dpotrf.f for details on the arguments.
extern FTN_FUNC void FTN_NAME(dpotrf)(char*, int*, double*, int*, int*);

// Computes the inverse of a real symmetric matrix A
// See dpotri.f for details on the arguments.
extern FTN_FUNC void FTN_NAME(dpotri)(char*, int*, double*, int*, int*);

// Computes the inverse of a complex matrix
// see zgetri.f for details on the arguments.
extern FTN_FUNC void FTN_NAME(zgetri)(int*, double*, int*, int*, double*, int*, int*);


// Computes the least squares solution of an over-determined Ax=b
// see dgels.f for details on the arguments
extern FTN_FUNC void FTN_NAME(dgels)(char*, int*, int*, int*, double*, int*, double*, int*, double*, int*, int*);

// Computes the LU factorization of a general matrix
// see dgetrf.f for details on the arguments
extern FTN_FUNC void FTN_NAME(dgetrf)(int*, int*, double*, int*, int*, int*);

// Computes the inverse of a general matrix using its LU factorization
// see dgetri.f for details on the arguments
extern FTN_FUNC void FTN_NAME(dgetri)(int*, double*, int*, int*, double*, int*, int*);

// Computes the solution to a real system of linear equations AX = B
// see dgesv.f for details on the arguments
extern FTN_FUNC void FTN_NAME(dgesv)(int*, int*, double*, int*, int*, double*, int*, int*);

//
extern FTN_FUNC void FTN_NAME(dstev)( char *, int *, double *, double *, double *, int *, double *, int* ) ;

//
extern FTN_FUNC void FTN_NAME(dsteqr)( char *, int *, double *, double *, double *, int *, double *, int* ) ;


#endif  /* LAPACK_H */
