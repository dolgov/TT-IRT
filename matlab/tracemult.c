#include "mex.h"
#include "lapack.h"
#include "blas.h"

/*
 * Computes "trace" matrix product and indexing operations
 * C(i,:) = A(i,:)*B(:,:,j(i)) and
 * C(i) = A(i,j(i))
 *
 *  Compilation: mex -v -largeArrayDims -O -lmwblas -lmwlapack tracemult.c
 */

/* Fortran housekeeping */
mwIndex ione = 1;
char cN = 'N';
double done = 1.0;
double dzero = 0.0;

/* Correctly determine three array dimensions */
void matlab_size3(const mxArray *prhs, mwIndex *dims, double **A)
{
  mwIndex dimcount, *rhsdims, i;

  *A = mxGetPr(prhs);
  dimcount=mxGetNumberOfDimensions(prhs);
  rhsdims = (mwIndex *)mxGetDimensions(prhs);
  for (i=0; i<dimcount; i++) dims[i]=rhsdims[i];
  for (i=dimcount; i<3; i++) dims[i]=1;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/*  RHS: A(n x m), j(n x 1), [ B(m x k x s) ]
 *  LHS: C(n x k), k==1 if no B
 */
{
  double *A, *B, *C, *dj;
  mwIndex n, m, k, s, dims[3], i, j;

  if (nrhs<2) { mexPrintf("Specify at least A and j\n"); return; }

  /* Fetch the data. We need to determine the sizes correctly */
  /* A */
  matlab_size3(prhs[0], dims, &A);
  n = dims[0];
  m = dims[1];
  /* double j */
  matlab_size3(prhs[1], dims, &dj);
  if (dims[0]!=n) { mexPrintf("size(j,1) differs from size(A,1)\n"); return; }
  k = 1;
  s = 0;
  /* B if present */
  if (nrhs>2) {
    matlab_size3(prhs[2], dims, &B);
    if (dims[0]!=m) { mexPrintf("size(A,2) differs from size(B,1)\n"); return; }
    k = dims[1];
    s = dims[2];
  }

  /* Create the output */
  plhs[0] = mxCreateDoubleMatrix(n, k, mxREAL);
  C = mxGetPr(plhs[0]);

  if (s>0) {
    for (i=0; i<n; i++) {
      j = (mwIndex)dj[i];
      j--; /* matlab->C indexing */
      dgemm(&cN,&cN, &ione,&k,&m, &done,&A[i],&n, &B[j*m*k],&m,&dzero,&C[i],&n);
    }
  }
  else {
    /* Just sample here */
    for (i=0; i<n; i++) {
      j = (mwIndex)dj[i];
      j--; /* matlab->C indexing */
      C[i] = A[i + j*n];
    }
  }
}

