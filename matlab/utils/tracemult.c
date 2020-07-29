#include "mex.h"
#include "lapack.h"
#include "blas.h"

/*
 * Computes "trace" matrix product and indexing operations
 * C(:,:,i) = A(:,:,i)*B(:,:,j(i)) or
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
void matlab_parse3(const mxArray *prhs, mwIndex *dims, double **A, bool iscomplex)
{
    mwIndex dimcount, *rhsdims, i, itwo;
    
    *A = mxGetPr(prhs);
    dimcount=mxGetNumberOfDimensions(prhs);
    rhsdims = (mwIndex *)mxGetDimensions(prhs);
    for (i=0; i<dimcount; i++) dims[i]=rhsdims[i];
    for (i=dimcount; i<3; i++) dims[i]=1;
    if (mxIsComplex(prhs) | iscomplex) {
        /* Pretend 2xDouble == Complex to pass this mess to BLAS */
        iscomplex = true;
        i = dims[0]*dims[1]*dims[2];
        *A = (double *)malloc(sizeof(double)*2*i);
        itwo = 2;
        dcopy(&i, mxGetPr(prhs), &ione, *A, &itwo);
        if (mxIsComplex(prhs)) {
            dcopy(&i, mxGetPi(prhs), &ione, &((*A)[1]), &itwo);
        }
        else {
            for (i=0; i<dims[0]*dims[1]*dims[2]; i++) (*A)[2*i+1] = 0.0;
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/*  RHS: A(n x s), j(n x 1), or
 *       A(p x m x n), j(n x 1),  B(m x k x s)
 *  LHS: C(p x k x n), no p and k if no B
 */
{
  double *A, *B, *C, *dj;
  mwIndex n, m, k, p, dims[3], i, j;
  bool iscomplex;
  double zone[2], zzero[2]; /* Fortran housekeeping */
          
  zone[0] = 1.0; zone[1] = 0.0;   zzero[0] = 0.0; zzero[1] = 0.0;

  if (nrhs<2) { mexPrintf("Specify at least A and j\n"); return; }
  
  iscomplex = mxIsComplex(prhs[0]);
  if (nrhs>2) iscomplex = iscomplex | mxIsComplex(prhs[2]);

  /* Fetch the data. We need to determine the sizes correctly */
  /* A */
  matlab_parse3(prhs[0], dims, &A, iscomplex);
  if (nrhs>2) {
      /* We have a 3D array for mult */
      p = dims[0];
      m = dims[1];
      n = dims[2];
  }
  else {
      /* Just 2D array for trace */
      n = dims[0];
      k = 1;
      p = 0;
  }
  /* double j */
  matlab_parse3(prhs[1], dims, &dj, false);
  if (dims[0]*dims[1]!=n) { mexPrintf("size(j,1) differs from size(A)\n"); return; }
  /* B if present */
  if (nrhs>2) {
      matlab_parse3(prhs[2], dims, &B, iscomplex);
      if (dims[0]!=m) { mexPrintf("size(A,2) differs from size(B,1)\n"); return; }
      k = dims[1];
  }
  
  /* Create the output */
  if (p>0) {
      dims[0] = p;
      dims[1] = k;
      dims[2] = n;
      i = 3;
      if (iscomplex) plhs[0] = mxCreateNumericArray(i, dims, mxDOUBLE_CLASS, mxCOMPLEX);
      else           plhs[0] = mxCreateNumericArray(i, dims, mxDOUBLE_CLASS, mxREAL);
  }
  else {
      if (iscomplex) plhs[0] = mxCreateDoubleMatrix(n, k, mxCOMPLEX);
      else           plhs[0] = mxCreateDoubleMatrix(n, k, mxREAL);
  }
  
  matlab_parse3(plhs[0], dims, &C, iscomplex);

  if (p>0) {
      /* tracemult op */
      if (iscomplex) {
          for (i=0; i<n; i++) {
              j = (mwIndex)dj[i];
              j--; /* matlab->C indexing */
              zgemm(&cN,&cN, &p,&k,&m, zone,&A[2*i*p*m],&p, &B[2*j*m*k],&m, zzero,&C[2*i*p*k],&p);
          }
          /* Copy stuff back */
          i = dims[0]*dims[1]*dims[2];
          j = 2;
          dcopy(&i, C, &j,     mxGetPr(plhs[0]), &ione);
          dcopy(&i, &C[1], &j, mxGetPi(plhs[0]), &ione);
          free(A);
          free(B);
          free(C);
      }
      else {
          for (i=0; i<n; i++) {
              j = (mwIndex)dj[i];
              j--; /* matlab->C indexing */
              dgemm(&cN,&cN, &p,&k,&m, &done,&A[i*p*m],&p, &B[j*m*k],&m,&dzero,&C[i*p*k],&p);
          }
      }
  }
  else {
      /* Just sample here */
      if (iscomplex) {
          for (i=0; i<n; i++) {
              j = (mwIndex)dj[i];
              j--; /* matlab->C indexing */
              C[2*i] = A[2*(i + j*n)];
              C[2*i+1] = A[2*(i + j*n)+1];
          }
          /* Copy stuff back */
          i = dims[0]*dims[1]*dims[2];
          j = 2;
          dcopy(&i, C, &j,     mxGetPr(plhs[0]), &ione);
          dcopy(&i, &C[1], &j, mxGetPi(plhs[0]), &ione);
          free(A);
          free(C);
      }
      else {
          for (i=0; i<n; i++) {
              j = (mwIndex)dj[i];
              j--; /* matlab->C indexing */
              C[i] = A[i + j*n];
          }
      }
  }
}

