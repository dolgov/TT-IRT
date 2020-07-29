#include "mex.h"
#include "lapack.h"
#include "blas.h"


/*
 * Solves a block-diagonal linear system, with a matrix given in a 3D TT format, with the first block full
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
/* RHS: XAX1 [rx1',rx1,rc1], CiXC2 [rc1,n,rx2], rhs [rx1*n*rx2]
   LHS: sol [rx1*n*rx2]
*/
{
  /* Todo: Make it adaptive to complex arithmetics! */
  double *XAX1, *CiXC2, *rhs, *sol, *Ai;
  mwIndex rx1, rc1, n, rx2, dims[3], i, tmpsize, *ipiv, info;

  if (nrhs<3) { mexPrintf("Specify at least XAX1,CXC2, rhs\n"); return; }

  /* Fetch the data. We need to determine the sizes correctly */
  /* XAX1 */
  matlab_size3(prhs[0], dims, &XAX1);
  rx1 = dims[0];
  rc1 = dims[2];
  if (dims[1]!=rx1) { mexPrintf("XAX1 is not square!\n"); return; }
  /* CiXC2 */
  matlab_size3(prhs[1], dims, &CiXC2);
  if (dims[0]!=rc1) { mexPrintf("rc1 in XAX1 and C are not consistent!\n"); return; }
  n = dims[1];
  rx2 = dims[2];
  /* rhs */
  matlab_size3(prhs[2], dims, &rhs);
  if (dims[0]*dims[1]*dims[2]!=rx1*n*rx2) { mexPrintf("RHS size is not consistent!\n"); return; }

  /* Create the output */
  plhs[0] = mxCreateDoubleMatrix(rx1*n, rx2, mxREAL);
  sol = mxGetPr(plhs[0]);

  /* rxr storage */
  Ai = (double *)malloc(sizeof(double)*rx1*rx1);

  /* Main loop -- generate and solve */
  tmpsize = rx1*n*rx2;
  dcopy(&tmpsize, rhs, &ione, sol, &ione);
  ipiv = (mwIndex *)malloc(sizeof(mwIndex)*rx1);
  tmpsize = rx1*rx1;
  for (i=0; i<n*rx2; i++) {
    dgemv(&cN, &tmpsize, &rc1, &done, XAX1, &tmpsize, &CiXC2[i*rc1], &ione, &dzero, Ai, &ione);
    dgesv(&rx1, &ione, Ai, &rx1, ipiv, &sol[i*rx1], &rx1, &info);
  }

  free(ipiv);
  free(Ai);
}
