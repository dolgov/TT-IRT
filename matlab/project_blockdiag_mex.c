#include "mex.h"
#include "lapack.h"
#include "blas.h"
#include <string.h>

/*
 * Projects a block-diagonal linear system, with a matrix given in a 3D TT format, with the first block full, onto a solution block (accumulates "Phi_left")
 */

/* Fortran housekeeping */
mwIndex ione = 1;
char cN = 'N';
char cT = 'T';
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
/* RHS: XAX1 [rx1',rx1,rc1], Ci [rc1,n,rc2], X [rx1,n,rx2]
   LHS: XAX1' [rx2',rx2,rc2]
*/
{
  /* Todo: Make it adaptive to complex arithmetics! */
  double *XAX1, *Ci, *x, *XAX1_new, *xperm, *Cperm;
  double *dtmp1, *dtmp2;
/*   double *XXAX1, *XXAX1C, *XXAX1CT, *XAX_new_i_1, *XAX_new_i_2; */
  mwIndex rx1, rc1, n, rc2, rx2, dims[3], i,j, tmpsize;

  if (nrhs<3) { mexPrintf("Specify at least XAX1,C,x\n"); return; }

  /* Fetch the data. We need to determine the sizes correctly */
  matlab_size3(prhs[0], dims, &XAX1);
  rx1 = dims[0];
  rc1 = dims[2];
  matlab_size3(prhs[1], dims, &Ci);
  if (dims[0]!=rc1) { mexPrintf("rc1 in XAX1 and C are not consistent!\n"); return; }
  n = dims[1];
  rc2 = dims[2];
  matlab_size3(prhs[2], dims, &x);
  if (dims[0]!=rx1) { mexPrintf("rx1 in XAX1 and x are not consistent!\n"); return; }
  if (dims[1]!=n)   { mexPrintf("n in C and x are not consistent!\n"); return; }
  rx2 = dims[2];

  /* Create the output */
  /*Replace by numericarray*/
  plhs[0] = mxCreateDoubleMatrix(rx2*rx2, rc2, mxREAL);
  XAX1_new = mxGetPr(plhs[0]);
  /* Clear it */
  memset(XAX1_new, 0, sizeof(double)*rx2*rx2*rc2);

  /* Permute x and C for easier blas */
  xperm = (double *)malloc(sizeof(double)*rx1*rx2*n);
  for (i=0; i<n; i++) {
    for (j=0; j<rx2; j++) {
      dcopy(&rx1, &x[i*rx1+j*rx1*n], &ione, &xperm[j*rx1+i*rx1*rx2], &ione);
    }
  }
  Cperm = (double *)malloc(sizeof(double)*rc1*rc2*n);
  for (i=0; i<n; i++) {
    for (j=0; j<rc2; j++) {
      dcopy(&rc1, &Ci[i*rc1+j*rc1*n], &ione, &Cperm[j*rc1+i*rc1*rc2], &ione);
    }
  }

  /* Temp storage */
  /* Determine the maximal size we will need */
  tmpsize = rx2*rx1;
  if (rx2>rx1) tmpsize = rx2*rx2;
  if (rc2>rc1) tmpsize = tmpsize*rc2;
  else         tmpsize = tmpsize*rc1;
  dtmp1 = (double *)malloc(sizeof(double)*tmpsize);
  dtmp2 = (double *)malloc(sizeof(double)*tmpsize);
  /* XXAX1       = (double *)malloc(sizeof(double)*rx2*rx1*rc1);
  XXAX1C      = (double *)malloc(sizeof(double)*rx2*rx1*rc2);
  XXAX1CT     = (double *)malloc(sizeof(double)*rx2*rx1*rc2);
  XAX_new_i_1 = (double *)malloc(sizeof(double)*rx2*rx2*rc2);
  XAX_new_i_2 = (double *)malloc(sizeof(double)*rx2*rx2*rc2); */

  /* Main loop -- reduce over n */
  for (i=0; i<n; i++) {
    tmpsize = rx1*rc1;
    /* tmp1 = X'*XAX1 */
    dgemm(&cT,&cN,&rx2,&tmpsize,&rx1, &done,&xperm[i*rx1*rx2],&rx1, XAX1,&rx1,&dzero,dtmp1,&rx2);
    tmpsize = rx2*rx1;
    /* tmp2 = tmp1*C */
    dgemm(&cN,&cN,&tmpsize,&rc2,&rc1, &done,dtmp1,&tmpsize, &Cperm[i*rc1*rc2],&rc1,&dzero,dtmp2,&tmpsize);
    /* Permute this partial projection */
    /* It was rx2, rx1*rc2, should become rx1*rc2, rx2 */
    tmpsize = rx1*rc2;
    for (j=0; j<rx2; j++) {
      dcopy(&tmpsize, &dtmp2[j], &rx2, &dtmp1[j*tmpsize], &ione);
    }
    /* tmp2 = tmp1*X */
    tmpsize = rc2*rx2;
    dgemm(&cT,&cN,&rx2,&tmpsize,&rx1, &done,&xperm[i*rx1*rx2],&rx1, dtmp1,&rx1,&dzero,dtmp2,&rx2);
    /* Permute tmp2 */
    /* Was rx2*rc2, rx2', needed rx2', rx2*rc2 */
    for (j=0; j<rx2; j++) {
      dcopy(&tmpsize, &dtmp2[j*tmpsize], &ione, &dtmp1[j], &rx2);
    }
    /* Add to XAX1_new */
    tmpsize = rx2*rx2*rc2;
    daxpy(&tmpsize, &done, dtmp1, &ione, XAX1_new, &ione);
  }

  free(dtmp1);
  free(dtmp2);
  free(xperm);
  free(Cperm);
}
