#include "mex.h"
#include "lapack.h"
#include "blas.h"

extern void tt_irt1(mwIndex d, mwIndex *n, double *xs, mwIndex *ttrank, double *ttcore, mwIndex M, double *q, double *z, double *lPz);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/*  RHS: n(d x 1), xs(sum(n) x 1), ttrank(d+1 x 1), ttcore(sum(r*n*r) x 1), q(M x d)
 *  LHS: Z(M x d), lPz (M x 1)
 */
{
  double *dtmp;
  mwIndex *itmp, d, *n, *ttrank, i, M;

  if (nrhs<5) { mexPrintf("RHS should contain n, xs, ttrank, ttcore and q\n"); return; }
  if (nlhs<2) { mexPrintf("LHS should contain z and Pz\n"); return; }

  itmp = (mwIndex *)mxGetDimensions(prhs[0]);
  d = itmp[0]*itmp[1];

  dtmp = mxGetPr(prhs[0]); /* n */
  n = (mwIndex *)malloc(sizeof(mwIndex)*d);
  for (i=0; i<d; i++) n[i]=(mwIndex)dtmp[i];

  dtmp = mxGetPr(prhs[2]); /* ttrank */
  ttrank = (mwIndex *)malloc(sizeof(mwIndex)*(d+1));
  for (i=0; i<=d; i++) ttrank[i]=(mwIndex)dtmp[i];

  /* number of points */
  itmp = (mwIndex *)mxGetDimensions(prhs[4]);
  M = itmp[0];

  /* Create the output */
  plhs[0] = mxCreateDoubleMatrix(M, d, mxREAL);
  i=1;
  plhs[1] = mxCreateDoubleMatrix(M, i, mxREAL);

  /* produce samples */
  tt_irt1(d, n, mxGetPr(prhs[1]), ttrank, mxGetPr(prhs[3]), M, mxGetPr(prhs[4]), mxGetPr(plhs[0]), mxGetPr(plhs[1]));

  free(n);
  free(ttrank);
}
