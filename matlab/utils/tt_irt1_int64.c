/*
 Compute the inverse Rosenblatt (cumulative) transform via linearly interpolated TT
 Args (sizes):
   lapackint d: dimension
   lapackint *n: mode sizes (d)
   double *xs: array of grid points (sum(n))
   lapackint *ttrank: tt ranks (d+1)
   double *ttcore: TT2.0 contiguous storage of tt cores (sum(r*n*r))
   lapackint M: number of points
   double *q: seed points (M*d)
   double *z: transformed points (M*d) !should be allocated
   double *lPz: log(sampling density) values at z (M) !should be allocated
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// this may change depending on the caller and blas/lapack distribution
#define lapackint long long // for MATLAB

// fortran
double done = 1.0;
double dzero = 0.0;
lapackint ione = 1;
char cN = 'N';

extern void dgemm_(char *, char *, lapackint *, lapackint *, lapackint *, double *, double *, lapackint *, double *, lapackint *, double *, double *, lapackint *);
extern void dscal_(lapackint *, double *, double *, lapackint *);
extern void daxpy_(lapackint *, double *, double *, lapackint *, double *, lapackint *);
extern void dcopy_(lapackint *, double *, lapackint *, double *, lapackint *);


void tt_irt1(lapackint d, lapackint *n, double *xs, lapackint *ttrank, double *ttcore, lapackint M, double *q, double *z, double *lPz)
{
  lapackint k,i,j,m,  rmax, nmax, i0, i1, i2, Mb;
  lapackint *pstt, *psxs; // position vector in tt and xs, resp.
  double **C, C1, C2, *h, *fk, *fk2, *fkm1, *fkm2, *fkm1i;
  double hq, Aq, Bq, Dq, x1, x2, xk;

  rmax = ttrank[0];
  nmax = n[0];
  psxs = (lapackint *)malloc(sizeof(lapackint)*(d+1));
  psxs[0]=0;
  for (k=1; k<=d; k++) {
    psxs[k]=psxs[k-1]+n[k-1];
    if ((k<d)&&(n[k]>nmax)) nmax = n[k];
    if (ttrank[k]>rmax) rmax = ttrank[k];
  }

  Mb = 64; // blocking is much faster
  if (rmax<Mb) rmax=Mb;
  if (nmax+1<Mb) nmax=Mb-1;

  pstt = (lapackint *)malloc(sizeof(lapackint)*(d+1));
  pstt[0]=0;
  for (k=1; k<=d; k++) pstt[k]=pstt[k-1]+ttrank[k-1]*n[k-1]*ttrank[k];

  C = (double **)malloc(sizeof(double *)*d); // integrated factors
  C[d-1] = (double *)malloc(sizeof(double)*1);
  C[d-1][0] = 1.0;
  h = (double *)malloc(sizeof(double)*nmax);
  if (rmax>nmax+1) fk = (double *)malloc(sizeof(double)*rmax*rmax);
  else fk = (double *)malloc(sizeof(double)*rmax*nmax);
  fk2 = (double *)malloc(sizeof(double)*Mb*nmax);

  // backward sweep, integration
  for (k=d-1; k>0; k--) {
    C[k-1] = (double *)malloc(sizeof(double)*ttrank[k]);
    i = ttrank[k]*n[k]; // size for BLAS
    // multiply fk = core[k]*C[k+1]
    dgemm_(&cN, &cN, &i, &ione, &ttrank[k+1], &done, &ttcore[pstt[k]], &i, C[k], &ttrank[k+1], &dzero, fk, &i);
    // grid steps
    for (i=0; i<n[k]-1; i++) h[i] = xs[psxs[k]+i+1] - xs[psxs[k]+i];
    // integrate fk
    memset(C[k-1], 0, sizeof(double)*ttrank[k]);
    for (j=0; j<n[k]-1; j++) {
      for (i=0; i<ttrank[k]; i++) {
        C[k-1][i] += (fk[i+j*ttrank[k]] + fk[i+(j+1)*ttrank[k]])*h[j]*0.5;
      }
    }
  }

  // sampling
  fkm1 = (double *)malloc(sizeof(double)*Mb*rmax);
  fkm2 = (double *)malloc(sizeof(double)*rmax*rmax);
  fkm1i = (double *)malloc(sizeof(double)*rmax);
  for (m=0; m<M; m+=Mb) {
    if (m+Mb>M) Mb = M-m;
    for (i=0; i<Mb; i++) {
      fkm1[i] = 1.0;
      lPz[m+i] = 0.0;
    }

    for (k=0; k<d; k++) {
      i = ttrank[k]*n[k];
      // multiply fk = core[k]*C[k+1]
      dgemm_(&cN, &cN, &i, &ione, &ttrank[k+1], &done, &ttcore[pstt[k]], &i, C[k], &ttrank[k+1], &dzero, fk, &i);

      // grid steps
      for (i=0; i<n[k]-1; i++) h[i] = xs[psxs[k]+i+1] - xs[psxs[k]+i];
      // sample from the left
      dgemm_(&cN, &cN, &Mb, &n[k], &ttrank[k], &done, fkm1, &Mb, fk, &ttrank[k], &dzero, fk2, &Mb);
      // make this nonnegative
      for (j=0; j<n[k]*Mb; j++) fk2[j]=fabs(fk2[j]);
      // integrate up to whole points
      memset(fk, 0, sizeof(double)*Mb);
      for (j=1; j<n[k]; j++) {
        dcopy_(&Mb, &fk[(j-1)*Mb], &ione, &fk[j*Mb], &ione);
        hq = h[j-1]*0.5;
        daxpy_(&Mb, &hq, &fk2[(j-1)*Mb], &ione, &fk[j*Mb], &ione);
        daxpy_(&Mb, &hq, &fk2[j*Mb], &ione, &fk[j*Mb], &ione);
      }
      // Now fk is the same as Ck in matlab
      // check it for zero
      for (i=0; i<Mb; i++) {
        if (fk[i+(n[k]-1)*Mb]==0.0) {
          for (j=0; j<n[k]; j++) {
            fk2[i+j*Mb] = 1.0;
            fk[i+j*Mb] = (double)j;
          }
          hq = 1.0/(n[k]-1);
          dscal_(&n[k], &hq, &fk[i], &Mb);
          dscal_(&n[k], &hq, &fk2[i], &Mb);
        }
        // normalize
        hq = 1.0/fk[i+(n[k]-1)*Mb];
        dscal_(&n[k], &hq, &fk[i], &Mb);
        dscal_(&n[k], &hq, &fk2[i], &Mb);
      }

      /****** Inversion *******/
      // Binary search for the closest to q points
      for (i=0; i<Mb; i++) {
        hq = q[m+i+M*k];
        i0 = 0;
        i2 = n[k]-1;
        while (i2-i0>1) {
          i1 = (int)((i0+i2)*0.5);
          if (hq>fk[i+i1*Mb]) i0 = i1; // go right
          else i2 = i1; // go left
        }

        // now we have i0<=i(qk)<=i0+1
        // Solve quadratic equation
        x1 = xs[psxs[k]+i0];
        x2 = xs[psxs[k]+i0+1];
        C1 = fk2[i+i0*Mb];
        C2 = fk2[i+(i0+1)*Mb];
        hq = x2-x1;
        Aq = 0.5*(C2-C1)/hq;
        Bq = (C1*x2 - C2*x1)/hq;
        Dq = 2.0*Aq*x1+Bq;
        Dq *= Dq;
        Dq += 4.0*Aq*(q[m+i+M*k]-fk[i+i0*Mb]);
        xk = 0.5*(-Bq + sqrt(fabs(Dq)))/Aq;
        // if we have exactly linear equation
        if (Aq==0.0) xk = x1 + (q[m+i+M*k]-fk[i+i0*Mb])/Bq;
        z[m+i+M*k] = xk;

        // interpolation weights
        C1 = (x2-xk)/hq;
        C2 = (xk-x1)/hq;
        // interpolate our actual conditional probability
        lPz[m+i] += log(fabs(fk2[i+i0*Mb]*C1 + fk2[i+(i0+1)*Mb]*C2));

        if (k<d-1) {
          // sample the k-th block onto xk
          for (j=0; j<ttrank[k+1]; j++) {
              dcopy_(&ttrank[k], &ttcore[pstt[k]+i0*ttrank[k]+j*ttrank[k]*n[k]], &ione, &fkm2[j*ttrank[k]], &ione);
              dscal_(&ttrank[k], &C1, &fkm2[j*ttrank[k]], &ione);
              daxpy_(&ttrank[k], &C2, &ttcore[pstt[k]+(i0+1)*ttrank[k]+j*ttrank[k]*n[k]], &ione, &fkm2[j*ttrank[k]], &ione);
          }
          dcopy_(&ttrank[k], &fkm1[i], &Mb, fkm1i, &ione);
          dgemm_(&cN, &cN, &ione, &ttrank[k+1], &ttrank[k], &done, fkm1i, &ione, fkm2, &ttrank[k], &dzero, &fkm1[i], &Mb);

        } //k<d-1
      } // i (Mb)

    } // k
  } // M

  free(fkm1);
  free(fkm2);
  free(fkm1i);
  free(h);
  free(fk);
  free(fk2);
  free(psxs);
  free(pstt);
  for (k=0; k<d; k++) free(C[k]);
  free(C);
}
