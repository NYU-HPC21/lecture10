/* Jacobi smoothing to solve -u''=f
 * Global vector has N inner unknowns.
 * Author: Georg Stadler
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "utils.h"

/* compuate global residual, assuming ghost values are updated */
double compute_residual(double *u, int N, double invhsq){
  int i;
  double tmp, res = 0.0;

  for (i = 1; i <= N; i++){
    tmp = ((2.0*u[i] - u[i-1] - u[i+1]) * invhsq - 1);
    res += tmp * tmp;
  }
  return sqrt(res);
}


int main(int argc, char * argv[]){
  int i, N, iter, max_iters;

  sscanf(argv[1], "%d", &N);
  sscanf(argv[2], "%d", &max_iters);

  /* timing */
  Timer tt;
  tt.tic();

  /* Allocation of vectors, including left and right ghost points */
  double * u    = (double *) calloc(sizeof(double), N+2);
  double * unew = (double *) calloc(sizeof(double), N+2);
  double * utmp;

  double h = 1.0 / (N + 1);
  double hsq = h * h;
  double invhsq = 1./hsq;
  double res, res0, tol = 1e-5;

  /* initial residual */
  res0 = compute_residual(u, N, invhsq);
  res = res0;
  u[0] = u[N+1] = 0.0;

  for (iter = 0; iter < max_iters && res/res0 > tol; iter++) {

    /* Jacobi step for all the inner points */
    for (i = 1; i <= N; i++){
      unew[i]  = 0.5 * (hsq + u[i - 1] + u[i + 1]);
    }

    /* flip new_u and u (avoids copy) */
    utmp = u; u = unew; unew = utmp;

    if (0 == (iter % 10)) {
      res = compute_residual(u, N, invhsq);
      printf("Iter %d: Residual: %g\n", iter, res);
    }
  }

  /* Clean up */
  free(u);
  free(unew);

  /* timing */
  double elapsed = tt.toc();
  printf("Time elapsed is %f seconds.\n", elapsed);
  return 0;
}
