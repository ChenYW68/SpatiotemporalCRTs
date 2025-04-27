#include <cmath>
#include <Rcpp.h>
#include <string>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include <limits>


#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

void zeros(double *a, int n){
  for(int i = 0; i < n; i++)
    a[i] = 0.0;
}


// [[Rcpp::export]]
SEXP Empical_Cov_taper_c(SEXP x_r, SEXP y_r,
                         SEXP mu_x_r, SEXP mu_y_r,
                         SEXP Index_r,
                         SEXP n_r,
                         SEXP m_r,
                         SEXP E_r,
                         SEXP l_r,
                         SEXP nThreads_r){

  const int inc = 1;
  const double one = 1.0;
  const double zero = 0.0;
  char const *lower = "L";
  char const *ntran = "N";
  char Trans = 'T';

  int n = INTEGER(n_r)[0];
  int m = INTEGER(m_r)[0];
  int E = INTEGER(E_r)[0];
  int l = INTEGER(l_r)[0];
  int s, i, j0, j1, e, info, p = 1;

  int nThreads = INTEGER(nThreads_r)[0];
  int threadID = 0;

  int *Index = INTEGER(Index_r);
  double *x = REAL(x_r);
  double *y = REAL(y_r);

  double *mu_x = REAL(mu_x_r);
  double *mu_y = REAL(mu_y_r);


  int nProtect=0;
  SEXP R_r;
  PROTECT(R_r = Rf_allocMatrix(REALSXP, l, 1)); nProtect++; //Rf_allocVector
  //double *alpha = REAL(alpha_r);


  // double *R = (double *) R_alloc(n*m*nThreads, sizeof(double));
  double  temp_x = 0.0;

#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif

  //private(i, j, Dist, K, X, KerX, XtX, alpha_temp, tmp_p)

#ifdef _OPENMP
#pragma omp parallel for private(i, j0, j1, temp_x, threadID)
#endif
  for(i = 0; i < l; i++){
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif
    j0 = Index[i*2] - 1;
    j1 = Index[i*2 + 1] - 1;
    temp_x = 0.0;
    for(e = 0; e < E; e++){
      temp_x = temp_x + (x[e*n + j0] - mu_x[j0])*(y[e*m + j1] - mu_y[j1]);
    }
    temp_x = temp_x/(E - 1);
    F77_NAME(dcopy)(&p, &temp_x, &inc, &REAL(R_r)[i], &inc);


  }

  //make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 1;
  PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

  SET_VECTOR_ELT(result_r, 0, R_r);
  SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("R"));

  Rf_namesgets(result_r, resultName_r);

  //unprotect
  UNPROTECT(nProtect);
  // SEXP out;
  // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
  // B = REAL(out);

  return(result_r);

}





// [[Rcpp::export]]
SEXP Empical_Cov_c(SEXP x_r, SEXP y_r,
                   SEXP mu_x_r, SEXP mu_y_r,
                   SEXP n_r,
                   SEXP m_r,
                   SEXP E_r,
                   SEXP nThreads_r){

  const int inc = 1;
  const double one = 1.0;
  const double zero = 0.0;
  char const *lower = "L";
  char const *ntran = "N";
  char Trans = 'T';

  int n = INTEGER(n_r)[0];
  int m = INTEGER(m_r)[0];
  int E = INTEGER(E_r)[0];
  int s, i, j, e, info, p = 1;

  int nThreads = INTEGER(nThreads_r)[0];
  int threadID = 0;


  double *x = REAL(x_r);
  double *y = REAL(y_r);

  double *mu_x = REAL(mu_x_r);
  double *mu_y = REAL(mu_y_r);


  int nProtect=0;
  SEXP R_r;
  PROTECT(R_r = Rf_allocMatrix(REALSXP, n, m)); nProtect++; //Rf_allocVector
  //double *alpha = REAL(alpha_r);


  double *R = (double *) R_alloc(n*m*nThreads, sizeof(double));
  double  temp_x = 0.0;

#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif

  //private(i, j, Dist, K, X, KerX, XtX, alpha_temp, tmp_p)

#ifdef _OPENMP
#pragma omp parallel for private(i, j, temp_x, threadID)
#endif
  for(i = 0; i < n; i++){
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif

#ifdef _OPENMP
#pragma omp parallel for private(j, temp_x, threadID)
#endif
    for(j = 0; j < m; j++){
#ifdef _OPENMP
      threadID = omp_get_thread_num();
#endif
      // My[m*threadID] <- M[i*n + j]*y[j];
      // }
      //  F77_NAME(dgemv)(&Trans, &p, &n, &one, &M[i*n], &p, y, &inc, &zero, &R[threadID], &inc);
      temp_x = 0.0;
      for(e = 0; e < E; e++){
        temp_x = temp_x + (x[e*n + i] - mu_x[i])*(y[e*m + j] - mu_y[j]);
      }
      temp_x = temp_x/(E - 1);
      F77_NAME(dcopy)(&p, &temp_x, &inc, &REAL(R_r)[j*n + i], &inc);
    }

  }

  //make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 1;
  PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

  SET_VECTOR_ELT(result_r, 0, R_r);
  SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("R"));

  Rf_namesgets(result_r, resultName_r);

  //unprotect
  UNPROTECT(nProtect);
  // SEXP out;
  // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
  // B = REAL(out);

  return(result_r);

}










