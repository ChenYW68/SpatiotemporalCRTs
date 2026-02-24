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




// // [[Rcpp::export]]
// SEXP Muti_y(SEXP y_r, SEXP M_r,
//             SEXP n_r,
//             SEXP m_r,
//             SEXP E_r,
//             SEXP nThreads_r){
//
//   const int inc = 1;
//   const double one = 1.0;
//   const double zero = 0.0;
//   char const *lower = "L";
//   char const *ntran = "N";
//   char Trans = 'T';
//
//   int n = INTEGER(n_r)[0];
//   int m = INTEGER(m_r)[0];
//   int E = INTEGER(E_r)[0];
//   int s, i, e, info, p = 1;
//
//   int nThreads = INTEGER(nThreads_r)[0];
//   int threadID = 0;
//
//
//   double *y = REAL(y_r);
//   double *M = REAL(M_r);
//
//
//
//   int nProtect=0;
//   SEXP R_r;
//   PROTECT(R_r = Rf_allocMatrix(REALSXP, m, E)); nProtect++; //Rf_allocVector
//   //double *alpha = REAL(alpha_r);
//
//
//   double *R = (double *) R_alloc(m*E*nThreads, sizeof(double));
//
//
// #ifdef _OPENMP
//   omp_set_num_threads(nThreads);
// #else
//   if(nThreads > 1){
//     warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
//     nThreads = 1;
//   }
// #endif
//
//   //private(i, j, Dist, K, X, KerX, XtX, alpha_temp, tmp_p)
//
// #ifdef _OPENMP
// #pragma omp parallel for private(i, threadID)
// #endif
//   for(e = 0; e < E; e++){
// #ifdef _OPENMP
//     threadID = omp_get_thread_num();
// #endif
//
// #ifdef _OPENMP
// #pragma omp parallel for private(threadID)
// #endif
//     for(i = 0; i < m; i++){
// #ifdef _OPENMP
//       threadID = omp_get_thread_num();
// #endif
//       // My[m*threadID] <- M[i*n + j]*y[j];
//       // }
//       //  F77_NAME(dgemv)(&Trans, &p, &n, &one, &M[i*n], &p, y, &inc, &zero, &R[threadID], &inc);
//
//       F77_NAME(dgemv)(ntran, &p, &n, &one, &M[i*n], &p, &y[e*n], &inc, &zero, &R[E*m*threadID + e*m + i], &inc);
//
//
//
//       //F77_NAME(dgemv)(&Trans, &n, &p, &one, &KerX[p*n*threadID], &n, y, &inc, &zero, &tmp_p[p*threadID], &inc);
//
//
//
//       F77_NAME(dcopy)(&p, &R[E*m*threadID + e*m + i], &inc, &REAL(R_r)[e*m + i], &inc);
//     }
//
//   }
//
//   //make return object
//   SEXP result_r, resultName_r;
//   int nResultListObjs = 1;
//   PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
//   PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
//
//   SET_VECTOR_ELT(result_r, 0, R_r);
//   SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("R"));
//
//   Rf_namesgets(result_r, resultName_r);
//
//   //unprotect
//   UNPROTECT(nProtect);
//   // SEXP out;
//   // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
//   // B = REAL(out);
//
//   return(result_r);
//
// }



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










// // [[Rcpp::export]]
// SEXP theta2_fun_C(SEXP m_r, SEXP n_r,
//                   SEXP theta2_r, SEXP ds_r,
//                   SEXP D_r, SEXP bandKernel_r,
//                   SEXP Q_r, SEXP S00_r, SEXP S01_r,
//                   SEXP theta1_r, SEXP sigTheta1_r,
//                   SEXP nThreads_r){
//   int m = INTEGER(m_r)[0];
//   int n = INTEGER(n_r)[0];
//   double *theta2 = REAL(theta2_r);
//   double ds = REAL(ds_r)[0];
//   double *D = REAL(D_r);
//   double *bandKernel = REAL(bandKernel_r);
//   double *Q = REAL(Q_r);
//   double *S00 = REAL(S00_r);
//   double *S01 = REAL(S01_r);
//   double theta1 = REAL(theta1_r)[0];
//   double sigTheta1 = REAL(sigTheta1_r)[0];
//   int nThreads = INTEGER(nThreads_r)[0];
//
//   int threadID = 0;
//
//
// #ifdef _OPENMP
//   omp_set_num_threads(nThreads);
// #else
//   if(nThreads > 1){
//     warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
//     nThreads = 1;
//   }
// #endif
//
//   const int inc = 1;
//   int i, j, k, l, nProtect=0;
//   double *M = (double *) R_alloc(m*m*nThreads, sizeof(double)); //zeros(M, m*m*nThreads);
//   double *Qt = (double *) R_alloc(m*m*nThreads, sizeof(double)); //zeros(Qt, m*m*nThreads);
//   double *Qm = (double *) R_alloc(m*m*nThreads, sizeof(double)); //zeros(Qm, m*m*nThreads);
//   double *M_Q_M = (double *) R_alloc(m*m*nThreads, sizeof(double)); //zeros(M_Q_M, m*m*nThreads);
//   double *Q_S0 = (double *) R_alloc(m*m*nThreads, sizeof(double)); //zeros(Q_S0, m*m*nThreads);
//   double *Q_S1 = (double *) R_alloc(m*m*nThreads, sizeof(double)); //zeros(Q_S1, m*m*nThreads);
//
//   SEXP logLik_r;
//   PROTECT(logLik_r = Rf_allocVector(REALSXP, n)); nProtect++;
//
//   //double *f = (double *) R_alloc(n*nThreads, sizeof(double)); zeros(f, n*nThreads);
//   char const *ntran = "N";
//   const double one = 1.0;
//   const double zero = 0.0;
//   char Trans = 'T';
//   double f1, f2;
//   int mm = m*m;
//   double f;
//
// #ifdef _OPENMP
// #pragma omp parallel for private(i, j, l, f, f1, f2, threadID)
// #endif
//   for(k = 0; k < n; k++){
// #ifdef _OPENMP
//     threadID = omp_get_thread_num();
// #endif
//     for(i = 0; i < m; i++){
//       for(j = i; j < m; j++){
//         M[m*m*threadID + i*m + j] = ds*exp(-pow(D[i*m + j], 2)/theta2[k])*bandKernel[i*m + j];
//         M[m*m*threadID + j*m + i] = M[m*m*threadID + i*m + j];
//       }
//     }
//     //F77_NAME(dcopy)(&mm, Q, &inc, &Qt[m*m*threadID], &inc);
//     //zeros(&Qm[m*m*threadID], m*m);
//
//     F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, Q, &m, &M[m*m*threadID], &m, &zero, &Qm[m*m*threadID], &m);
//     //Rprintf("f1 =  %3.10f; f2 =  %3.10f;  \n", Q_M[m*m*threadID], Q_M[m*m*threadID + 15]);
//     F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, &Qm[m*m*threadID], &m, S01, &m, &zero, &Q_S0[m*m*threadID], &m);
//     // zeros(M_Q_M, m*m*nThreads);
//     F77_NAME(dgemm)(&Trans, ntran, &m, &m, &m, &one, &M[m*m*threadID], &m, &Qm[m*m*threadID], &m, &zero, &M_Q_M[m*m*threadID], &m);
//     //zeros(Q_S1, m*m*nThreads);
//
//     F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, &M_Q_M[m*m*threadID], &m, S00, &m, &zero, &Q_S1[m*m*threadID], &m);
//
//     f1 = 0.0;
//     f2 = 0.0;
//     for(l = 0; l < m; l++){
//       // for(j = i; j < i + 1; j++){
//       f1 += theta1*Q_S0[m*m*threadID + l*m + l];
//       f2 -= (theta1 + pow(sigTheta1, 2))*Q_S1[m*m*threadID + l*m + l]/2;
//       //}
//     }
//     //Rprintf("f1 =  %3.5f; f2 =  %3.5f;  \n", f1, f2);
//     //f[k*threadID] = f1 + f2;
//     f =  f1 + f2;
//     F77_NAME(dcopy)(&inc, &f, &inc, &REAL(logLik_r)[k], &inc);
//
//   }
//
//
//   SEXP result_r, resultName_r;
//   int nResultListObjs = 1;
//
//   PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
//   PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
//   SET_VECTOR_ELT(result_r, 0, logLik_r);
//   SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("logLik"));
//   Rf_namesgets(result_r, resultName_r);
//
//   //unprotect
//   UNPROTECT(nProtect);
//
//   return(result_r);
// }


















