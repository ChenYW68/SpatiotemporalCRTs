// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace std;
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
double Trace_Matrix(arma::colvec &h, arma::mat &S)
{
  // trace()h %*% t(h) %*% S)
 // arma::mat C(x.n_rows, x.n_rows, arma::fill::ones);
  double C = arma::trace(h * trans(h)*S);
  // double C = (trans(h)*S * h);
  // double C = arma::trace(h*S);
  return(C);
}
// [[Rcpp::export]]
double Trace_Muti(arma::mat S1, arma::mat S2)
{
  // trace()h %*% t(h) %*% S)
  // arma::mat C(S.n_rows, S.n_rows, arma::fill::zeros);
  // arma::mat PhiB;
  // for(int b = 0; b < S.n_rows; b++)
  // {
  //   PhiB = diagmat(Ms.row(b));
  //   C = C + PhiB * S * PhiB;
  // }
  double C = arma::trace(S1*S2);
  return(C);
}
// [[Rcpp::export]]
arma::mat Solve(arma::mat A)
{
  mat B = arma::inv(A);
  return(B);
}

