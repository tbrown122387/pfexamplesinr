#include <RcppEigen.h>
#include "svol_sisr_hilb.h"

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;               	// 'maps' rather than copies 
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers

// 1e3 particles, 5 bits 
using svol_pfilter = pf::filters::svol_sisr_hilb<1000,5,pf::resamplers::sys_hilb_resampler<1000,1,5,double>,double> sisrsvol_hilb;
//(phi,beta,sigma);

// [[Rcpp::export]]
VectorXd approxLL(Map<VectorXd> y, Map<VectorXd> thetaProposal, Map<VectorXd> uProposal) {
  svol_pfilter pf(thetaProposal(0), thetaProposal(1) thetaProposal(2)); // order: phi, beta, sigma
  
  return es.eigenvalues();
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
