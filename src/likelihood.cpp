#include <RcppEigen.h>
#include "svol_sisr_hilb.h"
#include "resamplers.h"

// [[Rcpp::depends(RcppEigen)]]


// #include <Rcpp.h>
// using namespace Rcpp;


using Eigen::Map; 
using Eigen::MatrixXd;
using Eigen::VectorXd;

// 1e3 particles, 5 bits 
#define NP 1000
#define NB 5

// some type aliases
using svol_pfilter = svol_sisr_hilb<NP,NB,pf::resamplers::sys_hilb_resampler<NP,1,NB,double> , double>; //(phi,beta,sigma);

// uProposal will be (time X (NP+1))
// first NP columns will be used for state sampling 
// last column will be used for resampling at each time point

// [[Rcpp::export]]
VectorXd approxLL(const Map<VectorXd> y, const Map<VectorXd> thetaProposal, const Map<MatrixXd> uProposal) {

  // construct particle filter object
  svol_pfilter pf(thetaProposal(0), thetaProposal(1), thetaProposal(2)); // order: phi, beta, sigma

  // iterate over the data 
  double log_like(0.0);
  // Eigen::Matrix<double,1,1> yt;
  // std::array<Eigen::Matrix<double,1,1>, NP> uStateTransition;
  // Eigen::Matrix<double,1,1> uResample;
  for(int time = 0; time < y.rows(); ++time){
    
    // change types of inputs
    // yt(0) = y(time);
    // for(unsigned particle = 0; particle < NP; ++particle) {
    //   uStateTransition[particle] = uProposal.block(time,particle,1,1); 
    // }
    // uResample(0) = uProposal(time,NP);
    
    // update particle filter and log-likelihood
    // pf.filter(yt, uStateTransition, uResample);
    // log_like += pf.getLogCondLike();
    log_like += 1.0;
  }
    
  //return es.eigenvalues();
  VectorXd r;
  r(0) = log_like;
  return r;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.


/*** R
numTime <- 3
numParts <- 1000
u <- matrix(rnorm(numTime*(numParts+1)), ncol = numParts+1)
approxLL(1:numTime, 1:3, u)
*/
