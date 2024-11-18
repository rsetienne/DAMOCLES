//' @useDynLib DDD


#define STRICT_R_HEADERS
#include "config.h"
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include "odeint_helper.h"


using namespace Rcpp;


class ode_rhs
{
public:
  ode_rhs(NumericVector M)
  {
    const size_t lv = M.size;
    for (size_t i = 0; i < lv; ++i) {
      lavec[i] = parsvec[i];            // parsvec[1:lv]
	    muvec[i] = parsvec[lv + i];       // parsvec[(lv + 1):(2 * lv)]
	    nn[i] = parsvec[2 * lv + i];      // parsvec[(2 * lv + 1):(3 * lv)]
    }
    kk = static_cast<size_t>(parsvec[parsvec.size() - 1]);
  }
  
  void operator()(const std::vector<double>& xx, std::vector<double>& dx, double /* t */)
  {
    //  dp = M %*% p
    //return(list(dp))
  
    dx.front() = dx.back() = 0.0;
    const size_t lx = xx.size();
	  for (size_t i = 1; i < lx; ++i) {
	    for (size_t j = 1; j < lx; ++j) { 
        dx[i] = M[i,j] * xx[j];
    }
  }
  
};


// [[Rcpp::export]]
NumericVector DAMOCLES_integrate_odeint(NumericVector ry, 
                                  NumericVector times, 
                                  NumericVector pars, 
                                  double atol, 
                                  double rtol,
                                  std::string stepper) 
{
  std::vector<double> y(ry.size(), 0.0);
  auto rhs_obj = ode_rhs(pars);
  odeint_helper::integrate(stepper, std::ref(rhs_obj), y, times[0], times[1], 0.1 * (times[1] - times[0]), atol, rtol);
  return NumericVector(y.cbegin(), y.cend());
}
