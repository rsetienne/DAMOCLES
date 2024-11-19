//' @useDynLib DAMOCLES


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
  ode_rhs(const NumericMatrix& pars) : M(pars)
  {
    
  }
  
  void operator()(const std::vector<double>& xx, 
                  std::vector<double>& dx, 
                  double /* t */)
  {
    //  dp = M %*% p
    //return(list(dp))
  
    //dx.front() = dx.back() = 0.0;
    const size_t lx = xx.size();
	  for (size_t i = 0; i < lx; ++i) {
	    dx[i] = 0.0;
	    for (size_t j = 0; j < lx; ++j) { 
        dx[i] += M(i, j) * xx[j];
      }
    }
  }
  const NumericMatrix& M;
};


// [[Rcpp::export]]
NumericVector DAMOCLES_integrate_odeint(const NumericVector& ry, 
                                        const NumericVector& times, 
                                        const NumericMatrix& M, 
                                        double atol, 
                                        double rtol,
                                        std::string stepper) 
{
  std::vector<double> y(ry.size(), 0.0);
  std::copy(ry.begin(), ry.end(), y.begin());
  auto rhs_obj = ode_rhs(M);
  odeint_helper::integrate(stepper, std::ref(rhs_obj), y, times[0], times[1], 0.1 * (times[1] - times[0]), atol, rtol);
  return NumericVector(y.cbegin(), y.cend());
}
