#ifndef __slope_limiter_h_
#define __slope_limiter_h_

#include <cmath>

double minmod(double r){
  return std::max(0.0, std::min(r, 1.0));
}

double superbee(double r){
  return std::max(0., std::max(
          std::min(1., 2*r), std::min(2., r)));
}
  
double MC(double r){
  return std::max(0., std::min((1.+r)/2, std::min(2., 2.*r)));
}

double vanLeer(double r){
  return (r + fabs(r))/(1+fabs(r));
}

#endif 
