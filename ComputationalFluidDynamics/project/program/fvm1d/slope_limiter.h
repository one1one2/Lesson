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
  if (r<=0) return 0;
  return (2.*r)/(1.+r);
}

#endif 
