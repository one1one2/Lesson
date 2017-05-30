#ifndef __RK_h_
#define __RK_h_

#define TEMPLATE template<class vartype, class Function>

TEMPLATE
void RK4(double t, const vartype& y, const Function& f, double dt, vartype& y1){
  vartype K1, K2, K3, K4;
  f(t, y, K1);
  f(t + 1./2*dt, y + 1./2*dt*K1, K2);
  f(t + 1./2*dt, y + 1./2*dt*K2, K3);
  f(t + dt, y + dt*K3, K4);
  y1 = y + h/6*(K1 + 2*K2 + 2*K3 + K4);
}





#undef TEMPLATE
#endif 
