#include "fvm1d.h"

typedef FVM_1D::vartype vartype;
typedef FVM_1D::soltype soltype;

void flux(const vartype& u, vartype& f){
  f = u;
#ifdef Euler 
  double p;
  p = (gamma - 1)*(u[2] - 0.5*u[1]*u[1]/u[0]);
  f[0] = u[1];
  f[1] = u[1]*u[1]/u[0] + p; 
  f[2] = u[1]*(u[2] + p)/u[0];
#endif
}


double F(double x, double alpha, double a){
  if (fabs(alpha*(x-a)) <= 1.)
    return sqrt(1. - alpha*alpha*(x-a)*(x-a));
  else return 0;
}

double G(double x, double beta, double z){
  return exp(-beta*(x-z)*(x-z));
}
void u0(double x, vartype& u){
#ifndef Euler
  double a = 0.5, z = -0.7, delta = 0.005, alpha = 10., beta = log(2)/36/delta/delta;
  u = 0;
  if (x >= -0.8 && x <= -0.6)
    u = 1./6*(G(x,beta,z-delta) + G(x,beta,z+delta) + 4*G(x,beta,z));
  if (x >= -0.4 && x <= -0.2)
    u = 1;
  if (x >= 0    && x <=  0.2)
    u = 1 - 10.*fabs(x - 0.1);
  if (x >= 0.4  && x <=  0.6)
    u = 1./6*(F(x,alpha,a-delta) + F(x,alpha,a+delta) + 4*F(x,alpha,a));
#else 
  if (x < -4.0){
    u[0] = 3.857143;
    u[1] = 2.629369*u[0];
    u[2] = 10.33333/(gamma - 1) + 0.5*u[1]*u[1]/u[0];
  }
  else{  
    u[0] = 1 + 0.2*sin(5*x);
    u[1] = 0;
    u[2] = 1/(gamma - 1);
  }
#endif 
}

int main(int argc, char **argv){
  double t_begin = omp_get_wtime();
  int N = atoi(argv[1]);

  soltype U0(N);
#ifdef Euler 
  double xl = -5.0;
  double xr = 5.0;
#else 
  double xl = -1.0;
  double xr = 1.0;
#endif
  double h = (xr - xl)/N;

  for (int i = 0; i < N; ++i){
    double x = xl + i*h;
    u0(x, U0[i]);
  }

  FVM_1D pro;
  pro.SetRegion(xl, xr);
  pro.SetInitial(U0);
  pro.SetRealflux(flux);

  pro.Solve(atof(argv[2]));
	
  std::ofstream out("1.dat");
  pro.PrintSolution(out);
  double t_end = omp_get_wtime();
  std::cout << (t_end - t_begin) << std::endl;
  return 0;
}

  
