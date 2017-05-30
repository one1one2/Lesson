#include "fvm2d.h"

typedef FVM_2D::vartype vartype;
typedef FVM_2D::soltype soltype;

void flux_x(const vartype& u, vartype& f){
  double p;
  p = (gamma - 1)*(u[3] - 0.5*u[1]*u[1]/u[0] - 0.5*u[2]*u[2]/u[0]);
  f[0] = u[1];
  f[1] = u[1]*u[1]/u[0] + p; 
  f[2] = u[1]*u[2]/u[0];
  f[3] = u[1]/u[0]*(u[3] + p);
}

void flux_y(const vartype& u, vartype& g){
  double p;
  p = (gamma - 1)*(u[3] - 0.5*u[1]*u[1]/u[0] - 0.5*u[2]*u[2]/u[0]);
  g[0] = u[2];
  g[1] = u[1]*u[2]/u[0];
  g[2] = u[2]*u[2]/u[0] + p;
  g[3] = u[2]/u[0]*(u[3] + p);
}

void u0(double x, double y, vartype& u){
if (Problem==DoubleMachReflection){
  vartype ul, ur;
  ul[0] = 8.;
  ul[1] = 57.1597;
  ul[2] = -33.0012;
  ul[3] = 563.544;
  ur[0] = 1.4;
  ur[1] = 0.;
  ur[2] = 0.;
  ur[3] = 2.5;
  if (y >= sqrt(3.)*(x - 1./6)) u = ul;
  else u = ur;
}else if (Problem==ForwardStep){
  u[0] = 1.4;
  u[1] = 3*u[0];
  u[2] = 0;
  u[3] = 1./(gamma - 1) + 0.5*u[1]*u[1]/u[0];
}
}



int main(int argc, char **argv){
  int M = atoi(argv[1]);
  int N = atoi(argv[2]);

  soltype U0(M,std::vector<vartype>(N));
  double xl,xr;
if (Problem==DoubleMachReflection){
  xl = 0, xr = 4.0;
}else if (Problem==ForwardStep){
  xl = 0, xr = 3.0;
}
  double hx = (xr - xl)/M;
  double yd = 0, yu = 1.0;
  double hy = (yu - yd)/N;
  for(int i = 0; i < M; ++i){
    for(int j = 0; j < N; ++j){
      double x = xl + (i+0.5) * hx;
      double y = yd + (j+0.5) * hy;
      u0(x, y, U0[i][j]);
    }
  }

  FVM_2D pro;
  pro.SetRegion(xl, xr, yd, yu);
  pro.SetInitial(U0);
  pro.SetRealflux(flux_x, flux_y);

  pro.Solve(atof(argv[3]));
	
  std::ofstream out("rho.dat");
  pro.PrintSolution(out, 0);
  return 0;
}

  
