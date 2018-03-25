#include "hyper2d.h"

typedef Hyper_2D::vartype vartype;
typedef Hyper_2D::soltype soltype;

const double nu = 0.1;
double f(double t, double x, double y){
  return fabs(sin(x))*cos(y)*exp(-2*nu*t);
}

double g(double t, double x, double y){
  return -cos(x)*fabs(sin(y))*exp(-2*nu*t);
}

int main(int argc, char **argv){
  int M = atoi(argv[1]);
  int N = atoi(argv[2]);

  soltype U0(M,std::vector<vartype>(N));
  double xl = 0, xr = M_PI;
  double hx = (xr - xl)/M;
  double yd = 0, yu = M_PI;
  double hy = (yu - yd)/N;
  for(int i = 0; i < M; ++i){
    for(int j = 0; j < N; ++j){
      double x = xl + i * hx;
      double y = yd + j * hy;
      //U0[i][j] = 1-sin(x)*sin(y);
      //U0[i][j] = 1 - (cos(x)-cos(x+hx))*(cos(y)-cos(y+hy))/hx/hy;
      //U0[i][j] = 1 - sin(2*x)*sin(2*y);
      U0[i][j] = 1 - (x-M_PI/2)*(x-M_PI/2)*(y-M_PI/2)*(y-M_PI/2);
    }
  }

  Hyper_2D pro;
  pro.SetInitial(U0);
  pro.SetCoffFunction(f, g);

  pro.Solve(atof(argv[3]));
	
  std::ofstream out("1.dat");
  pro.PrintSolution(std::cout);
  return 0;
}

  
