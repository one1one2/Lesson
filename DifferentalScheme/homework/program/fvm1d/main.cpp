#include "fvm1d.h"

typedef FVM_1D::vartype vartype;
typedef FVM_1D::soltype soltype;

void flux(const vartype& u, vartype& f){
  f = 1./2*u*u;
}

int main(int argc, char **argv){
  double t_begin = clock();
  int N = atoi(argv[1]);
  double ul = atof(argv[2]);
  double um = atof(argv[3]);
  double ur = atof(argv[4]);

  soltype U0(N);
  
  double ratiol = 2./5;
  double ratior = 1 - ratiol;
  for(int i = 0; i < N; ++i){
    if (i < ratiol * N)
      U0[i] = ul;
    if (i >= ratiol * N && i < ratior * N)
      U0[i] = um;
    if (i >= ratior * N)
      U0[i] = ur;
  }

  FVM_1D pro;
  pro.SetInitial(U0);
  pro.SetRealflux(flux);
  pro.SetScheme(atoi(argv[5]));
  pro.Solve(atof(argv[6]));
	
  pro.PrintSolution(std::cout);
  double t_end = clock();
  std::cerr << (t_end - t_begin)/CLOCKS_PER_SEC << std::endl;
  return 0;
}

  
