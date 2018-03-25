#ifndef __FVM_1D_H
#define __FVM_1D_H

#include <iostream>
#include <cmath>
#include <ostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <thread>
#include <functional>
#include <omp.h>
#include <cassert>

class FVM_1D{
public:
  typedef double vartype;
  typedef std::vector<vartype> soltype;
  typedef std::function<void(const vartype&, vartype&)> RealFlux;
private:
  const double CFL = 0.6;
  const double xl = -1.0;
  const double xr = 1.0;
  int N;
  double h, dt, t_now;
  double maxspeed;
  int Scheme;

  RealFlux realflux;
  soltype U, U_temp;

public:
  FVM_1D() = default;
  FVM_1D(const FVM_1D&) =delete;
  FVM_1D& operator=(const FVM_1D&) = delete;

  void SetInitial(const soltype&);
  void SetRealflux(const RealFlux&);
  void SetScheme(int _Scheme){ Scheme = _Scheme;} 

  void Solve(double);
  template<typename stream>
  void PrintSolution(stream&);
private:
  void Initialize();
  void ForwardOnestep();
  void ForwardOnestep(int);
  void ForwardOnestep(const vartype&, const vartype&, const vartype&, vartype&);
  void NumericalFlux(const vartype&, const vartype&, vartype&);
  void GetMaxspeed();
  void GetTimestep();
  int BoundaryFlag(int);
  void BoundaryCondition(int, vartype&, vartype&, vartype&);

};

template<typename stream>
void FVM_1D::PrintSolution(stream& os){
  for (int i = 0; i < N; i++)
    os << U[i] << " ";
}




#endif 
