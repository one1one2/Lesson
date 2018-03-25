#ifndef __Hyper_2D_H
#define __Hyper_2D_H

#include <iostream>
#include <cmath>
#include <ostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <cassert>
#include <functional>

class Hyper_2D{
public:
  typedef double vartype;
  typedef std::vector<std::vector<vartype>> soltype;
  typedef std::function<double(double,double,double)> CoffFunction;

private:
  const double CFL = 0.3;
  const double xl = 0;
  const double xr = M_PI;
  const double yu = M_PI;
  const double yd = 0;

  const int WENO_order = 5;
  int M, N;
  double hx, hy, dt, t_now;
  double maxspeed;

  CoffFunction coff_x, coff_y;
  soltype U, U_temp, UL, UR, UU, UD;

public:
  Hyper_2D() = default;
  Hyper_2D(const Hyper_2D&) =delete;
  Hyper_2D& operator=(const Hyper_2D&) = delete;

  void SetInitial(const soltype&);
  void SetCoffFunction(const CoffFunction&, const CoffFunction&);
  
  void Solve(double);
  template<typename stream>
    void PrintSolution(stream&);
private:
  void Initialize();
  void ForwardOnestep();
  void GetNumericalFlux(const soltype&, soltype&);
  void Upwind(int, int);
  void GetLR(const soltype&, int, int);
  void GetMaxspeed();
  void GetTimestep();
};

template<typename stream>
void Hyper_2D::PrintSolution(stream& os){
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j){
      os << U[i][j] << " ";
    }
    os << "\n";
  }
}



#endif 
