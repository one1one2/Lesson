#ifndef __FVM_2D_H
#define __FVM_2D_H

#include "config.h"
/**
 * @brief Finite Vloume Method for 2D conservation laws 
 *    U_t + (F(U))_x +(G(U))_y = 0  
 */


class FVM_2D{
public:
  typedef Vector<double,4> vartype;
  typedef std::vector<std::vector<vartype>> soltype;
  typedef std::function<void(const vartype&, vartype&)> Realflux;

private:
  double xl, xr, yd, yu;
  
  int M, N;
  double hx, hy, dt, t_now;
  double maxspeed;

  Realflux flux_x, flux_y;
  soltype U, U_temp;
  soltype UL, UR, UD, UU;
public:
  FVM_2D() = default;
  FVM_2D(const FVM_2D&) =delete;
  FVM_2D& operator=(const FVM_2D&) = delete;

  void SetRegion(double, double, double, double);
  void SetInitial(const soltype&);
  void SetRealflux(const Realflux&, const Realflux&);
  
  void Solve(double);
  template<typename stream>
    void PrintSolution(stream&, int);
private:
  void Initialize();
  void ForwardOnestep();
  void Reconstruction(const soltype&);
  void Reconstruction(int, int, const soltype&);
  void GetNumericalFlux(const soltype&, soltype&, soltype&);
  void ForwardOnestep(int, int);
  void GetMaxspeed();
  void GetTimestep();

  void NumericalFlux_x(const vartype&, const vartype&, vartype&);
  void NumericalFlux_y(const vartype&, const vartype&, vartype&);
};

template<typename stream>
void FVM_2D::PrintSolution(stream& os, int index){
  for (int j = 0; j < N; ++j){
    for (int i = 0; i < M; ++i){
      os << U[i][j][index] << " ";
    }
    os << "\n";
  }
}



#endif 
