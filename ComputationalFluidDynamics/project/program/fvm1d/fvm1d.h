/**
 * @file fvm1d.h
 * @brief Finite Volume Method 1D
 * @author  lczheng, lczheng@pku.edu.cn 
 * 
 * @date 2017-05-31
 */
#ifndef __FVM_1D_H
#define __FVM_1D_H

#include "config.h"

class FVM_1D{
public:
#ifndef Euler 
  typedef Vector<double,1> vartype;
#else 
  typedef Vector<double,3> vartype;
#endif
  typedef std::vector<vartype> soltype;
  typedef std::function<void(const vartype&, vartype&)> RealFlux;
private:
  double xl, xr;
  int N;
  double h, dt, t_now;
  double maxspeed;

  RealFlux realflux;
  soltype U, UL, UR;

public:
  FVM_1D() = default;
  FVM_1D(const FVM_1D&) =delete;
  FVM_1D& operator=(const FVM_1D&) = delete;

  void SetRegion(double, double);
  void SetInitial(const soltype&);
  void SetRealflux(const RealFlux&);
  
  void Solve(double);
  template<typename stream>
  void PrintSolution(stream&);
private:
  void Initialize();
  void Reconstruction(const soltype&);
  void Reconstruction(int, const soltype&);
  void ForwardOnestep();
  void ForwardOnestep(int);
  void GetNumericalFlux(const soltype&, soltype&);
  void NumericalFlux(const vartype&, const vartype&, vartype&);
  void GetMaxspeed();
  void GetTimestep();
  void BoundaryCondition(int, vartype&, vartype&, vartype&);

};

template<typename stream>
void FVM_1D::PrintSolution(stream& os){
  for (int i = 0; i < N; i++)
    os << U[i] << "\n";
}




#endif 
