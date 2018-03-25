#include "fvm1d.h"

#define This FVM_1D

void This::SetInitial(const soltype& _U){
  U = _U;
  U_temp = _U;
  N = _U.size();
  h = (xr-xl)/N;
}

void This::SetRealflux(const RealFlux& _realflux){
  realflux = _realflux;
}

void This::Solve(double t_end){
  Initialize(); 
  while (t_now < t_end){
    GetMaxspeed();
    GetTimestep();
    dt = std::min(dt, t_end - t_now);
    std::cerr << "t = " << t_now << ", dt = " << dt << std::endl;
    ForwardOnestep();
    t_now += dt;
  }
  std::cerr << "Finished!" << std::endl;
}

void This::Initialize(){
  std::cerr << "Initializing..." << std::endl;
  t_now = 0;
}

void This::GetMaxspeed(){
  maxspeed = 2;
}

void This::GetTimestep(){
  dt = CFL*h/maxspeed;
}

#define N_thread 8

void This::ForwardOnestep(){
#pragma omp parallel for num_threads(N_thread)
  for (int i = 0; i < N; ++i){
    ForwardOnestep(i);
  }

#pragma omp parallel for num_threads(N_thread)
  for (int i = 0; i < N; ++i){
    U[i] = U_temp[i];
  }
}

#undef N_thread

void This::ForwardOnestep(int i){
  vartype ur, u, ul;
  BoundaryCondition(i, ul, u, ur);
  ForwardOnestep(ul, u, ur, U_temp[i]);
}

void This::ForwardOnestep(const vartype& ul, const vartype& u, 
      const vartype& ur, vartype& u1){
  vartype Fl, Fr;
  NumericalFlux(ul, u, Fl);
  NumericalFlux(u, ur, Fr);
  u1 = u - dt/h * (Fr - Fl);
}

void This::NumericalFlux(const vartype& ul, const vartype& ur, vartype& F){
  vartype fl, fr;
  realflux(ul, fl);
  realflux(ur, fr);
//Lax-Friedich
  if (Scheme == 1)
  F = 1./2 * (fl + fr - maxspeed/CFL*(ur - ul));
//Lax-Wendroff 仅针对Burgers方程 
  else if (Scheme == 2){
  //if (ul == ur)
  //F = fl;
  //else 
  //F = 0.5 * (fl + fr) - 0.5*dt/h *(ur+ul)* (fr - fl); 
//Richtmyer 
  vartype u_temp = 0.5*(ul + ur) - 0.5*dt/h*(fr - fl);
  realflux(u_temp, F);}
//Godunov 格式 仅针对Burgers方程 
  else{
  if (ul < ur){
    if (ul > 0) F = fl;
    else if (ur < 0) F = fr;
    else F = 0;
  }else{
    if ((ul + ur)/2 > 0) F = fl;
    else F = fr;
  } 
  }
}

void This::BoundaryCondition(int i, vartype& ul, vartype& u, vartype& ur){
  u = U[i];
  if (BoundaryFlag(i) == 1)
    ul = u;
  else{
    ul = U[i - 1];
  }
  if (BoundaryFlag(i) == 2)
    ur = u;
  else{	
    ur = U[i + 1];
  }
}

int This::BoundaryFlag(int i){
  if (i == 0) 
    return 1;
  if (i == N - 1)
    return 2;
  return 0;
}




#undef This 
