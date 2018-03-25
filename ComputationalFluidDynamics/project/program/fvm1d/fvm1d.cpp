/**
 * @file fvm1d.cpp
 * @brief 
 * @author  lczheng, lczheng@pku.edu.cn 
 * 
 * @date 2017-05-31
 */
#include "config.h"
#include "fvm1d.h"
#include "slope_limiter.h"
#include "numericalflux.h"

#define This FVM_1D

#define N_thread 6

void This::SetRegion(double _xl, double _xr){
  xl = _xl;
  xr = _xr;
}

void This::SetInitial(const soltype& _U){
  U = _U;
  UL = _U;
  UR = _U;
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
#ifdef Euler
  double p;
  maxspeed = 0;
  for (int i = 0; i < N; ++i){
    vartype u = U[i];
    p = (gamma - 1)*(u[2] - 0.5*u[1]*u[1]/u[0]);
    maxspeed = std::max(maxspeed, fabs(u[1]/u[0]) + sqrt(gamma*p/u[0]));
  }
#else 
  maxspeed = 1;
#endif
}

void This::GetTimestep(){
  dt = CFL*h/maxspeed;
}



void This::Reconstruction(const soltype& _U){
#pragma omp parallel for num_threads(N_thread)
  for (int i = 0; i < N; ++i){
    Reconstruction(i, _U);
  }
}

void This::Reconstruction(int i, const soltype& _U){
  vartype ur, u, ul;
  int l, r;
#ifdef Euler 
  l = (i == 0) ? i : i - 1;
  r = (i == N - 1) ? i : i + 1;
#else 
  l = (i == 0) ? N - 1 : i - 1;
  r = (i == N - 1) ? 0 : i + 1;
#endif 
  ul = _U[l]; u = _U[i]; ur = _U[r];

  for (int k = 0; k < u.size(); ++k){
    double ratio = 0;
    if (fabs(u[k] - ur[k]) > 1e-16){
    ratio = (u[k] - ul[k])/(ur[k] - u[k]);
    ratio = limiter(ratio)*(ur[k] - u[k]);
  }
  UL[i][k] = u[k] - ratio/2;
  UR[i][k] = u[k] + ratio/2;
  }
}


void This::GetNumericalFlux(const soltype& _U, 
      soltype& F){
#ifdef Reconstruct
  Reconstruction(_U); 
#else 
  UL = _U;
  UR = _U;
#endif 
  soltype FL(N);
  F.resize(N);
#pragma omp parallel for num_threads(N_thread)
  for (int i = 0; i < N; ++i){
#ifdef Euler 
    int l = (i == 0) ? i : i - 1;
#else 
    int l = (i == 0) ? N - 1 : i - 1;
#endif 
    NumericalFlux(UR[l], UL[i], FL[i]);
  }
#pragma omp parallel for num_threads(N_thread)
  for (int i = 0; i < N; ++i){
#ifdef Euler 
    int r = (i == N - 1) ? i : i + 1;
#else 
    int r = (i == N - 1) ? 0 : i + 1;
#endif
    F[i] = (FL[i] - FL[r])/h;
  }
}

void This::ForwardOnestep(){
  soltype F1, F2, F3, F4;
  soltype U1(U), U2(U), U3(U);
  GetNumericalFlux(U, F1); 

#ifdef RK 
//#pragma omp parallel for num_threads(N_thread)
  //for (int i = 0; i < N; ++i){
    //U1[i] += dt*F1[i];
  //}
  //GetNumericalFlux(U1, F2);

//#pragma omp parallel for num_threads(N_thread)
  //for (int i = 0; i < N; ++i){
    //U[i] += dt/2*(F1[i] + F2[i]);
  //}
//#ifdef RK3
#pragma omp parallel for num_threads(N_thread)
  for (int i = 0; i < N; ++i){
    U1[i] += dt/3*F1[i];
  }
  GetNumericalFlux(U1, F2);
#pragma omp parallel for num_threads(N_thread)
  for (int i = 0; i < N; ++i){
    U2[i] += 2.*dt/3*F2[i];
  }
  GetNumericalFlux(U2, F3);
#pragma omp parallel for num_threads(N_thread)
  for (int i = 0; i < N; ++i){
    U[i] += dt/4*(F1[i] + 3.0*F3[i]);
  }  
#else  
#pragma omp parallel for num_threads(N_thread)
  for (int i = 0; i < N; ++i){
    U[i] += dt*F1[i]; 
  }
#endif
}

void This::NumericalFlux(const vartype& ul, const vartype& ur, vartype& F){
  LW(ul, ur, realflux, dt/h, F);
  //GForce(ul, ur, realflux, dt/h, 1./(1.+CFL), F);
//#ifdef  Euler
  //double pl = (gamma - 1)*(ul[2] - 0.5*ul[1]*ul[1]/ul[0]);
  //double pr = (gamma - 1)*(ur[2] - 0.5*ur[1]*ur[1]/ur[0]);
  //double SL, SR;  
  //SL = std::min(ur[1]/ur[0] - sqrt(gamma*pr/ur[0]), ul[1]/ul[0] - sqrt(gamma*pl/ul[0]));
  //SR = std::max(ur[1]/ur[0] + sqrt(gamma*pr/ur[0]), ur[1]/ul[0] + sqrt(gamma*pl/ul[0]));
  //HLL(ul, ur, realflux, SL, SR, F); 
//#endif 
}



#undef N_thread
#undef This 
