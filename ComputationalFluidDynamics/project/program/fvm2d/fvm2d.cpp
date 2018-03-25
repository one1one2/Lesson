#include "fvm2d.h"
#include "numericalflux.h"
#include "slope_limiter.h"

#define This FVM_2D

#define N_thread 4

void This::SetRegion(double _xl, double _xr,
      double _yd, double _yu){
  xl = _xl;
  xr = _xr;
  yd = _yd;
  yu = _yu;
}            

void This::SetInitial(const soltype& _U){
  U = _U;
  UL = U; UR = U; UU = U; UD = U;
  U_temp = _U;
  M = _U.size();
  N = U[0].size();
}

void This::SetRealflux(const Realflux& _flux_x,
      const Realflux& _flux_y){
  flux_x = _flux_x;
  flux_y = _flux_y;
}

void This::Solve(double t_end){
  Initialize(); 
  while (t_now < t_end){
    GetMaxspeed();
    GetTimestep();
    dt = std::min(dt, t_end - t_now);
    ForwardOnestep();
    t_now += dt;
  }
  std::cerr << "Finished!" << std::endl;
}

void This::Initialize(){
  std::cerr << "Initializing..." << std::endl;
  t_now = 0;
  hx = (xr - xl)/M;
  hy = (yu - yd)/N;
}

void This::GetMaxspeed(){
  double p;
  maxspeed = 0;
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j){
      vartype u = U[i][j];
      p = (gamma - 1)*(u[3] - 0.5*u[1]*u[1]/u[0] - 0.5*u[2]*u[2]/u[0]);
      maxspeed = std::max(maxspeed, fabs(u[1]/u[0]) + sqrt(gamma*p/u[0]));
      maxspeed = std::max(maxspeed, fabs(u[2]/u[0]) + sqrt(gamma*p/u[0]));
    }
  }
  //maxspeed *= 2; 
  std::cout << "MaxSpeed = " << maxspeed<<std::endl; 
}

void This::GetTimestep(){
  dt = CFL*std::min(hx, hy)/maxspeed;
}

void This::ForwardOnestep(){
  std::cerr << "t = " << t_now << ", dt = " << dt << std::endl;
  soltype F1, G1, F2, G2, F3, G3, F4, G4, U1, U2;
  U1 = U; U2=U;
  GetNumericalFlux(U, F1, G1);
#ifdef RK
#pragma omp parallel for num_threads(N_thread)
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j){
      U1[i][j] += dt/3*F1[i][j] + dt/3*G1[i][j];
    }
  }
  GetNumericalFlux(U1, F2, G2);
#pragma omp parallel for num_threads(N_thread)
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j){
      U2[i][j] += 2*dt/3*F2[i][j] + 2*dt/3*G2[i][j];
    }
  }
  GetNumericalFlux(U2, F3, G3);
#pragma omp parallel for num_threads(N_thread)
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j){
      U_temp[i][j] += dt/4*(F1[i][j] + G1[i][j] + 3.*F3[i][j] + 3.*G3[i][j]);
    }
  }
#else  
#pragma omp parallel for num_threads(N_thread)
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j){
      U_temp[i][j] += dt*(F1[i][j] + G1[i][j]);
    }
  }
#endif   
#pragma omp parallel for num_threads(N_thread)
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j){
      if (Problem != ForwardStep || (xl + i*hx <= 0.6 || yd + j*hy >= 0.2)) 
        U[i][j] = U_temp[i][j];
    }
  }
}

void This::Reconstruction(const soltype& _U){
#pragma omp parallel for num_threads(N_thread)
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j){
      Reconstruction(i, j, _U);
    }
  }
}

void This::Reconstruction(int i, int j, const soltype& _U){
  vartype ur, u, ul, ud, uu;
  int l, r, d, _u;
  l = (i == 0) ? i : i - 1;
  r = (i == M - 1) ? i : i + 1;
  d = (j == 0) ? j : j - 1;
  _u = (j == N - 1) ? j : j + 1;
  ul = _U[l][j]; 
  u = _U[i][j]; 
  ur = _U[r][j]; 
  ud = _U[i][d]; 
  uu = _U[i][_u];
  for (int k = 0; k < u.size(); ++k){
    double ratio = 0;
    if (fabs(u[k] - ur[k]) > 1e-16){
      ratio = (u[k] - ul[k])/(ur[k] - u[k]);
      ratio = limiter(ratio)*(ur[k] - u[k]);
    }
    UL[i][j][k] = u[k] - ratio/2;
    UR[i][j][k] = u[k] + ratio/2;
    ratio = 0;
    if (fabs(u[k] - uu[k]) > 1e-16){
      ratio = (u[k] - ud[k])/(uu[k] - u[k]);
      ratio = limiter(ratio)*(uu[k] - u[k]);
    }
    UD[i][j][k] = u[k] - ratio/2;
    UU[i][j][k] = u[k] + ratio/2;
  }
}

void This::GetNumericalFlux(const soltype& _U, soltype& F, soltype& G){
#ifdef Reconstruct
  Reconstruction(_U);
#else 
  UL = _U;  UR = _U;
  UU = _U;  UD = _U;
#endif 
  vartype ul, ur;

if (Problem==DoubleMachReflection){
  ul[0] = 8.;
  ul[1] = 57.1597;
  ul[2] = -33.0012;
  ul[3] = 563.544;
  ur[0] = 1.4;
  ur[1] = 0.;
  ur[2] = 0.;
  ur[3] = 2.5;
}else if (Problem==ForwardStep){
  ul[0] = 1.4;
  ul[1] = 3*ul[0];
  ul[2] = 0;
  ul[3] = 1./(gamma - 1) + 0.5*ul[1]*ul[1]/ul[0];
}
  vartype Fl, Fr, Gd, Gu;
  vartype u_temp;
  F = _U; G = _U;
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j){
      //Left   入流边界条件  
      if (i == 0) u_temp = ul; 
      else u_temp = UR[i-1][j];
      NumericalFlux_x(u_temp, UL[i][j], Fl);
      //Right 出流边界条件 
      if (i == M - 1) u_temp = UR[i][j];
      else u_temp = UL[i+1][j];
//前台阶问题的反射边界
if (Problem==ForwardStep){
      if (xl + i*hx == 0.6 && yd + j*hy <= 0.2){
        u_temp = UR[i][j];
        u_temp[1] *= -1; 
      }
}
      NumericalFlux_x(UR[i][j], u_temp, Fr);
      F[i][j] = (Fl - Fr)/hx;
//Down
if (Problem==DoubleMachReflection){
      if (j == 0){
        if (xl + i*hx >= 1./6){
          u_temp = UD[i][j];
          u_temp[2] *= -1;
          //flux_y(u_temp,Gd);
        }else
          u_temp = ul;
      }
      else u_temp = UU[i][j-1]; 
}else if (Problem==ForwardStep){
    if (j == 0|| (yd + j*hy == 0.2 && xl + i*hx >= 0.6)){
      u_temp = UD[i][j];
      u_temp[2] *= -1;
    }
    else u_temp = UU[i][j-1];
}  
    NumericalFlux_y(u_temp, UD[i][j], Gd);
//Up
if (Problem==DoubleMachReflection){
      if (j == N - 1){
        if (1 > sqrt(3.)*(xl + i*hx - 1./6) - 20*t_now)
          u_temp = ul;
        else u_temp = ur;
      } 
      else u_temp = UD[i][j+1]; 
}else if (Problem==ForwardStep){
      if (j == N - 1){
        u_temp = UU[i][j];
        u_temp[2] *= -1;
      }else u_temp = UD[i][j+1];
}
      NumericalFlux_y(UU[i][j], u_temp, Gu);
      G[i][j] = (Gd - Gu)/hy;
    }
  }
}


void This::NumericalFlux_x(const vartype& ul, const vartype& ur, vartype& F){
  //LF(ul, ur, flux_x, dt/hx, F);
  double pl = (gamma - 1)*(ul[3] - 
        0.5*ul[1]*ul[1]/ul[0] - 0.5*ul[2]*ul[2]/ul[0]);
  double pr = (gamma - 1)*(ur[3] - 
        0.5*ur[1]*ur[1]/ur[0] - 0.5*ur[2]*ur[2]/ur[0]);
  double SL, SR;
  SL = std::min(ur[1]/ur[0] - sqrt(gamma*pr/ur[0]), ul[1]/ul[0] - sqrt(gamma*pl/ul[0]));
  SR = std::max(ur[1]/ur[0] + sqrt(gamma*pr/ur[0]), ul[1]/ul[0] + sqrt(gamma*pl/ul[0]));
  HLL(ul, ur, flux_x, SL, SR, F);
}
  
void This::NumericalFlux_y(const vartype& ul, const vartype& ur, vartype& F){
  //LF(ul, ur, flux_y, dt/hy, F);
  double pl = (gamma - 1)*(ul[3] - 
        0.5*ul[1]*ul[1]/ul[0] - 0.5*ul[2]*ul[2]/ul[0]);
  double pr = (gamma - 1)*(ur[3] - 
        0.5*ur[1]*ur[1]/ur[0] - 0.5*ur[2]*ur[2]/ur[0]);
  double SL, SR;
  SL = std::min(ur[2]/ur[0] - sqrt(gamma*pr/ur[0]), ul[2]/ul[0] - sqrt(gamma*pl/ul[0]));
  SR = std::max(ur[2]/ur[0] + sqrt(gamma*pr/ur[0]), ul[2]/ul[0] + sqrt(gamma*pl/ul[0]));
  HLL(ul, ur, flux_y, SL, SR, F);
}

#undef This 
