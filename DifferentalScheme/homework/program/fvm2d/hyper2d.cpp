#include "hyper2d.h"
#include "weno.h"

#define This Hyper_2D

#define _WENO_Reconstruct

void This::SetInitial(const soltype& _U){
  U = _U;
  U_temp = _U;
  UL = U;
  UR = U;
  UU = U;
  UD = U;
  M = _U.size();
  N = U[0].size();
  hx = (xr - xl)/M;
  hy = (yu - yd)/N;
}

void This::SetCoffFunction(const CoffFunction& _coff_x,
      const CoffFunction& _coff_y){
  coff_x = _coff_x;
  coff_y = _coff_y;
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
}

void This::GetMaxspeed(){
  maxspeed = 1;
}

void This::GetTimestep(){
  dt = CFL*std::min(hx, hy)/maxspeed;
}

void This::ForwardOnestep(){
  std::cerr << "t = " << t_now << ", dt = " << dt << std::endl;

#ifndef _WENO_Reconstruct
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j){
      GetLR(U,i,j);
    }
  }
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j){
      Upwind(i,j);
    }
  }
  U = U_temp;
#endif 
#ifdef _WENO_Reconstruct
  soltype F1, F2, F3;
  soltype U1, U2;
  U1 = U; U2 = U;
  GetNumericalFlux(U, F1);
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      U1[i][j] += dt/3*F1[i][j];
  GetNumericalFlux(U1, F2);
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      U2[i][j] += 2.*dt/3*F2[i][j];
  GetNumericalFlux(U2, F3);
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j){
      U[i][j] += dt/4*(F1[i][j] + 3.*F3[i][j]);
    }
  }
#endif 
}

void This::GetLR(const soltype& _U, int i, int j){

#ifndef _WENO_Reconstruct
  UL[i][j] = _U[(i-1+M)%M][j];
  UR[i][j] = _U[(i+1)%M][j];
  UD[i][j] = _U[i][(j-1+N)%N];
  UU[i][j] = _U[i][(j+1)%N];
#else 
  int order = (WENO_order-1)/2;
  std::vector<vartype> ux(2*order + 1), uy(2*order + 1);
  
  for (int index = 0; index < 2*order + 1; ++index){
    ux[index] = _U[(i + index - order + M)%M][j];
    uy[index] = _U[i][(j + index - order + N)%N];
  }
  WENO_Reconstruct(order, ux, UL[i][j], UR[i][j]);
  WENO_Reconstruct(order, uy, UD[i][j], UU[i][j]);
#endif 
}

void This::GetNumericalFlux(const soltype& _U, soltype& F){
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j){
      GetLR(_U,i,j);
    }
  }
  F = _U;
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j){
      double x = xl + i  * hx;
      double y = yd + j  * hy;
      vartype FL, FR, FD, FU;
      FL = 0.5*(coff_x(t_now, x - 0.*hx, y)*(UL[i][j]+UR[(i-1+M)%M][j]) 
             //- 0.5*dt/hx*coff_x(t_now,x,y)*coff_x(t_now,x,y)
             - 0.0*hx/dt * (UL[i][j] - UR[(i-1+M)%M][j]));
      FR = 0.5*(coff_x(t_now, x + 0.*hx, y)*(UR[i][j]+UL[(i+1+M)%M][j])
             //- 0.5*dt/hx*coff_x(t_now,x,y)*coff_x(t_now,x,y)
            - 0.0*hx/dt* (UL[(i+1+M)%M][j] - UR[i][j]));
      FD = 0.5*(coff_y(t_now, x, y - 0.*hy)*(UD[i][j]+UU[i][(j-1+N)%N]) 
             //- 0.5*dt/hy*coff_y(t_now,x,y)*coff_y(t_now,x,y)
            - 0.0*hy/dt * (UD[i][j] - UU[i][(j-1+N)%N]));
      FU = 0.5*(coff_y(t_now, x, y + 0.*hy)*(UU[i][j]+UD[i][(j+1+N)%N])
             //- 0.5*dt/hy*coff_y(t_now,x,y)*coff_y(t_now,x,y)
             - 0.0*hy/dt* (UD[i][(j+1+N)%N] - UU[i][j]));
      F[i][j] = (FL-FR)/hx + (FD-FU)/hy;
    }
  }
}


void This::Upwind(int i, int j){
  double x = xl + i  * hx;
  double y = yd + j  * hy;

  double a = coff_x(t_now, x, y);
  double v = a * dt / hx;
  if (v <= 0)
    U_temp[i][j] += v * (U[i][j] - UR[i][j]);
  else
    U_temp[i][j] += v * (UL[i][j] - U[i][j]);

  a = coff_y(t_now, x, y);
  v = a * dt / hy;
  if (v <= 0)
    U_temp[i][j] += v * (U[i][j] - UU[i][j]);
  else
    U_temp[i][j] += v * (UD[i][j] - U[i][j]);
}

void WENO_Reconstruct(int order, const std::vector<double>& u, double& ul, double& ur){
  assert(u.size() == 2*order + 1); 
  switch (2*order+1){
    case 3:
      ul = WENO3CellL(u[0], u[1], u[2]);
      ur = WENO3CellR(u[0], u[1], u[2]);
      break;
    case 5:
      ul = WENO5CellL(u[0], u[1], u[2], u[3], u[4]);
      ur = WENO5CellR(u[0], u[1], u[2], u[3], u[4]);
      break;
    case 7: 
      ul = WENO7CellL(u[0], u[1], u[2], u[3], u[4], u[5], u[6]);
      ur = WENO7CellR(u[0], u[1], u[2], u[3], u[4], u[5], u[6]);
      break;
    case 9:
      ul = WENO9CellL(u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8]);
      ur = WENO9CellR(u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8]);
      break;
    default:
      break;
  }
}







#undef This 
