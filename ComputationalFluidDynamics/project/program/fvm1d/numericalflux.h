#ifndef __numericalflux_h_
#define __numericalflux_h_

#include <cmath> 
#include "fvm1d.h"

#define TEMPLATE template<class vartype, class RealFlux>

TEMPLATE
void LF(const vartype& ul, const vartype& ur, const RealFlux& f,
      double lambda, vartype& F){
  vartype fl, fr;
  f(ul, fl);
  f(ur, fr);
  F = 1./2 * (fl + fr - 1/lambda*(ur - ul));
}

TEMPLATE
void LW(const vartype& ul, const vartype& ur, const RealFlux& f, 
      double lambda, vartype& F){
  vartype fl, fr;
  f(ul, fl);
  f(ur, fr);
  F = 0.5 * (fl + fr) - 0.5 * lambda * (fr - fl); 
}

TEMPLATE
void LW_twostep(const vartype& ul, const vartype& ur, const RealFlux& f,
      double lambda, vartype& F){
  vartype fl, fr;
  f(ul, fl);
  f(ur, fr);
  vartype u_temp = 0.5*(ul + ur) - 0.5*lambda*(fr - fl);
  f(u_temp, F);
}

TEMPLATE
void Force(const vartype& ul, const vartype& ur, const RealFlux& f,
      double lambda, vartype& F){
  vartype F1, F2;
  LF(ul, ur, f, lambda, F1);
  LW_twostep(ul, ur, f, lambda, F2);
  F = 1./2* (F1 + F2);
}

TEMPLATE
void GForce(const vartype& ul, const vartype& ur, const RealFlux& f,
      double lambda, double w, vartype& F){
  vartype F1, F2;
  LF(ul, ur, f, lambda, F1);
  LW_twostep(ul, ur, f, lambda, F2);
  F = (1-w)*F1 + w*F2;
}

TEMPLATE
void Upwind(const vartype& ul, const vartype& ur, const RealFlux& f,
      double lambda, vartype& F){
  vartype fl, fr;
  f(ul, fl);
  f(ur, fr);
  F = fl;
  for (int i = 0; i < ul.size(); ++i){
    if (fabs(ur[i]-ul[i]) > 1e-10 && (fr[i]-fl[i])/(ur[i]-ul[i]) < 0){
      F[i] = fr[i];
    }
  }
}

//TEMPLATE
//void Roe(const vartype& ul, const vartype& ur, const RealFlux& f,
      //vartype& F){
  //double rhoL = ul[0], rhoR = ur[0];
  //double uL = ul[1]/ul[0], uR = ur[1]/ur[0];
  //double pL = (gamma - 1)*(ul[2] - 0.5*uL*ul[1]);
  //double pR = (gamma - 1)*(ur[2] - 0.5*uR*ur[1]);
  //double HL = (ul[2] + pL)/rhoL;
  //double HR = (ur[2] + pR)/rhoR;
  //double rho = sqrt(rhoL*rhoR);
  //double u = (sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR));
  //double H = (sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR));
  //double a = sqrt((gamma-1)*(H-0.5*u*u));
  //double lambda1L = uL - sqrt(gamma*pL/rhoL);
  //double lambda1R = uR - sqrt(gamma*pR/rhoR);
  //double lambda1 = u - sqrt(gamma*p/rho);
  //double lambda1bar = lambda1L*(lambda1R-lambda1)/(lambda1R-lambda1L);
  //vartype K1;
  //K1[0] = 1;
  //K1[1] = u - a;
  //K1[2] = H - u*a;
  //F = fl + lambda1bar*alpha1*K1;
//}

TEMPLATE 
void HLL(const vartype& ul, const vartype& ur, const RealFlux& f,
      double Sl, double Sr, vartype& F){
  vartype fl, fr;
  f(ul, fl); 
  f(ur, fr);
  if (Sl >= 0) F = fl;
  else if (Sr <= 0) F = fr;
  else F = 1./(Sr - Sl)*(Sr*fl - Sl*fr + Sl*Sr*(ur-ul));
}

TEMPLATE 
void HLLC(const vartype& ul, const vartype& ur, const RealFlux& f,
      double Sl, double Sr, vartype& F){
// For Euler Equation
  vartype fl, fr;
  f(ul, fl); 
  f(ur, fr);
  double rhoL = ul[0], rhoR = ur[0];
  double uL = ul[1]/ul[0], uR = ur[1]/ur[0];
  double pL = (gamma - 1)*(ul[2] - 0.5*uL*ul[1]);
  double pR = (gamma - 1)*(ur[2] - 0.5*uR*ur[1]);
  double S;
  
  S = (pR-pL+rhoL*uL*(Sl-uL)-rhoR*uR*(Sr-uR))/(rhoL*(Sl-uL)-rhoR*(Sr-uR));

  vartype D(0.);
  D[1] = 1.0; D[2] = S;

  if (Sl >= 0) F = fl;
  else if (Sr <= 0) F = fr;
  else if (S >= 0)
    F = (S*(Sl*ul-fl)+Sl*(pL+rhoL*(Sl-uL)*(S-uL))*D)/(Sl-S);
  else 
    F = (S*(Sr*ur-fr)+Sr*(pR+rhoR*(Sr-uR)*(S-uR))*D)/(Sr-S);

}


TEMPLATE
void Godunov(const vartype& ul, const vartype& ur, const RealFlux& f,
      double lambda, vartype& F){
  vartype fl, fr;
  f(ul, fl);
  f(ur, fr);
  F = fl;
}
 


#undef TEMPLATE

#endif
