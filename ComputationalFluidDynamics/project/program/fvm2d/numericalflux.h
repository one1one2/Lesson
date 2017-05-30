#ifndef __numericalflux_h_
#define __numericalflux_h_

#include <cmath> 

#define TEMPLATE template<class vartype, class RealFlux>

TEMPLATE
void LF(const vartype& ul, const vartype& ur, const RealFlux& f,
      double lambda, vartype& F){
  vartype fl, fr;
  f(ul, fl);
  f(ur, fr);
  F = 1./2 * (fl + fr - 1./lambda*CFL*(ur - ul));
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
// For Euler Equation 2D
  vartype fl, fr;
  f(ul, fl); 
  f(ur, fr);
  double rhoL = ul[0], rhoR = ur[0];
  double uL = ul[1]/ul[0], uR = ur[1]/ur[0];
  double pL = (gamma - 1)*(ul[3] - 0.5*uL*uL*rhoL);
  double pR = (gamma - 1)*(ur[3] - 0.5*uR*uR*rhoR);
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
