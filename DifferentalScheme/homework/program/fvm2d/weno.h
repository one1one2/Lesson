#ifndef __weno_h_
#define __weno_h_

#include <vector>

//WENO部分代码来自师姐提供
double WENO3CellR(double v3,double v2,double v1);
double WENO3CellL(double v1,double v2,double v3); 
double WENO5CellL(double v5,double v4,double v3,double v2,double v1);
double WENO5CellR(double v1,double v2,double v3,double v4,double v5);
double WENO7CellL(double v7,double v6,double v5,double v4,double v3,
      double v2,double v1);
double WENO7CellR(double v1,double v2,double v3,double v4,double v5,
      double v6,double v7);
double WENO9CellR(double v1,double v2,double v3,double v4,double v5,
      double v6,double v7,double v8,double v9);
double WENO9CellL(double v9,double v8,double v7,double v6,double v5,
      double v4,double v3,double v2,double v1);

//
void WENO_Reconstruct(int order, const std::vector<double>& u, double& ul, double& ur);

#endif
