/**
 * @file Heat.h
 * @brief 
 * @author  lczheng, lczheng@pku.edu.cn 
 * 
 * @date 2016-11-08
 */
#ifndef __Heat_h__
#define __Heat_h__

#include <iostream>
#include <vector>
#include "mpi.h"

class Heat{
private:
//	enum Boundrary { DIRICHLET =1 , NEUMANN = 2};
	typedef double (*RHS)(double, double, double, double); 
	typedef double (*RHF)(double, double, double);
	typedef std::vector<double> Point;
	RHS f,u,g;
	RHF u0;
    int N;
	double h;
	Vector<double> sol;
public:
	Heat(){};
	void set_f(const RHS &fun);
	void set_Initial(const RHF &fun);
	void set_Boundrary(int flag , const RHS &fun);
	void set_N(int N);
	void solve(double t);
private:
	Point transform(int N);
    void init();
	void onestep();
};


#endif
