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
#include <cmath>
#include <cstdlib>
#include "mpi.h"

enum {DIRICHLET = 1, NEUMANN = 2};
class Heat{
private:
	typedef double (*RHS)(double, double, double, double); 
	typedef double (*RHF)(double, double, double);
	typedef std::vector<std::vector<int> > INDEX;
    typedef std::vector<double> SOL;
	int begin, end, length, send_forward, send_backward;
	int recv_forward, recv_backward;
	
	INDEX ind;
	SOL sol;
	RHS f, u, g_up, g_down;
	RHF u0;
    int N, M;
	int size, rank;
	double h, CFL, t_end, t;
public:
	Heat(){};
	void set_size(int);
	void set_rank(int);
	void set_f(const RHS&);
	void set_Initial(const RHF&);
	void set_Boundrary(int, const RHS&);
	void set_N(int);
	void set_t(double);
	void set_CFL(double);
	// set the real solution, only for test
	void set_Solution(const RHS&);
	// calsulate the L2 norm error
	double error();
	std::vector<double> solve();
private:
	int transform(int, int, int);
	void prepare();
    void init();
	void onestep(double);
};


#endif
