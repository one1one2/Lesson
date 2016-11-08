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

typedef std::vector<std::vector<int> > INDEX;
typedef std::vector<double> SOL;

double u(double x, double y, double z, double t);
double f(double x, double y, double z, double t);
double u0(double x, double y, double z);
int transform(int, int, int, int);
void prepare();
void init();
void init(int rank, int size, int N, double h, int M, int &begin, int &end, int &length, int &recv_forward,
			int &recv_backward, int &send_forward, int &send_backward, SOL &sol, INDEX &ind); 
void onestep(double dt, double &t, int rank, int size, int N, double h, int M, int begin, int end, int length, int recv_forward,
			int recv_backward, int send_forward, int send_backward, SOL &sol, INDEX ind, double t_end);  
SOL solve(double CFL, double t_end, int rank, int size, int N);

double error(SOL sol, int N, double t_end);


#endif
