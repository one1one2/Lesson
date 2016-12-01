/**
 * @file hw3.h
 * @brief 
 * @author  lczheng, lczheng@pku.edu.cn 
 * 
 * @date 2016-11-27
 */

#ifndef __hw3_h__
#define __hw3_h__

#include<iostream>
#include<cmath>
#include<cstdlib>
#include<complex>
#include<fftw3.h>
#include<vector>
#include"mpi.h"

class hw3{
private:
	typedef double (*RHS)(double, double, double);
	typedef std::complex<double> COMPLEX;
	typedef std::vector<std::complex<double> > SOL;
	RHS f, u;
	int size, rank, Np;
	double x0, y0, z0;
	int Nx, Ny, Nz;
	int N, step, M;
	double h;
	SOL sol, F, soltemp;
	std::vector<fftw_plan> p;
	double lambda;
    MPI_Comm xcomm, ycomm, zcomm;
public:
	hw3(){};
	void set_size(int);
	void set_rank(int);
	void set_f(const RHS&);
	void set_N(int);
	void set_step(int);
	SOL solve();

	void set_solution(const RHS&);
	double error();

private:
	double u0(double, double, double);
	void initial();
	void fft_1d(MPI_Comm, int);
    void fft_3d(int);
	void onestep(); 
};

#endif
