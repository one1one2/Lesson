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
	fftw_plan p;

/*	std::vector<int> x_sendcnt, x_sdispls, x_recvcnt, x_rdispls;
	std::vector<int> y_sendcnt, y_sdispls, y_recvcnt, y_rdispls;
	std::vector<int> z_sendcnt, z_sdispls, z_recvcnt, z_rdispls;
	MPI_Datatype xvector, yvector, zvector;
	MPI_Datatype xvector_single, yvector_single, zvector_single;
	MPI_Datatype recvvector, recvvector_single; */
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
	void fft(MPI_Comm);
	void onestep(); 
};

#endif
