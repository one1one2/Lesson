/**
 * @file hw3.h
 * @brief 
 * @author  lczheng, lczheng@pku.edu.cn 
 * 
 * @date 2016-11-27
 */

#ifndef __hw3_h__
#define __hw3_h__

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <fftw3.h>
#include <vector>
#include "mpi.h"

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
/**
 * @brief 设置进程数目
 *
 * @param  进程数目
 */
	void set_size(int);
	/**
	 * @brief 设置进程编号 
	 *
	 * @param int 进程编号
	 */
	void set_rank(int);
	/**
	 * @brief 设置右端项
	 *
	 * @param RHS 右端项函数
	 */
	void set_f(const RHS&);
	/**
	 * @brief 设置网格密度 
	 *
	 * @param int 网格密度
	 */
	void set_N(int);
	/**
	 * @brief 设置lambda参数
	 *
	 * @param double lambda 
	 */
	void set_lambda(double);
	/**
	 * @brief 设置迭代步数
	 *
	 * @param int 迭代步数
 	 */
	void set_step(int);
	/**
	 * @brief 求解
	 *
	 * @return 解
	 */
	SOL solve();

	/**
	 * @brief 设置真解
	 *
	 * @param RHS 真解
	 */
	void set_solution(const RHS&);
	/**
	 * @brief 设置误差
	 *
	 * @return 误差 
	 */
	double error();

private:
	double u0(double, double, double);
	void initial();
	void fft_1d(MPI_Comm, int);
    void fft_3d(int);
	void onestep(); 
};

#endif
