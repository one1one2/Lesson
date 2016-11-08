/**
 * @file Heat.cpp
 * @brief 
 * @author  lczheng, lczheng@pku.edu.cn 
 * 
 * @date 2016-11-08
 */

#include "Heat.h"

void Heat::set_f(const RHS &fun){
	f = fun;
}

void Heat::set_Initial(const RHF &fun){
	u0 = fun;
}

void Heat::set_Boundrary(int flag , const RHS &fun){
	if (flag == DIRICHLET)
		u = fun;
	else if (flag == NEUMANN) 
		g = fun;
}

void Heat::set_N(){
	N = N1; 
	h = 1./(double)N;
}

Point Heat::transform(int N){


}

void Heat::init(){
	sol.resize(N * N * N);
	for (int i = 0; i < N ; i++)
}

void Heat::onestep(){
	sol += dt*
}
