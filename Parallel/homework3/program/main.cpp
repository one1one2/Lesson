/**
 * @file main.cpp
 * @brief 
 * @author  lczheng, lczheng@pku.edu.cn 
 * 
 * @date 2016-11-18
 */

#include "hw3.h"

double u(double x, double y, double z){
	return sin(x+1);
	return exp(exp(sin(x)));
//	return sin(x)*sin(y)*sin(z);
}
double f(double x, double y, double z){
	return sin(x+1)+pow(u(x,y,z),3.0);
      return 
		 exp(exp(sin(x))+sin(x))*(sin(x)-cos(x)*cos(x)
					  *(1+exp(sin(x))))
		  + pow (u(x,y,z),3.0);
//	return -exp(sin(x))*(cos(x)*cos(x)-sin(x))+pow(u(x,y,z),3.0);
//	return 3*sin(x)*sin(y)*sin(z)+pow(u(x,y,z),3);
}
int main(int argc, char **argv){
	MPI_Init(&argc, &argv);
    int rank, size;
	int N, Np;
	std::vector<std::complex<double> > sol;
	int step;
	double t1, t2;
	double error_each, error = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (rank == 0){	
		if (argc < 3){
			std::cerr<<"Usage: mpirun -n <size> ./main <N> <step>"<<std::endl;
			MPI_Abort(MPI_COMM_WORLD,0); 
			return 0;
		}
	}
	N = atoi(argv[1]);
	step = atoi(argv[2]);
	if (rank == 0){
		Np = (int)pow(size, 1./3);
		if (Np*Np*Np != size || N%(Np*Np) !=0){
			std::cerr<<"请确保进程数是立方数 Np^3 且网格密度N的2*Np^2的倍数。"<<std::endl;
			MPI_Abort(MPI_COMM_WORLD,1); 
		}
	}

	hw3 Qu;
	Qu.set_size(size);
	Qu.set_rank(rank);
	Qu.set_f(f);
	Qu.set_step(step);
	Qu.set_N(N);
	Qu.set_solution(u);

	t1 = MPI_Wtime();
	sol = Qu.solve();
	t2 = MPI_Wtime();

	error_each = Qu.error();
//	std::cout<<"rank = "<<rank<<" error= "<<error_each<<std::endl;
	MPI_Reduce(&error_each, &error, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	error =sqrt(error/N/N/N);
    if (rank == 0){	
		std::cerr<<"L2 error is "<<error<<std::endl;
		std::cerr<<"Finished!"<<std::endl;
		std::cerr<<"clock time is "<<t2 -t1<<std::endl;
	}  
	MPI_Finalize();
/*	typedef std::complex<double> COMPLEX; 
	int N = 10;
	std::vector<COMPLEX> in(N),out(N);
	*/
	return 0;
}

