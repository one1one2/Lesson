/**
 * @file main.cpp
 * @brief 
 * @author  lczheng, lczheng@pku.edu.cn 
 * 
 * @date 2016-11-18
 */

#include "hw3.h"

double f(double x, double y, double z){
	return 0;
}
double u(double x, double y, double z){
	return 0;
}
int main(int argc, char **argv){
    int rank, size;
	int N;
	std::vector<std::complex<double> > sol;
	int step = 1;
	double t1, t2;
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (rank == 0){	
		if (argc < 2){
			std::cerr<<"Usage: mpirun -n <size> ./main <N>"<<std::endl;
			std::abort(); 
			return 0;
		}
		t1 = MPI_Wtime();
	}
	N = atoi(argv[1]);
	hw3 Qu;
	Qu.set_size(size);
	Qu.set_rank(rank);
	Qu.set_f(f);
	Qu.set_step(step);
	Qu.set_N(N);
	Qu.set_solution(u);
	sol = Qu.solve();
    if (rank == 0){	
//		std::cerr<<"L2 error is "<<Qu.error()<<std::endl;
		std::cerr<<"Finished!"<<std::endl;
		t2 = MPI_Wtime();
		std::cerr<<"clock time is "<<t2 -t1<<std::endl;
	}
	MPI_Finalize();
/*	typedef std::complex<double> COMPLEX; 
	int N = 10;
	std::vector<COMPLEX> in(N),out(N);
	*/
	return 0;
}

