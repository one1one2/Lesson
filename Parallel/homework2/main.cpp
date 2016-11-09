/**
 * @file main.cpp
 * @brief 
 * @author  lczheng, lczheng@pku.edu.cn 
 * 
 * @date 2016-11-08
 */

#include"Heat.h"


double u(double x, double y, double z, double t){
//	return sin((x*x+y*y+z*z)*t);
    return (x*x+y*y+z*z)*exp(-t);
}

double f(double x, double y, double z, double t){
//	return cos((x*x+y*y+z*z)*t)*(x*x+y*y+z*z-6*t)
//		+sin((x*x+y*y+z*z)*t)*4*t*t*(x*x+y*y+z*z); 
    return (-6-(x*x+y*y+z*z))*exp(-t);
}

double u0(double x, double y, double z){
	return u(x,y,z,0);
}

double g_up(double x, double y, double z, double t){
	return u(x,y,z,t);
}

double g_down(double x, double y, double z, double t){
//	return cos((x*x+y*y+z*z)*t)*(2*t/sqrt(3.0));
    return 2*exp(-t)/sqrt(3.0);
}

int main(int argc, char *argv[]){
    int rank, size;
	int N;
	double CFL = 0.1;
	double t_end;
	double t1, t2;
	std::vector<double> sol;
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (rank == 0){	
		if (argc < 3){
			std::cerr<<"Usage: mpirun -n <size> ./main <N> <t_end>"<<std::endl;
			std::abort(); 
			return 0;
		}
		t1 = MPI_Wtime();
	}
	N = atoi(argv[1]);
	t_end = atof(argv[2]);
	Heat Qu;
	Qu.set_size(size);
	Qu.set_rank(rank);
	Qu.set_f(f);
	Qu.set_Initial(u0);
	Qu.set_Boundrary(DIRICHLET, g_up);
	Qu.set_Boundrary(NEUMANN, g_down);
	Qu.set_CFL(CFL);
	Qu.set_N(N);
	Qu.set_t(t_end);
	Qu.set_Solution(u);
	sol = Qu.solve();
    if (rank == 0){	
		std::cerr<<"L2 error is "<<Qu.error()<<std::endl;
		std::cerr<<"Finished!"<<std::endl;
		t2 = MPI_Wtime();
		std::cerr<<"clock time is "<<t2 -t1<<std::endl;
	}
	MPI_Finalize();
	return 0;

}
