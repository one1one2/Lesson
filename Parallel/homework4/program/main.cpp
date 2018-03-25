/**
 * @file main.cpp
 * @brief 
 * @author  lczheng, lczheng@pku.edu.cn 
 * 
 * @date 2016-12-31
 */

#include"Heat.h"

double u(double x, double y, double z, double t){
	return sin((x*x+y*y+z*z)*t);
//    return (x*x+y*y+z*z);
}

double f(double x, double y, double z, double t){
	return cos((x*x+y*y+z*z)*t)*(x*x+y*y+z*z-6*t)
		+sin((x*x+y*y+z*z)*t)*4*t*t*(x*x+y*y+z*z); 
//    return (-6);
}

double u0(double x, double y, double z){
	return u(x,y,z,0);
}

double g_up(double x, double y, double z, double t){
	return u(x,y,z,t);
}

double g_down(double x, double y, double z, double t){
	return cos((x*x+y*y+z*z)*t)*(2*t/sqrt(3.0));
//    return 2.0/sqrt(3.0);
}

int transform(int i, int j, int k,int N){
	if (i < 0 || j < 0 || i + j + k > N || i + j - k > N)
		return -1;
	if (k <= 0)
		return (N + k) * (N + k + 1) * (N + k + 2) / 6 
		  + j * (2 * N + 2 * k - j + 3) / 2 + i;  
	return  (N + 1) * (N + 2) * (2 * N + 3) / 6 
			-  (N - k + 1) * (N - k + 2) * (N - k + 3) / 6
			+ j * (2 * N - 2 * k - j + 3) / 2 + i;  
}
int main(int argc, char *argv[]){
    int rank, size;
	int N, M;
	double CFL = 2;
	double t_end;
	double t1, t2, err, err_each;
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
	}
		t1 = MPI_Wtime();
	N = atoi(argv[1]);
	M = (N + 1) * (N + 2) * (2 * N + 3) / 6;
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
    err_each = Qu.error();
    MPI_Reduce(&err_each, &err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0){    
        err = sqrt(err/M);
        std::cerr<<"L2 error is "<<err<<std::endl;
		std::cerr<<"Finished!"<<std::endl;
		t2 = MPI_Wtime();
		std::cerr<<"clock time is "<<t2 -t1<<std::endl;
	}
	MPI_Finalize();
	return 0;

}
