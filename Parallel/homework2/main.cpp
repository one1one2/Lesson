/**
 * @file main.cpp
 * @brief 
 * @author  lczheng, lczheng@pku.edu.cn 
 * 
 * @date 2016-11-08
 */

#include"Heat.h"



int main(int argc, char *argv[]){
	int begin, end, length, send_forward, send_backward;
	int recv_forward, recv_backward, ruler;

	std::vector<std::vector<int> > ind;
    int N = 10;
	int size, rank;
	double CFL = 0.1, t_end = 1;
	SOL sol;
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    sol = solve(0.1, 0.1, rank, size, 10);
	std::cout<<"rank: "<<rank<<std::endl;
    if (rank == 0) 
		std::cout<<"L2 error is "<<error(sol , N,t_end)<<std::endl;
	std::cerr<<"Finished!"<<std::endl;
	MPI_Finalize();
	return 0;

}
