/**
 * @file Heat.cpp
 * @brief 
 * @author  lczheng, lczheng@pku.edu.cn 
 * 
 * @date 2016-11-08
 */

#include "Heat.h"

double u(double x, double y, double z, double t){
	return t;
}
double f(double x, double y, double z, double t){
	return 1;
}
double u0(double x, double y, double z){
	return u(x,y,z,0);
}

int transform(int i, int j, int k, int N){
	if (i < 0 || j < 0 || i + j + k > N || i + j - k > N)
	  return -1;
	if (k <= 0)
	  return (N + k) * (N + k + 1) * (N + k + 2) / 6 
		  + j * (2 * N + 2 * k - j + 3) / 2 + i;  
	return  (N + 1) * (N + 2) * (2 * N + 3) / 6 
		-  (N - k + 1) * (N - k + 2) * (N - k + 3) / 6
		+ j * (2 * N - 2 * k - j + 3) / 2 + i;  
}

void prepare(int rank, int size, int N, int M, INDEX &ind){
	ind.resize(9);
	for (int i = 0; i < 9; i++)
	  ind[i].resize(M);

	for (int i = 0; i <= N; i++)
	  for (int j = 0; j <= N - i; j++)
		for (int k = i + j - N; k <= N - i - j; k++){
			int index = transform(i, j, k, N);
			ind[0][index] = i;    
			ind[1][index] = j;    
			ind[2][index] = k;
			ind[3][index] = transform(i, j + 1, k, N);
			ind[4][index] = transform(i, j - 1, k, N);
			ind[5][index] = transform(i - 1, j, k, N);
			ind[6][index] = transform(i + 1, j, k, N);
			ind[7][index] = transform(i, j, k + 1, N);
			ind[8][index] = transform(i, j, k - 1, N);
		}
	for (int i = 1; i < size; i++)
	  for (int j = 0; j < 9; j++)
		MPI_Send(&ind[j][(i * M) / size], (int)(((i + 1) * M) / size) - (int)((i * M) / size), 
					MPI_INT, i, 0, MPI_COMM_WORLD);
}

void init(int rank, int size, int N, double h, int M, int &begin, int &end, int &length, int &recv_forward,
			int &recv_backward, int &send_forward, int &send_backward, SOL &sol, INDEX &ind){ 
	begin = (rank * M) / size; 
	end = ((rank + 1) * M) / size;
	length = end - begin;
	if (rank != 0){
		 ind.resize(9);
		 for (int j = 0; j < 9; j++){
			ind[j].resize(length);
			MPI_Recv(&ind[j][0], length, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	recv_forward = 0;
	recv_backward = 0;
	for (int j = 3; j < 9; j++)
		for (int index = 0; index < length; index++){
			if (begin - ind[j][index] > recv_forward&& ind[j][index] != -1) 
			  recv_forward = begin - ind[j][index];
			if (ind[j][index] - end > recv_backward) 
			  recv_backward = ind[j][index] - end;
		}

	for (int j = 3; j < 9; j++)
		for (int index = 0; index < length; index++)
			ind[j][index] -= begin - recv_forward;
    if (rank > 0){
		MPI_Send(&recv_forward, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
	}
    if (rank < size - 1){
		MPI_Send(&recv_backward, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
	}
    if (rank > 0){
		MPI_Recv(&send_forward, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}else{ 
		send_forward = 0;
		recv_forward = 0;
	}
    if (rank < size - 1){
		MPI_Recv(&send_backward, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}else{ 
		send_backward = 0;
		recv_backward = 0;
	}
	for (int j = 3; j < 9; j++)
		for (int index = 0; index < length; index++)
			ind[j][index] -= begin + recv_forward;

	if (rank == 0) 
	  sol.resize(M);
	else 
	  sol.resize(length);
	for (int index = 0; index < length; index++){
		int i, j, k;
		i = ind[0][index]; j = ind[1][index]; k = ind[2][index];
		sol[index] = u0(i * h, j * h, k * h);
	}
}

void onestep(double dt, double &t, int rank, int size, int N, double h, int M, int begin, int end, int length, int recv_forward,
			int recv_backward, int send_forward, int send_backward, SOL &sol, INDEX ind, double t_end){	
	SOL temp;
	temp.resize(recv_forward + length + recv_backward);
	if (t > 0){
		if (rank > 0)
		 MPI_Recv(&temp[0], recv_forward, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (rank < size - 1) 
		 MPI_Recv(&temp[recv_forward + length], recv_backward, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	for (int index = 0; index < length; index++)
	  temp[recv_forward + index] = sol[index];
	for (int index = 0; index < length; index++){
		int i, j, k;
		i = ind[0][index];
		j = ind[1][index];
		k = ind[2][index];
		if (i + j + k == N) 
			sol[index] = u(i * h, j * h, k * h, t);
		else if (i == 0 || j == 0)
			sol[index] = u(i * h, j * h, k * h, t);
		else if (i + j - k == N) 
			sol[index] = u(i * h, j * h, k * h, t);
		else 
			sol[index] = temp[index + recv_forward] + dt * f(i * h, j * h, k * h, t)
					+ dt / h / h * (temp[ind[3][index]] + temp[ind[4][index]] + temp[ind[5][index]] 
					+	temp[ind[6][index]] + temp[ind[7][index]] + temp[ind[8][index]] 
					- 6 * temp[index + recv_forward]);	
	} 
	t += dt;
    if (t < t_end){
		if (rank > 0)
		 MPI_Send(&sol[0], send_forward, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
		if (rank < size - 1) 
		 MPI_Send(&sol[length - send_backward], send_backward, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
	} 
}

std::vector<double> solve(double CFL, double t_end, int rank, int size, int N){
	SOL sol;
	double h = 1.0 / N;
	int M = (N + 1) * (N + 2) * (2 * N + 3) / 6;
    int begin, end, length, recv_forward, recv_backward, send_forward, send_backward;
    INDEX ind;
	if (rank == 0){
		std::cout<<"Prepareing... form rank "<<rank <<std::endl;
       prepare(rank, size, N, M, ind);
		std::cout<<"Prepareing finished form rank "<<rank <<std::endl;
	}
	std::cout<<"Initializing... form rank "<<rank <<std::endl;
	init(rank, size, N, h, M, begin, end, length, recv_forward,
			recv_backward, send_forward, send_backward, sol, ind); 
	std::cout<<"Initializing finished form rank "<<rank <<std::endl<<std::flush;
	double t = 0;
	double tau = CFL * h * h;
    std::cerr<<"Calculating... from rank "<<rank<<std::endl;
	do{
		onestep(tau, t, rank, size, N, h, M, begin, end, length, recv_forward, recv_backward, send_forward, send_backward, sol, ind, t_end);
	}while ( t + tau < t_end);
	onestep(t_end - t, t, rank, size, N, h, M, begin, end, length, recv_forward, recv_backward, send_forward, send_backward, sol, ind, t_end);
    std::cerr<<"Calculation finished from rank "<<rank<<std::endl;
	if (rank != 0)
	  MPI_Send(&sol[0], length, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    if (rank == 0)
		for (int i = 1; i < size - 1; i++){
			MPI_Recv(&sol[(i * M) / size], (int)(((i + 1) * M) / size) - (int)((i * M) / size), 
					MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	return sol;
}


double error(SOL sol, int N, double t_end){
	double h = 1.0 / N;
	double sum = 0;
	for (int i = 0; i <= N; i++)
		for (int j = 0; j <= N - i; j++)	
			for (int k = i + j - N; k <= N - i - j ; k++)
			  sum += (u(i*h,j*h,k*h,t_end) - sol[transform(i,j,k,N)]) * 
				  (u(i*h,j*h,k*h,t_end) - sol[transform(i,j,k,N)]);
	return std::sqrt(sum);
}  
