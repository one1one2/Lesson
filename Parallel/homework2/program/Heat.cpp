/**
 * @file Heat.cpp
 * @brief 
 * @author  lczheng, lczheng@pku.edu.cn 
 * 
 * @date 2016-11-08
 */

#include "Heat.h"

/**
 * @brief 设置进程总数
 *
 * @param size1
 */
void Heat::set_size(int size1){
	size = size1;
}

/**
 * @brief 设置进程编号
 *
 * @param rank1
 */
void Heat::set_rank(int rank1){
	rank = rank1;
}

/**
 * @brief 设置右端项f
 *
 * @param fun
 */
void Heat::set_f(const RHS &fun){
	f = fun;
}

/**
 * @brief 设置初值
 *
 * @param fun
 */
void Heat::set_Initial(const RHF &fun){
	u0 = fun;
}

/**
 * @brief 设置边值
 *
 * @param flag  标志， DIRICHLET 或者 NEUMANN
 * @param fun
 */
void Heat::set_Boundrary(int flag , const RHS &fun){
	if (flag == DIRICHLET)
		g_up = fun;
	else if (flag == NEUMANN) 
		g_down = fun;
}

/**
 * @brief 设置网格密度
 *
 * @param N1
 */
void Heat::set_N(int N1){
	N = N1; 
	h = 1./(double)N;
	M = (N + 1) * (N + 2) * (2 * N + 3) / 6;
}


/**
 * @brief 设置计算终止时间
 *
 * @param t1
 */
void Heat::set_t(double t1){
	t_end = t1;
}

/**
 * @brief 设置CFL条件数
 *
 * @param CFL1
 */
void Heat::set_CFL(double CFL1){
	CFL = CFL1;
}

/**
 * @brief 设置真解，可以用来对有真解的情况测试误差.
 *
 * @param uu
 */
void Heat::set_Solution(const RHS &uu){
	u = uu;
}

/**
 * @brief 利用坐标取得点的编号的一个帮助函数。
 */
int Heat::transform(int i, int j, int k){
	if (i < 0 || j < 0 || i + j + k > N || i + j - k > N)
		return -1;
	if (k <= 0)
		return (N + k) * (N + k + 1) * (N + k + 2) / 6 
		  + j * (2 * N + 2 * k - j + 3) / 2 + i;  
	return  (N + 1) * (N + 2) * (2 * N + 3) / 6 
			-  (N - k + 1) * (N - k + 2) * (N - k + 3) / 6
			+ j * (2 * N - 2 * k - j + 3) / 2 + i;  
}

/**
 * @brief 由0号进程完成的一些初始化工作，包括给出所有编号的点的坐标，以及他们邻居的点
 * 的编号，并将这些信息发送给所需要的进程。
 */
void Heat::prepare(){ 
	ind.resize(9);
	for (int i = 0; i < 9; i++)
	  ind[i].resize(M);

	for (int i = 0; i <= N; i++)
	  for (int j = 0; j <= N - i; j++)
		for (int k = i + j - N; k <= N - i - j; k++){
			int index = transform(i, j, k);
			ind[0][index] = i;    
			ind[1][index] = j;    
			ind[2][index] = k;
			ind[3][index] = transform(i, j + 1, k);
			ind[4][index] = transform(i, j - 1, k);
			ind[5][index] = transform(i - 1, j, k);
			ind[6][index] = transform(i + 1, j, k);
			ind[7][index] = transform(i, j, k + 1);
			ind[8][index] = transform(i, j, k - 1);
		}
	for (int i = 1; i < size; i++)
	  for (int j = 0; j < 9; j++)
		MPI_Send(&ind[j][(i * M) / size], (int)(((i + 1) * M) / size) - (int)((i * M) / size), 
					MPI_INT, i, 0, MPI_COMM_WORLD);
}

/**
 * @brief 每个进程单独的初始化工作，包括计算需要向相邻进程索取的信息，
 * 和发送给相邻进程的信息，以及做t=0时刻解的初始化。
 */
void Heat::init(){ 
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
			if (ind[j][index] - end + 1 > recv_backward) 
			  recv_backward = ind[j][index] - end + 1;
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
	if (rank == 0) 
		sol.resize(M);
	else 
		sol.resize(length);
	for (int index = 0; index < length; index++){
		sol[index] = u0(ind[0][index] * h, ind[1][index] * h, ind[2][index] * h);
	}
	if (rank > 0)
		MPI_Send(&sol[0], send_forward, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
	if (rank < size - 1) 
		MPI_Send(&sol[length - send_backward], send_backward, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
}

/**
 * @brief 一步迭代
 *
 * @param dt 时间步长
 */
void Heat::onestep(double dt){
	SOL temp;
	temp.resize(recv_forward + length + recv_backward);
	if (rank > 0)
		MPI_Recv(&temp[0], recv_forward, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	if (rank < size - 1) 
		MPI_Recv(&temp[recv_forward + length], recv_backward, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	for (int index = 0; index < length; index++)
	  temp[recv_forward + index] = sol[index];

	for (int index = 0; index < length; index++){
		int i, j, k;
		i = ind[0][index];
		j = ind[1][index];
		k = ind[2][index];
		if (i + j + k == N){  //上边界——狄利克雷。 
			sol[index] = u(i * h, j * h, k * h, t + dt); 
		}else if (i + j - k != N){ // 下平面的点最后再处理。
			if (i == 0 && j == 0){ //对称性处理x=y=0的情况。
				sol[index] = temp[index + recv_forward] + dt * f(i * h, j * h, k * h, t)
					+ dt / h / h * (2 * temp[ind[3][index]] + 
					+ 2 * temp[ind[6][index]] + temp[ind[7][index]] + temp[ind[8][index]] 
					- 6 * temp[index + recv_forward]);
			}else if (i == 0){ //处理x=0。
				sol[index] = temp[index + recv_forward] + dt * f(i * h, j * h, k * h, t)
					+ dt / h / h * (temp[ind[3][index]] + temp[ind[4][index]] 
					+ 2 * temp[ind[6][index]] + temp[ind[7][index]] + temp[ind[8][index]] 
					- 6 * temp[index + recv_forward]);
			}else if (j == 0){ //处理y=0。
				sol[index] = temp[index + recv_forward] + dt * f(i * h, j * h, k * h, t)
					+ dt / h / h * (2 * temp[ind[3][index]] + temp[ind[5][index]] 
					+   temp[ind[6][index]] + temp[ind[7][index]] + temp[ind[8][index]] 
					- 6 * temp[index + recv_forward]);
			}else { // 内部的点，采用一步显式差分格式
				sol[index] = temp[index + recv_forward] + dt * f(i * h, j * h, k * h, t)
					+ dt / h / h * (temp[ind[3][index]] + temp[ind[4][index]] + temp[ind[5][index]] 
					+	temp[ind[6][index]] + temp[ind[7][index]] + temp[ind[8][index]] 
					- 6 * temp[index + recv_forward]);
			}
		}else {  //下边界——诺依曼
			if (i != 0 && j != 0){
				sol[index] = temp[index + recv_forward] + dt * f(i * h, j * h, k * h, t)
					+ dt / h / h * ( 2 * temp[ind[4][index]] + 2 * temp[ind[5][index]] 
					+ 2 * temp[ind[7][index]] + 2 * sqrt(3.0) * h * g_down(i * h, j * h, k * h, t)  
					- 6 * temp[index + recv_forward]); 
			} 
			if ( i == 0 && j != 0){
				sol[index] = temp[index + recv_forward] + dt * f(i * h, j * h, k * h, t)
					+ dt / h / h * ( + 2 * temp[ind[4][index]] 
					+ 2 * temp[ind[7][index]] + 2 * sqrt(3.0) * h * g_down(i * h, j * h, k * h, t)  
					- 4 * temp[index + recv_forward]);
			} 
			if (i != 0 && j == 0){
				sol[index] = temp[index + recv_forward] + dt * f(i * h, j * h, k * h, t)
					+ dt / h / h * ( 2 * temp[ind[5][index]] 
					+ 2 * temp[ind[7][index]] + 2 * sqrt(3.0) * h * g_down(i * h, j * h, k * h, t)  
					- 4 * temp[index + recv_forward]); 
			} 
			if (i == 0 && j == 0){
				sol[index] = temp[index + recv_forward] + dt * f(i * h, j * h, k * h, t)
					+ dt / h / h * ( 2 * temp[ind[7][index]] + 2 * sqrt(3.0) * h * g_down(i * h, j * h, k * h, t)  
					- 2 * temp[index + recv_forward]); 
			}
		}
	} 
	t += dt;
    if (t < t_end){
		if (rank > 0)
		 MPI_Send(&sol[0], send_forward, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
		if (rank < size - 1) 
		 MPI_Send(&sol[length - send_backward], send_backward, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
	} 
}


/**
 * @brief 求解函数
 *
 * @return 0号进程返回整个解，其余进程返回各自负责的区域的解。
 */
std::vector<double> Heat::solve(){
	if (rank == 0){
		std::cerr<<"Prepareing... "<<std::endl;
       prepare();
	}
	init(); 	
	double tau = CFL * h * h;
	if (rank == 0){
	  std::cerr<<"Calculating..."<<std::endl;
	}
	while (t + tau < t_end){
		onestep(tau);
	}
	if (t < t_end){
	  onestep(t_end - t);
	}
	if (rank != 0){
	  MPI_Send(&sol[0], length, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
    if (rank == 0){
		for (int i = 1; i < size; i++){
			MPI_Recv(&sol[(i * M) / size], (int)(((i + 1) * M) / size) - (int)((i * M) / size), 
					MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	return sol;
}

/**
 * @brief 计算误差的一个测试函数
 *
 * @return L2误差 
 */
double Heat::error(){
	double sum = 0;
	for (int i = 0; i <= N; i++)
		for (int j = 0; j <= N - i; j++)	
			for (int k = i + j - N; k <= N - i - j ; k++){
			  sum += (u(i*h,j*h,k*h,t_end) - sol[transform(i,j,k)]) * 
				  (u(i*h,j*h,k*h,t_end) - sol[transform(i,j,k)]);
			}
	return std::sqrt(sum/M);
}  
