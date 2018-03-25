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
			ind[3][index] = transform(i, j, k - 1);
			ind[4][index] = transform(i, j - 1, k);
			ind[5][index] = transform(i - 1, j, k);
			ind[6][index] = transform(i + 1, j, k);
			ind[7][index] = transform(i, j + 1, k);
			ind[8][index] = transform(i, j, k + 1);
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
	end = ((rank + 1) * M) / size - 1;
	length = end - begin + 1;

	if (rank != 0){
		 ind.resize(9);
		 for (int j = 0; j < 9; j++){
			ind[j].resize(length);
			MPI_Recv(&ind[j][0], length, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, begin, end, begin, end, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A);
 
	int nnz;
	double values[7];
	int cols[7];

	for (int index = begin; index <= end; index++) {
		nnz = 0;
		int i = ind[0][index - begin], 
			j = ind[1][index - begin], 
			k = ind[2][index - begin];

		if (i + j + k == N) {  
			cols[nnz] = index;
			values[nnz] = 1.0;
			nnz++;
		} else { 
		cols[nnz] = index;
			if (i + j - k == N){
				values[nnz] = 6.0*CFL;
				if (i == 0) values[nnz] -= 2.0*CFL;
				if (j == 0) values[nnz] -= 2.0*CFL;
			} else {  
				values[nnz] = 1.0 + 6.0*CFL;
			}
			nnz++;
			for (int l = 3; l <= 8; l++) {
				if (ind[l][index - begin] >= 0) {
					cols[nnz] = ind[l][index - begin];
					if (ind[11 - l][index - begin] == -1)
						values[nnz] = -2.0*CFL;
					else 
						values[nnz] = -1.0*CFL;
					nnz++;
				}
			}
		}
		HYPRE_IJMatrixSetValues(A, 1, &nnz, &index, cols, values);
	}
	
	HYPRE_IJMatrixAssemble(A);
	HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);

    HYPRE_IJMatrixPrint(A, "IJ.out.A");
	HYPRE_IJVectorCreate(MPI_COMM_WORLD, begin, end, &b);
	HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, begin, end, &x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);


	/* Create solver */
	HYPRE_BoomerAMGCreate(&solver);

      /* Set some parameters (See Reference Manual for more parameters) */
      HYPRE_BoomerAMGSetPrintLevel(solver, 3);  /* print solve info + parameters */
      HYPRE_BoomerAMGSetCoarsenType(solver, 6); /* Falgout coarsening */
      HYPRE_BoomerAMGSetRelaxType(solver, 3);   /* G-S/Jacobi hybrid relaxation */
      HYPRE_BoomerAMGSetNumSweeps(solver, 1);   /* Sweeeps on each level */
      HYPRE_BoomerAMGSetMaxLevels(solver, 20);  /* maximum number of levels */
      HYPRE_BoomerAMGSetTol(solver, 1e-7);      /* conv. tolerance */

			
	rhs_values =  (double*) calloc(length, sizeof(double));
	x_values =  (double*) calloc(length, sizeof(double));
	rows = (int*) calloc(length, sizeof(int));
}


/**
 * @brief 一步迭代
 *
 * @param dt 时间步长
 */
void Heat::onestep(double dt){

	if (t > 0) HYPRE_IJVectorGetValues(x, length, rows, rhs_values);

	HYPRE_IJVectorInitialize(b);
    HYPRE_IJVectorInitialize(x);

	for (int index = 0; index < length; index++) {
		int i = ind[0][index], j = ind[1][index], k = ind[2][index];
		if (t < dt / 2)	{
			rhs_values[index] = u0(i*h, j*h, k*h);
		}
		if (i + j + k == N) {
			rhs_values[index] = g_up(i*h, j*h, k*h, t + dt);
		} else if (i + j - k == N) {
			rhs_values[index] = g_down(i*h, j*h, k*h, t + dt) * h * 2 * sqrt(3.0) * CFL;
		} else {
			rhs_values[index] += f(i*h, j*h, k*h, t + dt) * dt;
		}
		x_values[index] = 0.0;
		rows[index] = begin + index;
	}

	HYPRE_IJVectorSetValues(b, length, rows, rhs_values);
	HYPRE_IJVectorSetValues(x, length, rows, x_values);

	HYPRE_IJVectorAssemble(b);
	HYPRE_IJVectorGetObject(b, (void **) &par_b);

	HYPRE_IJVectorAssemble(x);
	HYPRE_IJVectorGetObject(x, (void **) &par_x);


      int num_iterations;
      double final_res_norm;


      /* Now setup and solve! */
      HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
      HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);

      /* Run info - needed logging turned on */
      HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
      HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
      if (rank == 0)
      {
         printf("\n");
         printf("Iterations = %d\n", num_iterations);
         printf("Final Relative Residual Norm = %e\n", final_res_norm);
         printf("\n");
      }


	t += dt;
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
	while (t < t_end){
		onestep(tau);
		if (rank == 0) std::cout<<"t = "<<t<<std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	sol.resize(length);
	HYPRE_IJVectorGetValues(x, length, rows, &sol[0]);
	/* Clean up */
    HYPRE_BoomerAMGDestroy(solver);
	free(x_values);
	free(rhs_values);
	free(rows);
	HYPRE_IJMatrixDestroy(A);
	HYPRE_IJVectorDestroy(b);
	HYPRE_IJVectorDestroy(x);

	return sol;
}

/**
 * @brief 计算误差的一个测试函数
 *
 * @return L2误差 
 */
double Heat::error(){
	double sum = 0;
	for (int index = 0; index < length; index++){
		int i = ind[0][index], 
			j = ind[1][index],
			k = ind[2][index];
//		if (fabs(u(i*h,j*h,k*h,t_end) - sol[index])>1e-5)
//			  std::cout<<i<<" "<<j<<" "<<k<<" "<<
//				  u(i*h,j*h,k*h,t_end)<<"  "<< sol[index]<<std::endl;
			sum += (u(i*h,j*h,k*h,t_end) - sol[index]) * 
				(u(i*h,j*h,k*h,t_end) - sol[index]);
		}
	return sum;
}  
