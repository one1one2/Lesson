#include "hw3.h"

void hw3::set_size(int size1){
	size = size1;
	Np = (int)pow(size,1./3);
}

void hw3::set_rank(int rank1){
	rank = rank1;
}

void hw3::set_f(const RHS& f1){
	f = f1;
}

void hw3::set_N(int N1){
	N = N1;
	h = 2*M_PI/N;
}

void hw3::set_step(int step1){
	step = step1;
}

void hw3::set_solution(const RHS& u1){
	u = u1;
}



double hw3::u0(double x, double y, double z){
	return (x+100*y+10000*z);
}

void hw3::initial(){
	Nx = rank%Np;
	Ny = (rank/Np)%Np;
	Nz = rank/Np/Np;

	M = N/Np;
	
	MPI_Comm_split(MPI_COMM_WORLD, Nz*Np + Ny, Nx, &xcomm);
	MPI_Comm_split(MPI_COMM_WORLD, Nz*Np + Nx, Ny, &ycomm);
	MPI_Comm_split(MPI_COMM_WORLD, Ny*Np + Nx, Nz, &zcomm);
/*

	MPI_Type_vector(1, M, 1, MPI_DOUBLE_COMPLEX, &xvector_single);
	MPI_Type_hvector(M*M/Np, 1, sizeof(MPI_DOUBLE_COMPLEX)*M, xvector_single, &xvector);
	MPI_Type_commit(&xvector);
	MPI_Type_vector(M*M/Np, M, N, MPI_DOUBLE_COMPLEX, &recvvector_single);
	MPI_Type_hvector(1, 1, sizeof(MPI_DOUBLE_COMPLEX)*N, recvvector_single, &recvvector);
	MPI_Type_commit(&recvvector);
	
	MPI_Type_vector(M*M, M/Np, M, MPI_DOUBLE_COMPLEX, &yvector);
	MPI_Type_commit(&yvector);
	MPI_Type_vector(M*M, M/Np, M*M, MPI_DOUBLE_COMPLEX, &zvector);
	MPI_Type_commit(&zvector);   */
/*
	x_sendcnt.resize(Np,1);
	x_sdispls.resize(Np);
	x_recvcnt.resize(Np,M*M*M/Np);
	x_rdispls.resize(Np);

	x_sdispls[0] = 0;
	x_sdispls[1] = 1;
	x_rdispls[0] = 0;
	x_rdispls[1] = M*M*M/Np;
*/

	x0 = 2*M_PI/Np*Nx - M_PI;
	y0 = 2*M_PI/Np*Ny - M_PI;
	z0 = 2*M_PI/Np*Nz - M_PI;

	F.resize(M*M*M);   
	sol.resize(M*M*M);
	soltemp.resize(M*M*M);

	for (int i = 0; i < M*M*M; i++){
//		sol[i] = 0;
		F[i] = f(x0 + i%M*h, y0 + (i/M)%M*h, z0 + i/M/M*h);
		sol[i] = u0(x0 + i%M*h, y0 + (i/M)%M*h, z0 + i/M/M*h);
	}  


}

void hw3::fft(MPI_Comm comm){
	MPI_Alltoall((void*)&sol[0], M*M*M/Np, MPI_DOUBLE_COMPLEX, (void*)&soltemp[0], M*M*M/Np, MPI_DOUBLE_COMPLEX, comm);
	for (int i = 0; i < M; i++)
	  for (int j = 0; j < M; j++)
		for (int k = 0; k < M; k++)
		  sol[k*M*M + j*M + i] = soltemp[j*M*M + k*M + i];
	for (int i = 0; i < M*M/Np; i++){
		p = fftw_plan_dft_1d(N, 
				reinterpret_cast<fftw_complex*>(&sol[i*N]),
				reinterpret_cast<fftw_complex*>(&soltemp[i*N]),
				FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(p);
	}
	for (int i = 0; i < M; i++)
	  for (int j = 0; j < M; j++)
		for (int k = 0; k < M; k++)
		  sol[k*M*M + j*M + i] = soltemp[j*M*M + k*M + i];
	MPI_Alltoall((void*)&sol[0], M*M*M/Np, MPI_DOUBLE_COMPLEX, (void*)&soltemp[0], M*M*M/Np, MPI_DOUBLE_COMPLEX, comm); 
	sol = soltemp;
}


void hw3::onestep(){
	//计算 u^3-f
/*	for (int i = 0; i < M*M*M; i++){
		sol[i] = sol[i]*sol[i]*sol[i] - F[i];
	}  */ 
	//对u^3-f 作FFT
	fft(xcomm);  //x方向的FFT,包括数据传输和FFT
	for (int i = 0; i < M; i++)
	  for (int j = 0; j < M; j++)
		for (int k = 0; k < M; k++)
		  soltemp[k*M*M + j*M + i] = sol[k*M*M + i*M + j];
	sol = soltemp;
	fft(ycomm);
	for (int i = 0; i < M; i++)
	  for (int j = 0; j < M; j++)
		for (int k = 0; k < M; k++)
		  soltemp[k*M*M + j*M + i] = sol[k*M*M + i*M + j];
	sol = soltemp;
	for (int i = 0; i < M; i++)
	  for (int j = 0; j < M; j++)
		for (int k = 0; k < M; k++)
		  soltemp[k*M*M + j*M + i] = sol[i*M*M + j*M + k];
	sol = soltemp;
	fft(zcomm);
	for (int i = 0; i < M; i++)
	  for (int j = 0; j < M; j++)
		for (int k = 0; k < M; k++)
		  soltemp[k*M*M + j*M + i] = sol[i*M*M + j*M + k];
	sol = soltemp;
//	MPI_Alltoallv((void*)&sol[0], (int *)&x_sendcnt[0], (int *)&x_sdispls[0],  xvector, 
//				(void*)&soltemp[0], (int *)&x_recvcnt[0], (int *)&x_rdispls[0], MPI_DOUBLE_COMPLEX, xcomm);

	
	//计算u的FFT	
/*	for (int i = 0; i < M*M*M; i++){
		if (rank > 0 || i > 0){ 
			sol[i]/= (i%M + Nx)*(i%M + Nx) 
		 	+   ((i/M)%M + Ny)*((i/M)%M + Ny)
			+	(i/M/M + Nz)*(i/M/M + Nz);
		} else{
		}
	}  */
	//FFT逆变换得到u
}

std::vector<std::complex<double> > hw3::solve(){
	initial();
//	for (int stepiter = 0; stepiter < step; stepiter++){
	  onestep();
//	}
    
		  if (rank == 0)
	  for (int i = 0; i < M*M*M; i++){
			 std::cout<<sol[i]<<"  "<<soltemp[i]<<" "<<std::endl;
		}
//
/*	MPI_Type_free(&xvector);
	MPI_Type_free(&yvector);
	MPI_Type_free(&zvector); 
	MPI_Type_free(&recvvector);  
	fftw_destroy_plan(p);   */
	MPI_Comm_free(&xcomm);
	MPI_Comm_free(&ycomm);
	MPI_Comm_free(&zcomm);
	return sol;
}
