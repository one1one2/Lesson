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
  // return u(x,y,z);
	return 0;
}

void hw3::initial(){
	Nx = rank%Np;
	Ny = (rank/Np)%Np;
	Nz = rank/Np/Np;

	M = N/Np;
	
	MPI_Comm_split(MPI_COMM_WORLD, Nz*Np + Ny, Nx, &xcomm);
	MPI_Comm_split(MPI_COMM_WORLD, Nz*Np + Nx, Ny, &ycomm);
	MPI_Comm_split(MPI_COMM_WORLD, Ny*Np + Nx, Nz, &zcomm);

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
	lambda = 10;
	p.resize(2*M*M/Np);
	for (int i = 0; i < M*M/Np; i++){
		p[i] = fftw_plan_dft_1d(N, 
				reinterpret_cast<fftw_complex*>(&sol[i*N]),
				reinterpret_cast<fftw_complex*>(&soltemp[i*N]),
				FFTW_FORWARD, FFTW_ESTIMATE);
		p[M*M/Np+i] = fftw_plan_dft_1d(N, 
				reinterpret_cast<fftw_complex*>(&sol[i*N]),
				reinterpret_cast<fftw_complex*>(&soltemp[i*N]),
				FFTW_BACKWARD, FFTW_ESTIMATE);
	}

}

void hw3::fft_1d(MPI_Comm comm, int fft_flag){
	MPI_Alltoall((void*)&sol[0], M*M*M/Np, MPI_DOUBLE_COMPLEX, (void*)&soltemp[0], M*M*M/Np, MPI_DOUBLE_COMPLEX, comm);
	for (int k = 0; k < Np; k++)
	  for (int j = 0; j < M*M/Np; j++)
		for (int i = 0; i < M; i++)
		  sol[j*N + k*M + i] = soltemp[k*M*M*M/Np + j*M + i];
	for (int i = 0; i < M*M/Np; i++){
		if (fft_flag == FFTW_FORWARD) 
			fftw_execute(p[i]);
		else 
			fftw_execute(p[i+M*M/Np]);
	}
	for (int k = 0; k < Np; k++)
	  for (int j = 0; j < M*M/Np; j++)
		for (int i = 0; i < M; i++)
		  sol[k*M*M*M/Np + j*M + i] = soltemp[j*N + k*M + i];
	MPI_Alltoall((void*)&sol[0], M*M*M/Np, MPI_DOUBLE_COMPLEX, (void*)&soltemp[0], M*M*M/Np, MPI_DOUBLE_COMPLEX, comm); 
}

void hw3::fft_3d(int fft_flag){
	fft_1d(xcomm, fft_flag);  //x方向的FFT,包括数据传输和FFT
	for (int i = 0; i < M; i++)
	  for (int j = 0; j < M; j++)
		for (int k = 0; k < M; k++)
		  sol[k*M*M + j*M + i] = soltemp[k*M*M + i*M + j];
	fft_1d(ycomm, fft_flag);
	for (int i = 0; i < M; i++)
	  for (int j = 0; j < M; j++)
		for (int k = 0; k < M; k++)
		  sol[k*M*M + j*M + i] = soltemp[i*M*M + k*M + j];
	fft_1d(zcomm, fft_flag);
	for (int i = 0; i < M; i++)
	  for (int j = 0; j < M; j++)
		for (int k = 0; k < M; k++)
			if (fft_flag == FFTW_FORWARD)
				sol[k*M*M + j*M + i] = soltemp[i*M*M + j*M + k];
			else 
				sol[k*M*M + j*M + i] = std::pow(1./N,3)*soltemp[i*M*M + j*M + k];
	
}

void hw3::onestep(){
	//计算 u^3-f

	for (int i = 0; i < M*M*M; i++){
		sol[i] = F[i] - sol[i]*sol[i]*sol[i] + lambda * sol[i];
	}   
//	if (rank == 0){
//		for (int i = 0; i < M*M*M; i++)
///		  std::cout<<i<<" "<<sol[i]<<std::endl;
//	}  
 	//对u^3-f 作FFT
	fft_3d(FFTW_FORWARD);
//	if (rank == 0){
//		for (int i = 0; i < M*M*M; i++)
//		  std::cout<<i<<" "<<sol[i]<<std::endl;
//	}  

	//计算u的FFT 除以 k^2+l^2+m^2
	for (int i = 0; i < M*M*M; i++){
		int coff = lambda + (N - abs(2*(i%M + Nx*N/Np) - N ) ) * (N - abs(2*(i%M + Nx*N/Np) - N)) / 4 
		 	+ (N - abs(2*((i/M)%M + Ny*N/Np) - N) )  * (N - abs(2*((i/M)%M + Ny*N/Np) - N) ) / 4  
			+ (N - abs(2*(i/M/M + Nz*N/Np) - N) ) * (N - abs(2*(i/M/M + Nz*N/Np) - N)) / 4;
		if (coff  > 0.5) 
			sol[i] = 1./coff * sol[i];
		else{
			std::cout<<"ERROR!";
			sol[i] = 0;	
		}
	}       
	//FFT逆变换得到u
	fft_3d(FFTW_BACKWARD);
}

std::vector<std::complex<double> > hw3::solve(){
	initial();
	for (int stepiter = 0; stepiter < step; stepiter++){
	  onestep();
	  /*
	double error_each = error();
	double err = 0;
	MPI_Reduce(&error_each, &err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);  
	err =sqrt(err/N/N/N);
    if (rank == 0){	
		std::cerr<<"step = "<<stepiter<<" L2 error is "<<err<<std::endl;
	}  */
	}
   /* 
		  if (rank == 0)
	  for (int i = 0; i < M*M*M; i++){
			 std::cout<<sol[i]<<"  "<<soltemp[i]<<" "<<std::endl;
		}
	*/
	for (int i=0; i < 2*M*M/Np; i++)
		fftw_destroy_plan(p[i]);   
	MPI_Comm_free(&xcomm);
	MPI_Comm_free(&ycomm);
	MPI_Comm_free(&zcomm);
	return sol;
}


double hw3::error(){
	double sum = 0;
	for (int i = 0; i < M*M*M; i++){
	//  if  (pow(sol[i].real() - u(x0 + i%M*h, y0 + (i/M)%M*h, z0 + i/M/M*h) , 2) + pow(sol[i].imag(),2)> 1e-10 )  
	//	std::cout<<"rank = "<<rank<<", i = "<<i<<", error = "<<fabs(sol[i].real() - u(x0 + i%M*h, y0 + (i/M)%M*h, z0 + i/M/M*h))<<"  "<<std::endl; 
	  sum += pow(sol[i].real() - u(x0 + i%M*h, y0 + (i/M)%M*h, z0 + i/M/M*h) , 2) 
	  + pow(sol[i].imag(),2) ;
	}
	return sum;
}
