/**
* \file homework2.cpp
* \brief solve Possion equation in a hexagon with Neumman boundary condition,
* \author Ruotai Li, liruotai@yeah.net
* \version 1.0.0 
* \date 2017-11-19
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include<mpi.h>

#define eps 1e-6  ///迭代的误差精度；

using namespace std;  ///后面没有用到std域内定义的函数，不会引起冲突；
      
double Error(const vector<double>& r, const vector<double>& p, int begin, int end);///函数申明，最后有定义；

int main(int argc,char *argv[]){  
  int k,kmax;  ///k为网格剖分需要的整数，kmax最大迭代次数，是一个终止条件； 
  int rank,size; 
  double t1,t2; ///t1记录并行开始时间，t2记录程序结束时间
  const double c=-1.25; /// 计算得到的Neunmman边值条件常数；
  double f0=0.0, f1=1.0, f2=2.0, f3=3.0, f4=4.0, f5=5.0; ///方程右端项在区域内部为分片常系数，从0开始，逆时针方向将区域内部等分成6部分，fi=i，表示第i部分；
  double err;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if(rank==0){
    if(argc<3){
      cerr<<"Usage: mpirun -n <size> ./homework2 <k> <kmax>"<<endl;///如果输入数量小于3，程序会输出提示
      abort();
      return 0;
    }
    t1=MPI_Wtime();
  }
  k=atoi(argv[1]);///网格剖分情况：横向网格节点数：4*k+1,纵向网格节点数：2*k+1
  kmax=atoi(argv[2]);///最大迭代次数
  /**
   * 对正六边形局域进行矩形网格剖分，引入虚拟节点
   */
  int m = 4*k+1;  ///网格横向节点数
  int n = 2*k+1; /// 网格纵向节点数
  /**
   * 网格剖分后，共有（m*n)个节点，按从左下到右上顺序排列*
   */
  const double dx = sqrt(3.)/(3.*k); ///单位矩形网格的横向长度
  const double dy = 1./k;   ///单位矩形网格的纵向长度
  /**
   * 动态生成节点差分后的系数矩阵，维度（m*n,5）
   */

  vector<vector<double> > A(m*n,vector<double>(5));
  for(int i=0;i<m*n;i++)
  {
    for(int j=0;j<5;j++)
    {
      A[i][j]=0;
    }
  }
  /**
   * 生成右端项,与节点对应的列向量
   */
  vector<double> F(m*n,0);
  /**
   * 拼写系数矩阵用稀疏方式存储，以及右端项,节点标号以从下到上，从左到右排列 *
   */
  for(int i=0; i < n;i++)
  {
    for(int j=0;j<m;j++)
    {
      if( i==0 && (k<=j && j<=3*k)) ///判断点是否在下边界
      {  
        A[i*m+j][2]=1.0;
        A[i*m+j][4]=-1.0;
        F[i*m+j]=dy*c; 
      }
      if (i==2*k && (k<=j && j<=3*k)) ///判断点是否在上边界
      {  
        A[i*m+j][2]=1.0;
        A[i*m+j][0]=-1.0;
        F[i*m+j]=dy*c; 
      }
      if ((0<i && i<k) && (k-i<=j&&j<=3*k+i)) ///判断点是否在正六边形的下半部
      {    
        if(j==k-i) ///判断点是否在左下边界上 
        { 
          A[i*m+j][2]=4.0;
          A[i*m+j][3]=-3.0;
          A[i*m+j][4]=-1.0;
          F[i*m+j]=2*dy*c; }
        else if(j==3*k+i) ///判断点是否在右下边界上	  
        { 
          A[i*m+j][1]=-3.0;
          A[i*m+j][2]=4.0;
          A[i*m+j][4]=-1.0;
          F[i*m+j]=2*dy*c;}
        else {            /// 点在正六边形的下半部且为内点
          A[i*m+j][0]=-1.0;
          A[i*m+j][1]=-3.0;
          A[i*m+j][2]=8.0;
          A[i*m+j][3]=-3.0;
          A[i*m+j][4]=-1.0;
          if(k-i<j&&j<k+i) ///判断点是否属于第4部分
          F[i*m+j]=dy*dy*f3;
          else if (3*k-i<=j&&j<3*k+i) ///判断点是否属于第6部分
          F[i*m+j]=dy*dy*f5;
          else
          F[i*m+j]=dy*dy*f4; ///点在第5部分
        }
      } ///正六边形下半部点判断已经相应系数、右端项赋值完毕

      if (i==k) ///处理横向中心线上的点
      {            
        if(j==0) ///判断是否为左端点
        { 
          A[i*m+j][2]=3.0;
          A[i*m+j][3]=0.0;
          F[i*m+j]=dx*c;
        }
        else if (j==4*k) ///判断是否为右端点
        { 
          A[i*m+j][2]=1.0;
          A[i*m+j][1]=-1.0;
          F[i*m+j]=dx*c;
        }
        else{          
          A[i*m+j][0]=-1.0;
          A[i*m+j][1]=-3.0;
          A[i*m+j][2]=8.0;
          A[i*m+j][3]=-3.0;
          A[i*m+j][4]=-1.0;
          if(0<j&&j<=2*k)
          F[i*m+j]=dy*dy*f3;
          else 
          F[i*m+j]=dy*dy*f0; 
        }
      }
      if ((k<i && i<2*k) && (i-k<=j&&(j<=5*k-i))) ///判断点是否在正六边形的上半部
      {        if(j==i-k) ///判断点是否在左上边界上 
        { A[i*m+j][2]=4.0;
          A[i*m+j][3]=-3.0;
          A[i*m+j][0]=-1.0;
          F[i*m+j]=2*dy*c; }
        else if(j==5*k-i)///判断点是否在右上边界上
        {A[i*m+j][1]=-3.0;
          A[i*m+j][2]=4.0;
          A[i*m+j][0]=-1.0;
          F[i*m+j]=2*dy*c;}
        else {           /// 点在正六边形的上半部且为内点
          A[i*m+j][0]=-1.0;
          A[i*m+j][1]=-3.0;
          A[i*m+j][2]=8.0;
          A[i*m+j][3]=-3.0;
          A[i*m+j][4]=-1.0;
          if(i-k<j&&j<=3*k-i) ///判断点是否属于第3部分
          F[i*m+j]=dy*dy*f2;
          else if (3*k-i<j&&j<=k+i) ///判断点是否属于第2部分
          F[i*m+j]=dy*dy*f1;
          else
          F[i*m+j]=dy*dy*f0; ///点在第1部分
        }
      } ///正六边形上半部点判断已经相应系数、右端项赋值完毕

      /**
       * 节点对应的系数矩阵，以及右端项已经处理完毕 
       */

    }
  }

  /* 
       ofstream outfile1,outfile2;
       outfile1.open("xishu.txt");
       outfile2.open("youduan.txt");
       for(int i=0；i<m*n;i++)
       {      outfile2<<F[i]<<endl;
       for(int j=0;j<5;j++)
       {
       outfile1<<A[i][j]<<endl;
       }
       }
       outfile1.close();
       outfile2.close();
       */      ///输出右端项和系数矩阵，检查是否拼接错误
 
 /**
  * 生成节点列向量，方向从左下到右上排列
  */
  vector<double> U(m*n,0); ///生成初始向量
  vector<double> Uold(m*n,0);///存储上一步的值

  vector<int> begin(size,0); ///记录每个进程起始位置和终止位置的向量
  vector<int> end(size,0);
  for(int i = 0; i < size; i++)
  {
    begin[i]=m*n*i/size;
    end[i]=m*n*(i+1)/size;  ///每项进程， 进程初始位置为begin,终止位置为end-1
  }
  
  int N=0; ///计算迭代次数
  double rho=1.0;
  while((rho>eps)&&(N<kmax))  ///开始迭代
  {  
    /**
     * Jacob迭代
     */
    Uold = U;
    for(int i = begin[rank]; i < end[rank]; i++) ///对于每个进程独立计算部分迭代
    {   
      if (A[i][2]!=0)
      {   
        if (i < m) ///判断点是否在下边界上
        U[i]=(-(A[i][4]*Uold[i+m])+F[i])/A[i][2];///下边界点的迭代
        else if (i/m==2*k) ///判断点是否在上边界上
        U[i]=(-(A[i][0]*Uold[i-m])+F[i])/A[i][2];///上边界点的迭代
        else
        { 
          if (i==m*k) ///判断是否为左端点 
          U[i]=(-(A[i][3]*Uold[i+1])+F[i])/A[i][2];
          else if (i==(m*k+4*k))///判断是否为右端点
          U[i]=(-(A[i][1]*Uold[i-1])+F[i])/A[i][2];
          else 
          U[i]=(-(A[i][0]*Uold[i-m]+A[i][1]*Uold[i-1]+A[i][3]*Uold[i+1]
                  +A[i][4]*Uold[i+m])+F[i])/A[i][2];
        }
     }
    }  
    err = Error(U,Uold,begin[rank],end[rank]);///计算每个进程迭代前后两次误差的二范数
    MPI_Reduce(&err, &rho, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); ///将每个进程得到的误差归一到0号进程，得到前后两次迭代总误差
    if (rank > 0){
      MPI_Send(&U[begin[rank]], m, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
    }  ///若进程不是0号进程，向其前一个进程，传递从该进程起始位置开始，m个u的值；
    if (rank < size - 1){
      MPI_Send(&U[end[rank] -m ], m, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    }    ///若进程不是最后一号进程，向其后一个进程，传递从该进程末位置开始，前m个u的值；
    if (rank < size - 1){
      MPI_Recv(&U[begin[rank+1]], m, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }  ///若进程不是最后一号进程，接收其后一个进程传来的，从后一个进程起始位置起,m个u的值；  
    if (rank > 0){
      MPI_Recv(&U[end[rank-1] - m], m, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } ///若进程不是第0号进程，接收其前一个进程传来的，从前一个进程末尾起，前m个u的值

    N+=1;
    if(rank==0) {
      rho=rho/(m*n);
      rho = sqrt(rho);
      ///计算总误差的二范数
    }

    MPI_Bcast(&rho,1,MPI_DOUBLE,0,MPI_COMM_WORLD);///将计算的误差广播给其他进程，每个进程独立判断是否终止；

  }
  if(rank!=0){
      MPI_Send(&U[begin[rank]], end[rank] - begin[rank], 
            MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    } ///迭代结束后，其他所有进程向0号进程发送自己进程内计算的u的值;
  if(rank==0){
    for (int i = 1; i < size; i++){
        MPI_Recv(&U[begin[i]], end[i] - begin[i], MPI_DOUBLE, i,
              0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   
    }///接收其他进程传送的u值,得到完整的计算结果
    cout<<"L2 error is "<<rho<<endl;  /// 输出误差
    t2=MPI_Wtime();
    cout<<"Clock time is "<<t2-t1<<endl; ///输出运行时间 
  /**
   * 将计算结果u输出到名为result的文本文件，用于作图
   */  
   ofstream outfile;        
    outfile.open("result.txt"); 
    for(int i=0;i<m*n;i++)
    {
      outfile<<U[i]<<endl;
    }
    outfile.close();   
  }

  MPI_Finalize();
  return 0;
}

/* --------------------------------------------------------------------------*/
/**
* \brief Calculating the difference of two vector in 2-Norm
*
* \param r vector
* \param p vector 
* \param begin initial postion of the vector 
* \param end   end position of the vector
*
* \return  Error of the two vector
*/
/* ----------------------------------------------------------------------------*/
double Error(const vector<double>& r, const vector<double>& p, int begin, int end){
  /**
   * 计算 向量r与p的差的二范数函数，r与p维度要相同
   *
   * len 是向量r,p的维度
   */
  double s = 0.0;
  for(int i=begin;i<end;i++)
  {
    s=s+(r[i]-p[i])*(r[i]-p[i]);
  }
  return s;
}


