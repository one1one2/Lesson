
# brief: first compile the .cpp file and generate an executive file \
then execute the .exe file to get result
# output:the error and calculation time,including a .txt fike of the solution u;
# author: liruotai
# date: 2017/11/19


# CXX编译器
CXX=mpicxx

# RANK 进程数
RANK=4

# N 网格剖分所需节点值，注意 横向节点数：4N+1,纵向节点数：2N+1，原因请见报告secondhw.pdf
N=100

# KMAX 最大迭代次数
KMAX=10000

# EXE执行homework.cpp编译后生成的homewrok文件
EXE=mpirun

all:result

homework:homework.cpp
	$(CXX) homework.cpp -o homework -O3

result:homework
	$(EXE) -n $(RANK) ./homework $(N) $(KMAX)



