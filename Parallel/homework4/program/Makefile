CC = mpicxx
OBJ = Heat.o 
LIBS = -L/home/one1one2/Software/hypre/lib/ -lHYPRE
CXXINCLUDES = -I/home/one1one2/Software/hypre/include/
LD_FLAG = -lm -ldl -Wall  -O3 $(LIBS) $(CXXINCLUDES)

default : main

$(OBJ) : %.o : %.cpp %.h   
	$(CC) -c $< -o $@ $(LD_FLAG) 

main.o : main.cpp *.cpp *.h  
	$(CC) -c -o main.o main.cpp $(LD_FLAG) 

main : $(OBJ) main.o 
	$(CC) -o main $(OBJ) main.o $(LD_FLAG) 

clean :
	rm -f $(OBJ) main.o
	rm -f main

.PHONY : clean
