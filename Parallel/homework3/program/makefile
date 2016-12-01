CC = mpicxx
OBJ = hw3.o

default : main

$(OBJ) : %.o : %.cpp %.h
	$(CC) -c $< -o $@ -lm -ldl -Wall -lfftw3 -O2 

main.o : main.cpp *.cpp *.h
	$(CC) -c -o main.o main.cpp -lm -ldl -Wall -lfftw3 -O2 

main : $(OBJ) main.o
	$(CC) -o main $(OBJ) main.o -lm -ldl -Wall -lfftw3 -O2

clean :
	rm -f $(OBJ) main.o
	rm -f main

.PHONY : clean
