source = $(wildcard *.cpp)
object = $(patsubst %.cpp,%.o,$(source))
exe = $(patsubst %.cpp,%,$(source))
LDFLAGS = -std=c++11 -fopenmp -pthread -O3 -ffast-math -march=native
CXX = g++

all : main

main : main.cpp fvm1d.cpp *.h
	$(CXX) -o $@  $(source) $(LDFLAGS)
	@echo "Usage: ./main <N> <time>"

clean:
	-rm -f $(exe)
	-rm -f *.o
	-rm -f *.dat

.PHONY: clean
