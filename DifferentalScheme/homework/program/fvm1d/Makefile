source = $(wildcard *.cpp)
object = $(patsubst %.cpp,%.o,$(source))
exe = $(patsubst %.cpp,%,$(source))
LDFLAGS = -std=c++11 -fopenmp -pthread -O3
CXX = g++

all : main

main : main.cpp fvm1d.cpp *.h
	$(CXX) -o $@  $(source) $(LDFLAGS)

clean:
	-rm -f $(exe)
	-rm -f *.o
	-rm -f *.dat

.PHONY: clean
