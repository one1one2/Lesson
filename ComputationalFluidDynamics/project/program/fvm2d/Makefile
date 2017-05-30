source = $(wildcard *.cpp)
include = $(wildcard *.h)
object = $(patsubst %.cpp,%.o,$(source))
exe = $(patsubst %.cpp,%,$(source))
LDFLAGS = -std=c++11 -pthread -O3 -ffast-math -march=native -fopenmp 
CXX = g++

all : main


main : $(source) $(include)
	$(CXX) -o $@  $(source) $(LDFLAGS)

clean:
	-rm -f $(exe)
	-rm -f *.o
	-rm -f *.dat

.PHONY: clean
