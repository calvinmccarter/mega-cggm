#CXX=/usr/bin/g++
CXX=/usr/local/bin/g++-4.9
CXXFLAGS=-fopenmp -fPIC -O3 -DNDEBUG
INCLUDEFLAGS=-I/usr/include -Ieigen3 -I/usr/local/include
LIBFLAGS= 

all: cggmfast_run cggmfast.o

cggmfast_run: cggmfast_run.cc cggmfast.o
	$(CXX) ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} $^ -o cggmfast_run

cggmfast_run.o: cggmfast.cc cggmfast.h
	$(CXX) -c ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} cggmfast_run.cc -o cggmfast_run.o

cggmfast.o: cggmfast.h cggmfast.cc
	${CXX} -c ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} cggmfast.cc -o cggmfast.o

clean:
	rm -rf *.o cggmfast_run *.mexa64 *.mexmaci64

