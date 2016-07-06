#CXX=/usr/bin/g++
CXX=/usr/local/bin/g++-4.9
CXXFLAGS=-fopenmp -fPIC -O3 -DNDEBUG
INCLUDEFLAGS=-I/usr/include -I/usr/local/include/eigen3 -I/usr/local/include
LIBFLAGS= 

all: pseudo_run pseudo.o

pseudo_run: pseudo_run.cc pseudo.o
	$(CXX) ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} $^ -o pseudo_run

pseudo_run.o: pseudo.cc pseudo.h
	$(CXX) -c ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} pseudo_run.cc -o pseudo_run.o

pseudo.o: pseudo.h pseudo.cc
	${CXX} -c ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} pseudo.cc -o pseudo.o

clean:
	rm -rf *.o pseudo_run *.mexa64 *.mexmaci64

