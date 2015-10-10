CXX=/usr/local/bin/g++-4.9
#CXX=g++
CXXFLAGS=-fopenmp -fPIC -O3 -DNDEBUG -march=corei7 -msse4.2
INCLUDEFLAGS=-I/usr/include -I/usr/local/include -Imetis-5.1.0/include
LIBFLAGS=-L/usr/local/lib -Lmetis-5.1.0/lib -lm

all: hugecggm-run cggmfast.o sparse.o smat.o util.o

hugecggm-run: cggmfast-run.o cggmfast.o smat.o sparse.o util.o
	$(CXX) ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} $^ -lmetis -o hugecggm_run

cggmfast-run.o: cggmfast-run.cpp cggmfast.h
	$(CXX) -c ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} cggmfast-run.cpp -o cggmfast-run.o

cggmfast.o: cggmfast.h cggmfast.cpp sparse.h smat.h util.h
	${CXX} -c ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} cggmfast.cpp -o cggmfast.o

sparse.o: sparse.h sparse.cpp
	${CXX} -c ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} sparse.cpp -o sparse.o

smat.o: smat.h smat.cpp
	${CXX} -c ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} smat.cpp -o smat.o

util.o: util.h util.cpp
	${CXX} -c ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} util.cpp -o util.o

clean:
	rm -rf *.o hugecggm_warm hugecggm_run *.mexa64 *.mexmaci64

