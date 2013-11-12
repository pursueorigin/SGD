SOFTWARE_HOME=${HOME}/Software

MPICXX  = mpic++
MPIRUN  = mpirun
FORTRAN = gfortran

ifeq ($(shell uname -s),Darwin)
  EXTRAFLAG=-framework vecLib 
endif

DBGFLAG=-g -O0
OPTIFLAG=-O3

INCLIBS = -lm -Lmatrix_market -lmmio $(EXTRALIB)
INCFLAGS= -I${SOFTWARE_HOME}/include
CXXFLAGS=-Wno-write-strings -Wno-format -Wno-deprecated $(EXTRAFLAG)

HEADERS = utilities.hpp \
          parser.hpp \
          partitioner.hpp \
          modeler.hpp \
          solver.hpp \
          evaluator.hpp \
          omp_parallel_reduce.hpp \
          residual.hpp \
          output.hpp \
          lapack_prototypes.hpp \
          weighted_csr_graph.hpp \
          matrix_market.hpp \
          reader.hpp 

all: lib harness-opt

lib:
	cd matrix_market && make lib

harness-opt: harness.cpp $(HEADERS) lib
	${MPICXX} ${CXXFLAGS} ${OPTIFLAG} ${INCFLAGS} $< -o $@ ${INCLIBS}

harness-dbg: harness.cpp $(HEADERS) lib
	${MPICXX} ${CXXFLAGS} ${DBGFLAG} ${INCFLAGS} $< -o $@ ${INCLIBS}

.PHONY: clean

clean:
	cd matrix_market && make clean
	rm -rf *.o harness-* *.dSYM
