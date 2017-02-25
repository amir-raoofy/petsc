MPI=mpic++

CFLAGS = $(PETSC_INC) ${PETSC_CC_INCLUDES}
LFLAGS = $(PETSC_LIB)

all: compile

compile:
	$(MPI) $(CFLAGS) threePointStencil.cpp -o threePointStencil $(PETSC_LIB)

clean:
	rm -rf *.o
	rm -rf threePointStencil

run:
	mpiexec -n 2 ./threePointStencil


