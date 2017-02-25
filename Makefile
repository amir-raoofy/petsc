MPI=mpic++

CFLAGS = $(PETSC_INC) ${PETSC_CC_INCLUDES}
LFLAGS = $(PETSC_LIB)

all: compile

compile:
	$(MPI) $(CFLAGS) threePointStencil.cpp -o threePointStencil $(PETSC_LIB)
	$(MPI) $(CFLAGS) fivePointStencil.cpp -o fivePointStencil $(PETSC_LIB)

clean:
	rm -rf *.o
	rm -rf threePointStencil fivePointStencil

run:
	mpiexec -n 2 ./fivePointStencil


