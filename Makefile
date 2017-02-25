MPI=mpic++

CFLAGS = $(PETSC_INC) ${PETSC_CC_INCLUDES}
LFLAGS = $(PETSC_LIB)

all: compile

compile:
	$(MPI) $(CFLAGS) main.cpp -o test $(PETSC_LIB)

clean:
	rm -rf *.o
	rm -rf test

run:
	mpiexec -n 2 ./test


