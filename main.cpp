#include <iostream>
#include <mpi.h>
#include <petsc.h>

int main(int argc, char *argv[])
{
	//random generator seed
	//srand(time(NULL));

	//MPI vars
	MPI_Comm comm = MPI_COMM_WORLD;
	PetscInt rank;

	//PETSC vars
	KSP solver;
	Mat A;
  PetscErrorCode ierr;
	Vec x,y;
	PetscInt n =10000001,low,high,i,j;
	PetscScalar len,d,u,l,h;
	//PetscReal nrm;

	//initialize petsc and set ranks
	ierr = PetscInitialize(&argc,&argv,0,"USAGE:\n");CHKERRQ(ierr);
	MPI_Comm_rank(comm, &rank );

	//set the parameters
	len =1;
	h=len/(n-1);
	d= 2/(h*h);
	u=-1/(h*h);
	l=-1/(h*h);

	//setting up A
	ierr = MatCreate(comm,&A);;CHKERRQ(ierr);
	ierr = MatSetType(A,MATMPIAIJ);;CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);;CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
	ierr = MatGetOwnershipRange(A,&low,&high);CHKERRQ(ierr);

	for (i = low; i < high; i++) {
		//diag
		MatSetValues(A,1,&i,1,&i,&d,INSERT_VALUES);
		//lower diag
		j=i-1;
		if (j>=0) {
			MatSetValues(A,1,&i,1,&j,&l,INSERT_VALUES);
		}
		//upper diag
		j=i+1;
		if (j<n) {
			MatSetValues(A,1,&i,1,&j,&u,INSERT_VALUES);
		}
	}
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	//MatView(A,0);

	/*
	 *MatTranspose(A,MAT_REUSE_MATRIX,&A);
	 *MatView(A,0);
	 *MatNorm(A,NORM_INFINITY,&nrm);
	 *if (rank==0) {
	 *  std::cout << nrm << std::endl;
	 *}
	 */
	
	//set up the vectors; initial solution and rhs
	VecCreate(comm,&x);
	VecSetSizes(x,PETSC_DECIDE,n);
	VecSetFromOptions(x);
	VecSetUp(x);
	VecSet(x,0.0);
	//VecView(x,0);
	VecAssemblyBegin(x);
	VecAssemblyEnd(x);

	//rhs
	VecDuplicate(x,&y);
	MatMult(A,x,y);
	VecSet(y,2.0);
	//VecView(y,0);

	//ksp
	KSPCreate(comm,&solver);
	KSPSetType(solver,KSPGMRES);
	KSPSetOperators(solver,A,A);
	KSPSetUp(solver);
	KSPSolve(solver,y,x);
	//VecView(x,0);

	KSPConvergedReason reason;
	KSPGetConvergedReason(solver,&reason);
	
	if (reason<0) {
		printf("Divergence.\n");
	}

	if (rank==0) {
		printf("Done.\n");
	}
	

	MatDestroy(&A);
	VecDestroy(&x);
	VecDestroy(&y);
	KSPDestroy(&solver);
	ierr = PetscFinalize();CHKERRQ(ierr);
	return 0;
}
