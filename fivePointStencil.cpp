#include <petscksp.h>

class PetscSolver
{
public:
	PetscSolver(PetscInt M, PetscInt N);
	~PetscSolver();
	void updateMat();
	void updateRHS();
	void solve();	
private:
  Vec            x,b,u;  /* approx solution, RHS */
  Mat            A;      /* linear system matrix */
  KSP            ksp;    /* linear solver context */
  PetscReal      norm;   /* norm of solution error */
  PetscInt       i,j,Ii,J,Istart,Iend,m,n,its;
  PetscErrorCode ierr;
  PetscScalar    v;
};

PetscSolver::PetscSolver(PetscInt M, PetscInt N){
	m=M;n=N;
	PetscInitialize(0,NULL,0,NULL);

	//create matrix
	MatCreate(PETSC_COMM_WORLD,&A);
  MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m*n,m*n);
  MatSetFromOptions(A);
	MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);
  MatSeqAIJSetPreallocation(A,5,NULL);;
  MatSeqSBAIJSetPreallocation(A,1,5,NULL);
  MatGetOwnershipRange(A,&Istart,&Iend);

	VecCreate(PETSC_COMM_WORLD,&b);
  VecSetSizes(b,PETSC_DECIDE,m*n);
  VecSetFromOptions(b);
  VecDuplicate(b,&u);
  VecDuplicate(b,&x);

	KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetOperators(ksp,A,A);
  KSPSetTolerances(ksp,1.e-2/((m+1)*(n+1)),1.e-50,PETSC_DEFAULT,PETSC_DEFAULT);
  KSPSetFromOptions(ksp);
}

PetscSolver::~PetscSolver(){
	KSPDestroy(&ksp);
  VecDestroy(&u);
	VecDestroy(&x);
  VecDestroy(&b);
	MatDestroy(&A);
	PetscFinalize();
}

void PetscSolver::updateMat(){
 for (Ii=Istart; Ii<Iend; Ii++) {
    v = -1.0; i = Ii/n; j = Ii - i*n;
    if (i>0)   {J = Ii - n; MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);}
    if (i<m-1) {J = Ii + n; MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);}
    if (j>0)   {J = Ii - 1; MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);}
    if (j<n-1) {J = Ii + 1; MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);}
    v = 4.0; MatSetValues(A,1,&Ii,1,&Ii,&v,ADD_VALUES);
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);
}

void PetscSolver::updateRHS(){
	VecSet(u,1.0);
  MatMult(A,u,b);
}

void PetscSolver::solve(){
  KSPSolve(ksp,b,x);
	VecAXPY(x,-1.0,u);
  VecNorm(x,NORM_2,&norm);
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g iterations %D\n",(double)norm,its);
}

int main(){

  PetscSolver petscSolver(8,7);
	petscSolver.updateMat();
	petscSolver.updateRHS();
	petscSolver.solve();

  return 0;
}
