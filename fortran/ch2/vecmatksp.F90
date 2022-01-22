    program main
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    PetscErrorCode ierr
    Vec x, b
    Mat A
    KSP ksp
    PetscInt i, j(4)
    PetscReal ab(4)
    PetscReal aA(0:3,0:3)

    j = (/0,1,2,3/)
    ab = (/7.0, 1.0, 1.0, 3.0/)
    aA = reshape((/ 1.0,  2.0,  3.0,  0.0, &
                    2.0,  1.0, -2.0, -3.0, &
                   -1.0,  1.0,  1.0,  0.0, &
                    0.0,  1.0,  1.0, -1.0/), shape(aA))
    
    aA = transpose(aA)

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    if (ierr .ne. 0) then
        print*,'Unable to initialize PETSc'
        stop
    end if

    call VecCreate(PETSC_COMM_WORLD,b,ierr); CHKERRA(ierr)
    call VecSetSizes(b,PETSC_DECIDE,4,ierr); CHKERRA(ierr)
    call VecSetFromOptions(b,ierr); CHKERRA(ierr)
    call VecSetValues(b,4,j,ab,INSERT_VALUES,ierr); CHKERRA(ierr)
    call VecAssemblyBegin(b,ierr); CHKERRA(ierr)
    call VecAssemblyEnd(b,ierr); CHKERRA(ierr)

    call MatCreate(PETSC_COMM_WORLD,A,ierr); CHKERRA(ierr)
    call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,4,4,ierr); CHKERRA(ierr)
    call MatSetFromOptions(A,ierr); CHKERRA(ierr)
    call MatSetUp(A,ierr); CHKERRA(ierr)
    do i = 0, 3
        call MatSetValues(A,1,i,4,j,aA(i,:),INSERT_VALUES,ierr); CHKERRA(ierr)
    end do
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr); CHKERRA(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr); CHKERRA(ierr)

    call KSPCreate(PETSC_COMM_WORLD,ksp,ierr); CHKERRA(ierr)
    call KSPSetOperators(ksp,A,A,ierr); CHKERRA(ierr)
    call KSPSetFromOptions(ksp,ierr); CHKERRA(ierr)
    call VecDuplicate(b,x,ierr); CHKERRA(ierr)
    call KSPSolve(ksp,b,x,ierr); CHKERRA(ierr)
    call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRA(ierr)

    call KSPDestroy(ksp,ierr); CHKERRA(ierr)
    call MatDestroy(A,ierr); CHKERRA(ierr)
    call VecDestroy(x,ierr); CHKERRA(ierr)
    call VecDestroy(b,ierr); CHKERRA(ierr)

    call PetscFinalize(ierr); CHKERRA(ierr)
    stop
    end program main