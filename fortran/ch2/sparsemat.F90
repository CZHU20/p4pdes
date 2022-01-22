    program main
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    PetscErrorCode ierr
    Mat A
    PetscInt i1(3),j1(3),i2,j2(3),i3,j3
    PetscReal aA1(9), aA2(3), aA3

    i1 = (/0,1,2/)
    j1 = i1
    i2 = 3
    j2 = (/1,2,3/)
    i3 = 1
    j3 = 3

    aA1 = (/ 1.0,  2.0,  3.0, &
             2.0,  1.0, -2.0, &
            -1.0,  1.0,  1.0/)
    
    aA2 = (/1.0, 1.0, -1.0/)
    aA3 = -3.0

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    if (ierr .ne. 0) then
        print*,'Unable to initialize PETSc'
        stop
    endif
    
    call MatCreate(PETSC_COMM_WORLD,A, ierr); CHKERRQ(ierr)
    call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,4,4,ierr); CHKERRQ(ierr)
    call MatSetFromOptions(A,ierr); CHKERRQ(ierr)
    call MatSetUp(A,ierr); CHKERRQ(ierr)
    call MatSetValues(A,3,i1,3,j1,aA1,INSERT_VALUES,ierr); CHKERRQ(ierr)
    call MatSetValues(A,1,i2,3,j2,aA2,INSERT_VALUES,ierr); CHKERRQ(ierr)
    call MatSetValue(A,i3,j3,aA3,INSERT_VALUES,ierr); CHKERRQ(ierr) 
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)

    call MatDestroy(A,ierr); CHKERRQ(ierr)
    call PetscFinalize(ierr); CHKERRQ(ierr)

    stop
    end program main