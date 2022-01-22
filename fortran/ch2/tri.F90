    program main
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    PetscErrorCode ierr
    Vec x,b,xexact
    Mat A
    KSP ksp
    PetscInt m, i, Istart, Iend, j(0:2)
    PetscReal v(0:2), xval, errnorm
    Character(len=80) outString

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    if (ierr .ne. 0) then
        print *, "Unable to initialize petsc"
        stop
    end if

    m = 4
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-tri_m",m,PETSC_NULL_BOOL, ierr); CHKERRQ(ierr)

    call VecCreate(PETSC_COMM_WORLD,x,ierr); CHKERRQ(ierr)
    call VecSetSizes(x,PETSC_DECIDE,m,ierr); CHKERRQ(ierr)
    call VecSetFromOptions(x,ierr); CHKERRQ(ierr)
    call VecDuplicate(x,b,ierr); CHKERRQ(ierr)
    call VecDuplicate(x,xexact,ierr); CHKERRQ(ierr)

    call MatCreate(PETSC_COMM_WORLD,A, ierr); CHKERRQ(ierr)
    call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,m,ierr); CHKERRQ(ierr)
    call MatSetFromOptions(A,ierr); CHKERRQ(ierr)
    call MatSetUp(A,ierr); CHKERRQ(ierr)
    call MatGetOwnershipRange(A,Istart,Iend,ierr); CHKERRQ(ierr)
    Iend = Iend - 1
    do i = Istart, Iend
        if (i == 0 ) then
            v(0) = 3.0
            v(1) = -1.0
            j(0) = 0
            j(1) = 1
            call MatSetValues(A,1,i,2,j,v,INSERT_VALUES,ierr); CHKERRQ(ierr)
        else
            v(0) = -1.0
            v(1) = 3.0
            v(2) = -1.0
            j(0) = i-1
            j(1) = i
            j(2) = i+1
            if (i == m-1) then
                call MatSetValues(A,1,i,2,j,v,INSERT_VALUES,ierr); CHKERRQ(ierr)
            else
                call MatSetValues(A,1,i,3,j,v,INSERT_VALUES,ierr); CHKERRQ(ierr)
            end if
        end if
        xval = EXP(COS((REAL(i,kind=8))))
        call VecSetValue(xexact,i,xval,INSERT_VALUES,ierr); CHKERRQ(ierr)
    end do

    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
    call VecAssemblyBegin(xexact,ierr); CHKERRQ(ierr)
    call VecAssemblyEnd(xexact,ierr); CHKERRQ(ierr)
    call MatMult(A,xexact,b,ierr); CHKERRQ(ierr)

    call KSPCreate(PETSC_COMM_WORLD,ksp,ierr); CHKERRQ(ierr)
    call KSPSetOperators(ksp,A,A,ierr); CHKERRQ(ierr)
    call KSPSetFromOptions(ksp,ierr); CHKERRQ(ierr)
    call KSPSolve(ksp,b,x,ierr); CHKERRQ(ierr)

    call VecAXPY(x,-1.0_8,xexact,ierr); CHKERRQ(ierr)
    call VecNorm(x,NORM_2,errnorm,ierr); CHKERRQ(ierr)
    write(outString,*) "error for m =",m,"system is |x-xexact|_2 =",errnorm,"\n"
    call PetscPrintf(PETSC_COMM_WORLD,outString,ierr); CHKERRQ(ierr)

    call KSPDestroy(ksp,ierr); CHKERRQ(ierr)
    call MATDestroy(A,A,ierr); CHKERRQ(ierr)
    call VecDestroy(x,ierr); CHKERRQ(ierr)
    call VecDestroy(xexact,ierr); CHKERRQ(ierr)

    call PetscFinalize(ierr); CHKERRQ(ierr)
    stop
    end program main