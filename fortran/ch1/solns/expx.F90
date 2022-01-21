    program main
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none
    
    PetscErrorCode ierr
    PetscMPIInt rank
    PetscInt i
    PetscReal x, localval, globalsum
    character(len=80) :: outputString

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    if (ierr .ne. 0) then
      print*,'Unable to initialize PETSc'
      stop
    endif
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRA(ierr)

    ! command line option
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-x",x,PETSC_NULL_BOOL, ierr); CHKERRA(ierr)

    localval = 1.0
    do i = 1, rank
        localval = localval*x/i
    end do

    globalsum = 0.
    call MPI_Allreduce(localval,globalsum,1,MPIU_REAL,MPIU_SUM, &
    &   PETSC_COMM_WORLD, ierr); CHKERRA(ierr)

    write(outputString,*) "e^x is about", globalsum, "\n"
    call PetscPrintf(PETSC_COMM_WORLD,outputString,ierr); CHKERRA(ierr)

    call PetscFinalize(ierr)

    stop
    end program main