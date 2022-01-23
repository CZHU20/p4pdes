#include <petsc/finclude/petsc.h>

    program main
    use petsc
    implicit none
    
    PetscErrorCode ierr
    DM da
    Mat A
    Vec b,u,uexact
    KSP ksp
    PetscReal errnorm

    character(len=80) :: outputString

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    external formMatrix, formRHS, formExact

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    if (ierr .ne. 0) then
      print*,'Unable to initialize PETSc'
      stop
    endif

    call DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
    DMDA_STENCIL_STAR,9,9,PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL_INTEGER, &
    PETSC_NULL_INTEGER,da, ierr); CHKERRQ(ierr)

    ! create linear system matrix A
    call DMSetFromOptions(da, ierr); CHKERRQ(ierr)
    call DMSetUp(da,ierr); CHKERRQ(ierr)
    call DMCreateMatrix(da,A,ierr); CHKERRQ(ierr)
    call MatSetFromOptions(A,ierr); CHKERRQ(ierr)

    ! create RHS b, numerical solution u, exact solution uexact
    call DMCreateGlobalVector(da,b,ierr); CHKERRQ(ierr)
    call VecDuplicate(b,u,ierr); CHKERRQ(ierr)
    call VecDuplicate(b,uexact,ierr); CHKERRQ(ierr)

    ! fill vectors and assemble linear system
    call formExact(da,uexact,ierr); CHKERRQ(ierr)
    call formRHS(da,b,ierr); CHKERRQ(ierr)
    call formMatrix(da,a,ierr); CHKERRQ(ierr)

    write(outputString,*) "Test\n"
    call PetscPrintf(PETSC_COMM_SELF,outputString,ierr); CHKERRA(ierr)

    VecDestroy(u,ierr); CHKERRQ(ierr);
    VecDestroy(b,ierr); CHKERRQ(ierr)
    VecDestroy(uexact,ierr); CHKERRQ(ierr)
    MatDestroy(A,ierr); CHKERRQ(ierr)
    ! KSPDestroy(ksp,ierr); CHKERRQ(ierr)
    DMDestroy(da,ierr); CHKERRQ(ierr)
    call PetscFinalize(ierr); CHKERRQ(ierr)

    stop
    end program main

!   --------------------------------------------------------------------
    subroutine formExact(da,uexact,ierr)
    use petsc
    implicit none

    DM, intent(in) :: da
    Vec, intent(out) :: uexact
    PetscErrorCode, intent(out) :: ierr

    PetscInt i,j, mx, my, xs, xm, ys, ym
    PetscReal hx,hy,x,y
    PetscScalar,pointer :: auexact(:,:)=>null()

    call DMDAGetInfo(da,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER, &
           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)
    CHKERRQ(ierr)

    call DMDAGetCorners(da,xs,ys,PETSC_NULL_INTEGER,xm,ym,PETSC_NULL_INTEGER,ierr)
    CHKERRQ(ierr)

    hx = 1.0_8/REAL((mx-1),kind=8); hy = 1.0_8/REAL((my-1),kind=8)
    
    call DMDAVecGetArrayF90(da,uexact,auexact,ierr); CHKERRQ(ierr)
    do j = ys, ys+ym-1
      y = REAL(j,kind=8)*hy
      do i = xs, xs+xm-1
        x = REAL(i,kind=8)*hx
        auexact(i,j) = x*x * (1.0_8 - x*x) * y*y * (y*y - 1.0_8)
      end do
    end do

    call DMDAVecRestoreArray(da, uexact, auexact,ierr);CHKERRQ(ierr)

    return
    end subroutine formExact
!   --------------------------------------------------------------------
    subroutine formRHS(da,b,ierr)
    use petsc
    implicit none

    DM, intent(in) :: da
    Vec, intent(out) :: b
    PetscErrorCode, intent(out) :: ierr

    PetscInt       i, j, mx, my, xs, xm, ys, ym
    PetscReal      hx, hy, x, y, f
    PetscScalar,pointer :: ab(:,:)=>null()

    call DMDAGetInfo(da,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER, &
           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)
    CHKERRQ(ierr)

    call DMDAGetCorners(da,xs,ys,PETSC_NULL_INTEGER,xm,ym,PETSC_NULL_INTEGER,ierr)
    CHKERRQ(ierr)

    hx = 1.0_8/REAL((mx-1),kind=8); hy = 1.0_8/REAL((my-1),kind=8)
    
    call DMDAVecGetArrayF90(da,b,ab,ierr); CHKERRQ(ierr)
    do j = ys, ys+ym-1
      y = REAL(j,kind=8)*hy
      do i = xs, xs+xm-1
        x = REAL(i,kind=8)*hx
        if (i==0 .or. i==max-1 .or. j==0 .or. j==my-1) then
          ab(i,j) = 0.0
        else
          f = 2.0 * ( (1.0 - 6.0*x*x) * y*y * (1.0 - y*y) &
                  + (1.0 - 6.0*y*y) * x*x * (1.0 - x*x) )
          ab(i,j) = hx * hy * f
      end do
    end do

    call DMDAVecRestoreArray(da, b, ab,ierr);CHKERRQ(ierr)

    return
    end subroutine formRHS
!   --------------------------------------------------------------------
    subroutine formMatrix(da, A, ierr)
    use petsc
    implicit none

    DM, intent(in) :: da
    Mat, intent(out) :: A
    PetscErrorCode, intent(out) :: ierr

    MatStencil     :: row, col(0:4)
    PetscReal     :: hx, hy, v(0:4)
    PetscInt      :: i, j, ncols, mx, my, xs, xm, ys, ym

    call DMDAGetInfo(da,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER, &
           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)
    CHKERRQ(ierr)

    call DMDAGetCorners(da,xs,ys,PETSC_NULL_INTEGER,xm,ym,PETSC_NULL_INTEGER,ierr)
    CHKERRQ(ierr)

    hx = 1.0_8/REAL((mx-1),kind=8); hy = 1.0_8/REAL((my-1),kind=8)

    do j = ys, ys+ym-1
      do i = xs, xs+xm-1
        
      end do
    end do
    
    return
    end subroutine formMatrix