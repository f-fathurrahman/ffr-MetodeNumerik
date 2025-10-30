! Start of Main Program
! This is a program to solve a set of 6 algebraic equations using
! the algebraic multi-grid (AMG) method. Only two levels are used. 
! The smoother is Gauss-Seidel.
! Exercise 4.6

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: nff=6
    INTEGER(4), PARAMETER :: ncc=3
    INTEGER(4), PARAMETER :: max_iter=100
    INTEGER(4) :: iter,i,j,k
    REAL(8), DIMENSION(nff) :: phiff,rff,qff,errff
    REAL(8), DIMENSION(ncc) :: errcc,rcc
    REAL(8) :: af(nff,nff), ac(ncc,ncc)
    REAL(8) :: sumres,res,fde,sumoff
    REAL(8), PARAMETER :: tol = 1.0d-8
    

! Open output file
    OPEN(unit=10, file="Exercise4_6.out", status = "unknown")
    OPEN(unit=20, file="Exercise4_6.rsl", status = "unknown")

! Fill in coefficients and RHS of original (fine) equation
    af(1,1)=6.0d0; af(1,2)=-2.0d0; af(1,3)=0.0d0; af(1,4)=-1.0d0; af(1,5)=0.0d0; af(1,6)=0.0d0; qff(1)=3.0d0
    af(2,1)=-2.0d0; af(2,2)=7.0d0; af(2,3)=-1.0d0; af(2,4)=0.0d0; af(2,5)=-1.0d0; af(2,6)=0.0d0; qff(2)=3.0d0
    af(3,1)=0.0d0; af(3,2)=-1.0d0; af(3,3)=7.0d0; af(3,4)=-3.0d0; af(3,5)=0.0d0; af(3,6)=-1.0d0; qff(3)=2.0d0
    af(4,1)=-1.0d0; af(4,2)=0.0d0; af(4,3)=-3.0d0; af(4,4)=7.0d0; af(4,5)=0.0d0; af(4,6)=-2.0d0; qff(4)=1.0d0
    af(5,1)=0.0d0; af(5,2)=0.0d0; af(5,3)=-1.0d0; af(5,4)=-2.0d0; af(5,5)=7.0d0; af(5,6)=-3.0d0; qff(5)=1.0d0
    af(6,1)=0.0d0; af(6,2)=-1.0d0; af(6,3)=-1.0d0; af(6,4)=0.0d0; af(6,5)=-3.0d0; af(6,6)=6.0d0; qff(6)=1.0d0

! Fill in coefficients of coarse equation
    ac(1,1)=9.0d0; ac(1,2)=-2.0d0; ac(1,3)=-1.0d0
    ac(2,1)=-2.0d0; ac(2,2)=8.0d0; ac(2,3)=-3.0d0
    ac(3,1)=-1.0d0; ac(3,2)=-4.0d0; ac(3,3)=7.0d0

    iter = 0

    phiff(:) = 0.0d0

 11 iter = iter + 1   ! Start iteration

! Gauss-Seidel Smoothing for Fine Mesh
    DO i = 1,nff
      sumoff = 0.0d0
      DO j = 1,nff
        IF (i == j) CYCLE
        sumoff = sumoff + af(i,j)*phiff(j)
      ENDDO
      phiff(i) = (qff(i) - sumoff)/af(i,i)
    ENDDO
    
! Calculate Residual for fine mesh. Default of 1 sweep
    DO i = 1,nff
      rff(i) = qff(i)
      DO j = 1,nff
        rff(i) = rff(i) - af(i,j)*phiff(j)
      ENDDO
    ENDDO

! Transfer Residual to coarse mesh (restriction)
    rcc(1) = rff(1)+rff(2)
    rcc(2) = rff(3)+rff(4)
    rcc(3) = rff(5)+rff(6)

! Calculate correction on coarse mesh
    errcc(:) = 0.0d0
    errcc(1) = (rcc(1) - ac(1,2)*errcc(2) - ac(1,3)*errcc(3))/ac(1,1)
    errcc(2) = (rcc(2) - ac(2,1)*errcc(1) - ac(2,3)*errcc(3))/ac(2,2)
    errcc(3) = (rcc(3) - ac(3,1)*errcc(1) - ac(3,2)*errcc(2))/ac(3,3)

! Transfer correction to fine mesh (Prolongation)
    errff(1) = errcc(1); errff(2) = errcc(1)
    errff(3) = errcc(2); errff(4) = errcc(2)
    errff(5) = errcc(3); errff(6) = errcc(3)

! Correct fine grid solution
    phiff(:) = phiff(:) + errff(:)

!   Calculate L2NORM on fine mesh
    sumres = 0.0d0
    DO i = 1,nff
      sumres = sumres + rff(i)*rff(i)
    ENDDO
    res = sqrt(sumres)

    print *, iter, res
    WRITE(20,*) iter, res

    IF(res > tol .AND. iter < max_iter) GO TO 11    ! Convergence check
    
! Print results for fine mesh
    WRITE(10,10) (phiff(i), i = 1,nff)

 10 format(6(1x,F13.6))

    CLOSE(unit=10)
    CLOSE(unit=20)

    END Program main

