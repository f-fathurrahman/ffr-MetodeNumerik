! Start of Main Program
! This is a program to solve the 1D advection-diffusion equation using
! the CD scheme and the Gauss-Seidel method
! Example 4.3 

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_nodes = 11, max_iter = 5
    INTEGER(4) :: i,iter
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,pi = 3.1415927d0, &
                             & zero= 0.0d0,half=0.5d0,alpha=5.0d0
    REAL(8) :: ae,aw,ao,sumres,res
    REAL(8) :: phi(n_nodes) 
    REAL(8), PARAMETER :: p = -0.1d0, tol = 1.0D-6
    
! Open output file
    OPEN(unit=10, file="Example4_3.rsl", status = "unknown")

    phi(:) = zero
    phi(1) = one
    phi(n_nodes) = zero

    DO iter = 1,max_iter

! GS Update
      DO i = 2,n_nodes-1
        ae = (half-p)
        aw = -(half+p)
        ao = two*p
        phi(i) = (-ae*phi(i+1)-aw*phi(i-1))/ao
      ENDDO

! Residual
      sumres = zero
      DO i = 2,n_nodes-1
        ae = (half-p)
        aw = -(half+p)
        ao = two*p
        sumres = sumres +(ao*phi(i)+ae*phi(i+1)+aw*phi(i-1))**2
      ENDDO
      res = SQRT(MAX(sumres,zero))
      print *, iter, res
      WRITE(10,*) iter, res

      IF(res < tol) EXIT

    ENDDO

    CLOSE(unit=10)

    END Program main

