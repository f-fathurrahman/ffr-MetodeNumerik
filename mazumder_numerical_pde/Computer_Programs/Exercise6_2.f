! Start of Main Program
! This is a program to solve the one-dimensional linear 2nd order
! ODE using an unequal mesh and the FVM
! Exercise 6.2

    program main


! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_nodes = 81, nc = n_nodes-1
    INTEGER(4) :: iter,n_eq,i
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,pi = 3.1415927d0, &
                             & zero= 0.0d0,phi_left=0.0D0,phi_right=1.0D0, &
                             & length=1.0D0,st=1.02d0
    REAL(8), DIMENSION(n_nodes) :: x
    REAL(8), DIMENSION(nc) :: a,b,c,d,s,dx,phi,xc
    REAL(8) :: phi_an,dx0,dx2,st_fac,dx12,q_left,q_right
    REAL(8) :: d1,d2,den
    
! Open output file
    OPEN(unit=10, file="Exercise6_2.out", status = "unknown")

    n_eq = nc
    st_fac = (st - st**n_nodes)/(one-st)
    dx0 = length/st_fac
    DO i = 1,nc
      dx(i) = dx0*(st**i)
    ENDDO
! Nodes
    x(1) = zero
    DO i = 2,n_nodes-1
      x(i) = x(i-1)+dx(i-1)
    ENDDO
    x(n_nodes) = length
! Cell Centers
    DO i = 1,nc
      xc(i) = (x(i)+x(i+1))/2.0d0
    ENDDO

! Set up coefficient matrix and source vector
    a(:) = zero; b(:) = zero; c(:) = zero; d(:) = zero
    DO i = 2,nc-1
      b(i) = (2.0d0*xc(i)-one)*dx(i)
      a(i-1) = 2.0d0/(dx(i-1)+dx(i))
      c(i) = 2.0d0/(dx(i)+dx(i+1))
      d(i) = -(2.0d0/(dx(i)+dx(i-1))+2.0d0/(dx(i)+dx(i+1)))
    ENDDO

! Boundary conditions
    d1 = (dx(1)+dx(2)/2.0d0)**2
    d2 = (dx(1)/2.0d0)**2
    den = d1*dx(1)/2.0d0 - d2*(dx(1)+dx(2)/2.0d0)
    d(1) = -2.0d0/(dx(1)+dx(2))-d1/den
    c(1) = 2.0d0/(dx(1)+dx(2))+d2/den
    b(1) = (2.0d0*xc(1)-one)*dx(1) - phi_left*(d1-d2)/den

    d1 = (dx(nc)+dx(nc-1)/2.0d0)**2
    d2 = (dx(nc)/2.0d0)**2
    den = d1*dx(nc)/2.0d0 - d2*(dx(nc)+dx(nc-1)/2.0d0)
    d(nc) = -2.0d0/(dx(nc)+dx(nc-1))-d1/den
    a(nc-1) = 2.0d0/(dx(nc)+dx(nc-1))+d2/den
    b(nc) = (2.0d0*xc(nc)-one)*dx(nc) - phi_right*(d1-d2)/den

! Solve system of equations using the Thomas Algorithm
    CALL TRI(n_eq,a,d,c,b,phi)

! Print results
    DO i = 1,nc
      phi_an = (xc(i)**3)/3.0d0 - (xc(i)**2)/2.0d0 + 7.0d0*xc(i)/6.0d0
      WRITE(10,10) xc(i),phi(i),phi_an,(phi_an-phi(i))
      WRITE(*,10) xc(i),phi(i),phi_an,(phi_an-phi(i))
    ENDDO
    d1 = (dx(1)+dx(2)/2.0d0)**2
    d2 = (dx(1)/2.0d0)**2
    den = d1*dx(1)/2.0d0 - d2*(dx(1)+dx(2)/2.0d0)
    q_left = (d1*phi(1)-d2*phi(2)-(d1-d2)*phi_left)/den
    d1 = (dx(nc)+dx(nc-1)/2.0d0)**2
    d2 = (dx(nc)/2.0d0)**2
    den = d1*dx(nc)/2.0d0 - d2*(dx(nc)+dx(nc-1)/2.0d0)
    q_right = -(d1*phi(nc)-d2*phi(nc-1)-(d1-d2)*phi_right)/den
    print *, q_left,q_right,q_left-q_right

    CLOSE(unit=10)

 10 format(4(1x,F15.8))

    END Program main

!--------------------------------------------------------------------------

! Tridiagonal solver (uses Thomas Algorithm)

    SUBROUTINE TRI(n,a,d,c,b,x)

    INTEGER(4), INTENT(in) :: n
    INTEGER(4) :: i
    REAL(8) :: xmult
    REAL(8), DIMENSION(n) :: a,b,c,d,x

    DO i = 2,n
      xmult = a(i-1)/d(i-1)
      d(i) = d(i) - xmult*c(i-1)
      b(i) = b(i) - xmult*b(i-1)
    ENDDO

    x(n) = b(n)/d(n)
    DO i = n-1,1,-1
      x(i) = (b(i)-c(i)*x(i+1))/d(i)
    ENDDO

    END SUBROUTINE TRI
