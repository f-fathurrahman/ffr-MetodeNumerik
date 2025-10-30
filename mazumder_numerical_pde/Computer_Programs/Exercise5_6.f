! Start of Main Program
! This is a program to solve the unsteady 1D wave equation 
! using the implicit (backward Euler) method
! Exercise 5.6

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_nodes = 26, n_times = 16
    INTEGER(4) :: i,k,n
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,zero= 0.0d0, &
                             & length=1.0D0,r=2.0d0
    REAL(8), DIMENSION(n_nodes) :: x
    REAL(8), DIMENSION(n_times,n_nodes) :: phi,phi_deriv,phia
    REAL(8) :: dx,dx2,err_i,sum_k,term,tt,dt,dt_max,dtx,dt2,pi
    REAL(8) , PARAMETER :: phi_l = 0.0d0, phi_r = 0.0d0
    REAL(8), DIMENSION(n_nodes) :: a,b,c,d,sol
    
! Open output file
    OPEN(unit=10, file="Exercise5_6.out", status = "unknown")

    dx = length/(n_nodes-1)
    dx2 = dx*dx
    DO i = 1,n_nodes
      x(i) = (i-1)*dx
    ENDDO

! Maximum allowable time step
    dt_max = dx
    dt = dt_max*r  ! actual time step
    dt2 = dt*dt


! Numerical solution: Implicit (Backward Euler) Method

! Set initial conditions
    pi = 4.0d0*ATAN(one)
    DO i = 1,n_nodes
      phi(1,i) = SIN(pi*x(i))
      phi_deriv(1,i) = 0.25d0*SIN(2.0d0*pi*x(i))
    ENDDO

    dtx = dt2/dx2

! First time step
    d(:) = zero; c(:) = zero; a(:) = zero; b(:) = zero
    sol(:) = zero
    d(1) = one
    b(1) = phi_l
    DO i = 2,n_nodes-1
      d(i) = one+dtx
      c(i) = -0.5d0*dtx
      a(i-1) = -0.5d0*dtx
      b(i) = phi(1,i) + dt*phi_deriv(1,i)
    ENDDO
    d(n_nodes) = one
    b(n_nodes) = phi_r
    CALL TRI(n_nodes,a,d,c,b,sol)
    phi(2,:) = sol(:)

! Subsequent time steps
    DO n = 3, n_times
      d(:) = zero; c(:) = zero; a(:) = zero; b(:) = zero
      sol(:) = zero
      d(1) = one
      b(1) = phi_l
      DO i = 2,n_nodes-1
        d(i) = one+2.0d0*dtx
        c(i) = -dtx
        a(i-1) = -dtx
        b(i) = 2.0d0*phi(n-1,i) - phi(n-2,i)
      ENDDO
      d(n_nodes) = one
      b(n_nodes) = phi_r
      CALL TRI(n_nodes,a,d,c,b,sol)
      phi(n,:) = sol(:)
    ENDDO

! Print results
    DO n = 1,n_times
      tt = dt*(n-1)
      WRITE(*,*) "Time step no. =", n, "     Time =", tt 
      WRITE(10,*) "Time step no. =", n, "     Time =", tt 
      DO i = 1,n_nodes
        phia(n,i) = cos(pi*tt)*sin(pi*x(i)) + &
                    sin(2.0d0*pi*tt)*sin(2.0d0*pi*x(i))/(8.0d0*pi)
        WRITE(*,10) x(i),phia(n,i),phi(n,i),phia(n,i)-phi(n,i)
        WRITE(10,10) x(i),phia(n,i),phi(n,i),phia(n,i)-phi(n,i)
      ENDDO
    ENDDO

 10 format(4(1x,F14.6))

    CLOSE(unit=10)

    END Program main
!---------------------------------------

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
