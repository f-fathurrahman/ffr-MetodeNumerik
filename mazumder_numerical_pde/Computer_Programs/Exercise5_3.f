! Start of Main Program
! This is a program to solve the unsteady 1D diffusion equation 
! using the Crank-Nicolson method
! Exercise 5.3

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_nodes = 26, n_times = 2000
    INTEGER(4) :: i,k,n
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,zero= 0.0d0, &
                             & length=1.0D0
    REAL(8), DIMENSION(n_nodes) :: x
    REAL(8), DIMENSION(n_nodes) :: a,b,c,d,sol
    REAL(8), DIMENSION(n_times,n_nodes) :: phi,phi_an
    REAL(8), DIMENSION(4) :: l,CC
    REAL(8) :: dx,dx2,err_i,sum_k,term,tt,dt,dt_max,dtx
    
! Open output file
    OPEN(unit=10, file="Exercise5_3.out", status = "unknown")

    dx = length/(n_nodes-1)
    dx2 = dx*dx
    DO i = 1,n_nodes
      x(i) = (i-1)*dx
    ENDDO

! Maximum allowable time step
    dt_max = dx2/two
    dt = dt_max   ! actual time step

! Analaytical solution
    l(1) = 0.8603d0; l(2) = 3.4256d0; l(3) = 6.4373d0; l(4) = 9.5293d0
    DO k = 1,4
      CC(k) = 4.0d0*sin(l(k))/(2.0d0*l(k)+sin(two*l(k)))
    ENDDO
    DO n = 2, n_times
      tt = dt*(n-1)
      DO i = 1, n_nodes
        sum_k = zero
        DO k = 1,4
          term = CC(k)*exp(-l(k)*l(k)*tt)*cos(l(k)*x(i))
          sum_k = sum_k + term
        ENDDO
        phi_an(n,i) = sum_k
      ENDDO
    ENDDO
      

! Numerical solution: Implicit (backward Euler) Method

    phi(1,:) = one  ! Initial Condition
    dtx = dt/dx2

    DO n = 2, n_times
      tt = dt*(n-1)

      d(:) = 0.0d0; a(:) = 0.0d0; c(:) = 0.0d0; b(:) = 0.0d0
! Set up tridiagonal solver
      i = 1 ! Left node
      d(i) = one + dtx
      c(i) = -dtx
      b(i) = (one-dtx)*phi(n-1,i)+dtx*phi(n-1,i+1)
      DO i = 2, n_nodes-1   ! Interior nodes
        d(i) = one + dtx
        a(i-1) = -dtx/two
        c(i) = -dtx/two
        b(i) = (one-dtx)*phi(n-1,i)+dtx*phi(n-1,i-1)/two + &
                                    dtx*phi(n-1,i+1)/two
      ENDDO
      i = n_nodes ! Right node
      d(i) = one + dtx*(one+dx)
      c(i) = zero
      a(i-1) = -dtx
      b(i) = (one-(one+dx)*dtx)*phi(n-1,i)+dtx*phi(n-1,i-1)
           
      CALL TRI(n_nodes,a,d,c,b,sol)
      
      phi(n,:) = sol(:)

    ENDDO

! Print results
    n = INT(0.1d0/dt)
   WRITE(10,*) n, n*dt
    DO i = 1,n_nodes
      err_i = phi_an(n,i)-phi(n,i)
      WRITE(*,10) x(i),phi(n,i),phi_an(n,i),err_i
      WRITE(10,10) x(i),phi(n,i),phi_an(n,i),err_i
    ENDDO
    n = INT(0.2d0/dt)
    WRITE(10,*) n, n*dt
    DO i = 1,n_nodes
      err_i = phi_an(n,i)-phi(n,i)
      WRITE(*,10) x(i),phi(n,i),phi_an(n,i),err_i
      WRITE(10,10) x(i),phi(n,i),phi_an(n,i),err_i
    ENDDO
    n = INT(0.4d0/dt)
    WRITE(10,*) n, n*dt
    DO i = 1,n_nodes
      err_i = phi_an(n,i)-phi(n,i)
      WRITE(*,10) x(i),phi(n,i),phi_an(n,i),err_i
      WRITE(10,10) x(i),phi(n,i),phi_an(n,i),err_i
    ENDDO
    n = INT(0.8d0/dt)
    WRITE(10,*) n, n*dt
    DO i = 1,n_nodes
      err_i = phi_an(n,i)-phi(n,i)
      WRITE(*,10) x(i),phi(n,i),phi_an(n,i),err_i
      WRITE(10,10) x(i),phi(n,i),phi_an(n,i),err_i
    ENDDO

 10 format(4(1x,F14.6))

    CLOSE(unit=10)

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
