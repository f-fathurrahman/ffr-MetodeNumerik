! Start of Main Program
! This is a program to solve the unsteady 1D wave equation 
! using the explicit (forward Euler) method
! Exercise 5.5

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
    
! Open output file
    OPEN(unit=10, file="Exercise5_5.out", status = "unknown")

    dx = length/(n_nodes-1)
    dx2 = dx*dx
    DO i = 1,n_nodes
      x(i) = (i-1)*dx
    ENDDO

! Maximum allowable time step
    dt_max = dx
    dt = dt_max*r  ! actual time step
    dt2 = dt*dt


! Numerical solution: Explicit (forward Euler) Method

! Set initial conditions
    pi = 4.0d0*ATAN(one)
    DO i = 1,n_nodes
      phi(1,i) = SIN(pi*x(i))
      phi_deriv(1,i) = 0.25d0*SIN(2.0d0*pi*x(i))
    ENDDO

    dtx = dt2/dx2

! First time step
    phi(2,1) = zero     
    DO i = 2,n_nodes-1
      phi(2,i) = phi(1,i) + phi_deriv(1,i)*dt + &
                 0.5d0*dtx*(phi(1,i+1)-2.0d0*phi(1,i)+phi(1,i-1))
    ENDDO
    phi(2,n_nodes) = zero     

! Subsequent time steps
    DO n = 3, n_times
      phi(n,1) = zero
      DO i = 2, n_nodes-1   ! Interior nodes
        phi(n,i) = 2.0d0*phi(n-1,i) - phi(n-2,i) + &
                   dtx*(phi(n-1,i+1)-two*phi(n-1,i)+ &
                                     phi(n-1,i-1))
      ENDDO
      phi(n,n_nodes) = zero
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

