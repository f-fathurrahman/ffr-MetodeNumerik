! Start of Main Program
! Example 5.2: Implicit Method

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_nodes = 21, max_times = 10
    INTEGER(4) :: itimes,i
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0, &
                             & zero= 0.0d0,phi_right=1.0D0, &
                             & length=1.0D0,deriv_left = 1.0d0
    REAL(8), DIMENSION(max_times,n_nodes) :: phi,phia
    REAL(8), DIMENSION(n_nodes) :: x,a,b,c,d,sol
    REAL(8) :: dx,dx2,pi,dt,phi_l,phi_r
    
! Open output file
    OPEN(unit=10, file="Example5_2.out", status = "unknown")
    OPEN(unit=20, file="Example5_2_point.out", status = "unknown")

    pi = 4.0d0*atan(one)
    dx = pi/(n_nodes-1)
    dx2 = dx*dx
    DO i = 1,n_nodes
      x(i) = (i-1)*dx
    ENDDO
    dt = 8.0d0*dx2/two

! Set up boundary conditions

    DO itimes = 1,max_times
      phi(itimes,1) = exp(-(itimes-1)*dt)
      phi(itimes,n_nodes) = -exp(-(itimes-1)*dt)
    ENDDO

! Set up initial condition
    DO i = 1,n_nodes
      phi(1,i) = cos(x(i))
    ENDDO
    
! Implicit Time Advancement

    DO itimes = 1,max_times-1

      phi_l = exp(-(itimes-1)*dt)
      phi_r = -exp(-(itimes-1)*dt)

! Analytical Solution
      phia(itimes+1,1) = exp(-itimes*dt)
      DO i = 2,n_nodes-1
        !phi(itimes+1,i) = phi(itimes,i) + dt*(phi(itimes,i+1)-2.0d0*phi(itimes,i)+phi(itimes,i-1))/dx2
        phia(itimes+1,i) = exp(-itimes*dt)*cos(x(i))
      ENDDO
      phia(itimes+1,n_nodes) = -exp(-itimes*dt)

! Numerical Solution
      d(:) = zero; b(:) = zero; a(:) = zero; c(:) = zero
      d(1) = one
      b(1) = phi_l
      DO i = 2,n_nodes-1
        d(i) = one/dt + two/dx2
        c(i) = -one/dx2
        a(i-1) = -one/dx2
        b(i) = phi(itimes,i)/dt
      ENDDO
      d(n_nodes) = one
      b(n_nodes) = phi_r

      CALL TRI(n_nodes,a,d,c,b,sol)

      phi(itimes+1,:) = sol(:)
        
      WRITE(*,*) "Solution after time step", itimes, itimes*dt
      WRITE(10,*) "Solution after time step", itimes, itimes*dt
      DO i = 1,n_nodes
        WRITE(10,10) itimes, x(i), phi(itimes+1,i),phia(itimes+1,i),phia(itimes+1,i)-phi(itimes+1,i)
        WRITE(*,10) itimes, x(i), phi(itimes+1,i),phia(itimes+1,i),phia(itimes+1,i)-phi(itimes+1,i)
      ENDDO

      WRITE(20,*) itimes*dt, phia(itimes+1,5)-phi(itimes+1,5)

    ENDDO


 10 format(I5,1x,4(F13.6))

    CLOSE(unit=10)
    CLOSE(unit=20)

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
