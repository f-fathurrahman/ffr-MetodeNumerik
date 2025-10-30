! Start of Main Program
! Example 5.8 using MOL & RK2

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_nodes = 21, max_times = 100
    INTEGER(4) :: itimes,i
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0, &
                             & zero= 0.0d0,phi_right=1.0D0, &
                             & length=1.0D0,deriv_left = 1.0d0
    REAL(8), DIMENSION(max_times,n_nodes) :: phi,phia
    REAL(8), DIMENSION(n_nodes) :: x
    REAL(8) :: dx,dx2,pi,dt,f1o,f1e,f1w,f2
    
! Open output file
    OPEN(unit=10, file="Example5_8.out", status = "unknown")
    OPEN(unit=20, file="Example5_8_point.out", status = "unknown")

    pi = 4.0d0*atan(one)
    dx = pi/(n_nodes-1)
    dx2 = dx*dx
    DO i = 1,n_nodes
      x(i) = (i-1)*dx
    ENDDO
    dt = 0.8d0*dx2/two

! Set up boundary conditions

    DO itimes = 1,max_times
      phi(itimes,1) = exp(-(itimes-1)*dt)
      phi(itimes,n_nodes) = -exp(-(itimes-1)*dt)
    ENDDO

! Set up initial condition
    DO i = 1,n_nodes
      phi(1,i) = cos(x(i))
    ENDDO
    
! RK2 Integration

    DO itimes = 1,max_times-1

      phia(itimes+1,1) = exp(-itimes*dt)
      DO i = 2,n_nodes-1
        f1o = (phi(itimes,i+1)-2.0d0*phi(itimes,i)+phi(itimes,i-1))/dx2
        IF(i == n_nodes-1)THEN
          f1e = exp(-(itimes-1)*dt)
        ELSE
          f1e = (phi(itimes,i+2)-2.0d0*phi(itimes,i+1)+phi(itimes,i))/dx2
        ENDIF
        IF(i == 2)THEN
          f1w = -exp(-(itimes-1)*dt)
        ELSE
          f1w = (phi(itimes,i)-2.0d0*phi(itimes,i-1)+phi(itimes,i-2))/dx2
        ENDIF
        f2 = (phi(itimes,i+1)+dt*f1e-2.0d0*(phi(itimes,i)+dt*f1o)+phi(itimes,i-1)+dt*f1w)/dx2
        phi(itimes+1,i) = phi(itimes,i) + dt*(f1o+f2)/2.0d0
        phia(itimes+1,i) = exp(-itimes*dt)*cos(x(i))
      ENDDO
      phia(itimes+1,n_nodes) = -exp(-itimes*dt)

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
