! Start of Main Program
! This program solves a spatially 1D parabolic PDE
! using the Method of Lines in conjuction with RK2
! Exercise 5.7

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_nodes = 26, max_times = 2000
    INTEGER(4) :: itimes,i,k,n
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0, &
                             & zero= 0.0d0,length=1.0D0
    REAL(8), DIMENSION(max_times,n_nodes) :: phi,phi_an
    REAL(8), DIMENSION(n_nodes) :: x,f1,f2
    REAL(8), DIMENSION(4) :: l,CC
    REAL(8) :: dx,dx2,pi,dt,dt_max,dtx,sum_k,term,err_i,tt
    
! Open output file
    OPEN(unit=10, file="Exercise5_7.out", status = "unknown")

    dx = length/(n_nodes-1)
    dx2 = dx*dx
    DO i = 1,n_nodes
      x(i) = (i-1)*dx
    ENDDO

! Maximum allowable time step
    dt_max = dx2/two
    dt = dt_max/two   ! actual time step

! Analaytical solution
    l(1) = 0.8603d0; l(2) = 3.4256d0; l(3) = 6.4373d0; l(4) = 9.5293d0
    DO k = 1,4
      CC(k) = 4.0d0*sin(l(k))/(2.0d0*l(k)+sin(two*l(k)))
    ENDDO
    DO n = 2, max_times
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

    phi(1,:) = one  ! Initial Condition

    dtx = dt/dx2
    
! RK2 Integration

    DO itimes = 1,max_times-1

      f1(1) = two*(phi(itimes,2)-phi(itimes,1))/dx2
      DO i = 2,n_nodes-1
        f1(i) = (phi(itimes,i+1)-2.0d0*phi(itimes,i)+phi(itimes,i-1))/dx2
      ENDDO
      f1(n_nodes) = two*(phi(itimes,n_nodes-1)-(one+dx)*phi(itimes,n_nodes))/dx2

      f2(1) = two*(phi(itimes,2)+dt*f1(2)-phi(itimes,1)-dt*f1(1))/dx2
      DO i = 2,n_nodes-1
        f2(i) = (phi(itimes,i+1)+dt*f1(i+1)-2.0d0*(phi(itimes,i)+dt*f1(i))+ &
                 phi(itimes,i-1)+dt*f1(i-1))/dx2
      ENDDO
      f2(n_nodes) = two*(phi(itimes,n_nodes-1)+dt*f1(n_nodes-1)- &
                         (one+dx)*(phi(itimes,n_nodes)+dt*f1(n_nodes)))/dx2

      DO i = 1,n_nodes
        phi(itimes+1,i) = phi(itimes,i) + dt*(f1(i)+f2(i))/2.0d0
         !phi(itimes+1,i) = phi(itimes,i) + dt*f1(i)
      ENDDO

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
