! Start of Main Program
! This is a program to solve a set of two nonlinear PDEs
! in a segregated manner using Newton iterations (outer) and
! Gauss-Seuidel iterations (inner).
! Exercise 8.8 (a)

    program main

    USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 81, m = 81
    INTEGER(4), PARAMETER :: max_iter = 1000, max_nt = 100
    INTEGER(4) :: i,j,k,iter,nt
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), PARAMETER :: length = 1.0D0
    REAL(8) :: tol_in = 1.0d-6, tol_out = 1.0D-6
    REAL(8), DIMENSION(n) :: x,phi_bot,phi_top,psi_bot,psi_top
    REAL(8), DIMENSION(m) :: y,phi_left,phi_right,psi_left,psi_right
    REAL(8), DIMENSION(n,m) :: phi,psi,fun1,fun2,dphi,dpsi
    REAL(8) :: dx,dx2,dy,dy2
    REAL(8) :: sumnt1,resnt1,sumnt2,resnt2,adv,dif,res
    REAL(8) :: ae,aw,as,an,ao
    REAL(4) :: ta(2),time_start,time_end
    

! Open output file
    OPEN(unit=10, file="Exercise8_8_segregated.out", status = "unknown")
    OPEN(unit=20, file="Exercise8_8_segregated.rsl", status = "unknown")
    OPEN(unit=30, file="Exercise8_8_segregated_newton.rsl", status = "unknown")

    dx = length/(n-1); dy = length/(m-1)
    dx2 = dx*dx; dy2 = dy*dy
    
    DO i = 1,n
      x(i) = (i-1)*dx
      phi_bot(i) = zero
      phi_top(i) = zero
      psi_bot(i) = zero
      psi_top(i) = zero
    ENDDO
    DO j = 1,m
      y(j) = (j-1)*dy
      phi_left(j) = one
      phi_right(j) = zero
      psi_left(j) = one
      psi_right(j) = zero
    ENDDO
    
! Initial Guess
    phi(:,:) = half; psi(:,:) = half

! Boundary conditions
    DO j = 2,m
      phi(1,j) = phi_left(j)
      psi(1,j) = psi_left(j)
    ENDDO
    DO j = 2,m
      phi(n,j)= phi_right(j)
      psi(n,j)= psi_right(j)
    ENDDO
    DO i = 1,n
      phi(i,1) = phi_bot(i)
      psi(i,1) = psi_bot(i)
    ENDDO
    DO i = 2,n-1
      phi(i,m) = phi_top(i)
      psi(i,m) = psi_top(i)
    ENDDO


! Start of Newton's Iteration loop (Outer Loop)

    time_start = etime(ta) 

    DO nt = 1, max_nt

! Compute function for first equation
    DO j = 2,m-1
      DO i = 2,n-1
         adv = (phi(i+1,j)*phi(i+1,j)-phi(i-1,j)*phi(i-1,j))/(two*dx) + &
               (phi(i,j+1)*psi(i,j+1)-phi(i,j-1)*psi(i,j-1))/(two*dy)
         dif = (phi(i+1,j)-two*phi(i,j)+phi(i-1,j))/dx2 + &
               (phi(i,j+1)-two*phi(i,j)+phi(i,j-1))/dy2 
         fun1(i,j) = adv-dif
      ENDDO
    ENDDO

! Compute residual of first non-linear equation
    sumnt1 = zero
    DO j = 2,m-1
      DO i = 2,n-1
         sumnt1 = sumnt1 + fun1(i,j)**2
      ENDDO
    ENDDO

    resnt1 = SQRT(MAX(zero,sumnt1))

! Compute function for second equation
    DO j = 2,m-1
      DO i = 2,n-1
         adv = (phi(i+1,j)*psi(i+1,j)-phi(i-1,j)*psi(i-1,j))/(two*dx) + &
               (psi(i,j+1)*psi(i,j+1)-psi(i,j-1)*psi(i,j-1))/(two*dy)
         dif = (psi(i+1,j)-two*psi(i,j)+psi(i-1,j))/dx2 + &
               (psi(i,j+1)-two*psi(i,j)+psi(i,j-1))/dy2 
         fun2(i,j) = adv-dif
      ENDDO
    ENDDO

! Compute residual of second non-linear equation
    sumnt2 = zero
    DO j = 2,m-1
      DO i = 2,n-1
         sumnt2 = sumnt2 + fun2(i,j)**2
      ENDDO
    ENDDO

    resnt2 = SQRT(MAX(zero,sumnt2))

    WRITE(*,*) "Newton Iteration", nt, resnt1, resnt2
    WRITE(20,*) "Newton Iteration Number", nt
    WRITE(30,*) nt,resnt1,resnt2

    IF(resnt1 < tol_out .and. resnt2 < tol_out) EXIT

! Start of Gauss-Seidel Sweep for first equation
    WRITE(20,*) "First Equation"
    WRITE(*,*) "First Equation"
    dphi(:,:) = zero

    DO iter = 1,max_iter

      DO i = 2, n-1
        DO j = 2,m-1

          aw = -phi(i-1,j)/dx - one/dx2
          ae = +phi(i+1,j)/dx - one/dx2
          as = -psi(i,j-1)/(two*dy) - one/dy2
          an = +psi(i,j+1)/(two*dy) - one/dy2
          ao = (two/dx2+two/dy2)

          dphi(i,j) = (-fun1(i,j) - ae*dphi(i+1,j) - aw*dphi(i-1,j) &
                                  - an*dphi(i,j+1) - as*dphi(i,j-1))/ao
        ENDDO
      ENDDO

! Compute inner iteration residual
      res = zero
      DO i = 2, n-1
        DO j = 2,m-1

          aw = -phi(i-1,j)/dx - one/dx2
          ae = +phi(i+1,j)/dx - one/dx2
          as = -psi(i,j-1)/(two*dy) - one/dy2
          an = +psi(i,j+1)/(two*dy) - one/dy2
          ao = (two/dx2+two/dy2)

          res = res + (-fun1(i,j) - ae*dphi(i+1,j) - aw*dphi(i-1,j) &
                                  - an*dphi(i,j+1) - as*dphi(i,j-1) &
                                  - ao*dphi(i,j))**2
        ENDDO
      ENDDO
      res = SQRT(MAX(zero,res))
      WRITE(20,*) iter,res
      WRITE(*,*) iter,res

      IF(res < tol_in) EXIT

    ENDDO ! GS Iteration for first equation

    phi(:,:) = phi(:,:) + dphi(:,:)  ! Update solution for first equation

! Start of Gauss-Seidel Sweep for second equation
    WRITE(20,*) "Second Equation"
    WRITE(*,*) "Second Equation"
    dpsi(:,:) = zero

    DO iter = 1,max_iter

      DO i = 2, n-1
        DO j = 2,m-1

          as = -psi(i,j-1)/dy - one/dy2
          an = +psi(i,j+1)/dy - one/dy2
          aw = -phi(i-1,j)/(two*dx) - one/dx2
          ae = +phi(i+1,j)/(two*dx) - one/dx2
          ao = (two/dx2+two/dy2)

          dpsi(i,j) = (-fun2(i,j) - ae*dpsi(i+1,j) - aw*dpsi(i-1,j) &
                                  - an*dpsi(i,j+1) - as*dpsi(i,j-1))/ao
        ENDDO
      ENDDO

! Compute inner iteration residual
      res = zero
      DO i = 2, n-1
        DO j = 2,m-1

          as = -psi(i,j-1)/dy - one/dy2
          an = +psi(i,j+1)/dy - one/dy2
          aw = -phi(i-1,j)/(two*dx) - one/dx2
          ae = +phi(i+1,j)/(two*dx) - one/dx2
          ao = (two/dx2+two/dy2)

          res = res + (-fun2(i,j) - ae*dpsi(i+1,j) - aw*dpsi(i-1,j) &
                                  - an*dpsi(i,j+1) - as*dpsi(i,j-1) &
                                  - ao*dpsi(i,j))**2
        ENDDO
      ENDDO
      res = SQRT(MAX(zero,res))
      WRITE(20,*) iter,res
      WRITE(*,*) iter,res

      IF(res < tol_in) EXIT

    ENDDO ! GS Iteration for first equation

    psi(:,:) = psi(:,:) + dpsi(:,:)  ! Update solution for first equation

    ENDDO ! Newton Iteration loop
  
    time_end = etime(ta)

! Print results
    WRITE(10,*) 'VARIABLES = "X", "Y", "PHI", "PSI"'
    WRITE(10,*) 'ZONE I=', n, ', J=', m, ', DATAPACKING=POINT'
    DO j = 1,m
      DO i = 1,n
         WRITE(10,10) x(i),y(j),phi(i,j),psi(i,j)
      ENDDO
    ENDDO

    print *, time_end-time_start

 10 format(4(1x,F13.6))

    CLOSE(unit=10)
    CLOSE(unit=20)
    CLOSE(unit=30)

    END Program main

