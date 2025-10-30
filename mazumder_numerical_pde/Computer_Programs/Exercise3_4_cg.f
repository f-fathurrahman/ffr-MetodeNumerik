! Start of Main Program
! This is a program to solve a second order linear elliptic PDE
! using the finite-difference method and the Conjugate Gradient method
! Exercise 3.4(f)

    program main

    USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 161, m = 161, neqn = n*m
    INTEGER(4), PARAMETER :: max_iter = 1000000
    INTEGER(4) :: i,j,k,iter
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), PARAMETER :: length = 1.0D0
    REAL(8), PARAMETER  :: tol = 1.0d-6
    REAL(8), DIMENSION(n) :: x,phi_bot,phi_top
    REAL(8), DIMENSION(m) :: y,phi_left,phi_right
    REAL(8), DIMENSION(:), ALLOCATABLE :: BBB,DDD,EEE,FFF,HHH,QQQ,SSS,YY,DEL,RR
    REAL(8), DIMENSION(:), ALLOCATABLE :: q,d
    REAL(8), DIMENSION(n,m) :: S,phi,phia
    REAL(8) :: dx,dx2,dy,dy2
    REAL(8) :: ae,aw,as,an,ao,source,res,fde,delta,dtq,alpha,beta,delta0
    REAL(4) :: ta(2),time_start,time_end
    

! Open output file
    OPEN(unit=10, file="Exercise3_4_cg.out", status = "unknown")
    OPEN(unit=20, file="Exercise3_4_cg.rsl", status = "unknown")

    dx = length/(n-1); dy = length/(m-1)
    dx2 = dx*dx; dy2 = dy*dy
    
    DO i = 1,n
      x(i) = (i-1)*dx
      phi_bot(i) = 100.0d0*x(i)+500.0d0*exp(-50.0*(one-x(i))**2) !bottom
      phi_top(i) = 500.0d0*exp(-50.0*((one-x(i))**2+one)) !top
    ENDDO
    DO j = 1,m
      y(j) = (j-1)*dy
      phi_left(j) = 500.0d0*exp(-50.0*(one+y(j)**2)) !left
      phi_right(j) = 100*(one-y(j))+500.0d0*exp(-50.0*y(j)**2) !right
    ENDDO
    DO i = 1,n
      DO j = 1,m
        S(i,j) = 5.0D4*(100.0d0*((one-x(i))**2+y(j)**2)-two)* &
                 exp(-50.0d0*((one-x(i))**2+y(j)**2))
      ENDDO
    ENDDO

    phi(:,:) = zero

! Boundary conditions
    DO j = 2,m-1
      phi(1,j) = phi_left(j)
    ENDDO
    DO j = 2,m-1
      phi(n,j)= phi_right(j)
    ENDDO
    DO i = 1,n
      phi(i,1) = phi_bot(i)
    ENDDO
    DO i = 1,n
      phi(i,m) = phi_top(i)
    ENDDO


! Link Coefficients
    aw = one/dx2
    ae = one/dx2
    as = one/dy2
    an = one/dy2
    ao = -two/dx2-two/dy2

! Start of iteration loop

    time_start = etime(ta)
    
! Define coefficient matrix [A], and Right hand Vector [Q]

    ALLOCATE(BBB(neqn),DDD(neqn),EEE(neqn),FFF(neqn),HHH(neqn),QQQ(neqn),SSS(neqn))
    ALLOCATE(RR(neqn))
    ALLOCATE(q(neqn),d(neqn))

    EEE(:) = zero; BBB(:) = zero; DDD(:) = zero; FFF(:) = zero; HHH(:) = zero
    QQQ(:) = zero

! Interior Nodes
    DO i = 2, n-1
      DO j = 2,m-1
        k = (j-1)*n + i
        EEE(k) = ao
        BBB(k) = as
        HHH(k) = an
        DDD(k) = aw
        FFF(k) = ae
        QQQ(k) = -S(i,j)  ! Copy to 1D array
        SSS(k) = phi(i,j)  ! Copy to 1D array
      ENDDO
    ENDDO

! Left Boundary
    i = 1
    DO j = 1,m
      k = (j-1)*n + i
      EEE(k) = one
      QQQ(k) = phi_left(j)
      SSS(k) = phi_left(j)
    ENDDO
! Right Boundary
    i = n
    DO j = 1,m
      k = (j-1)*n + i
      EEE(k) = one
      QQQ(k) = phi_right(j)
      SSS(k) = phi_right(j)
    ENDDO
! Bottom Boundary
    j = 1
    DO i = 1,n
      k = (j-1)*n + i
      EEE(k) = one
      QQQ(k) = phi_bot(i)
      SSS(k) = phi_bot(i)
    ENDDO
! Top Boundary
    j = m
    DO i = 1,n
      k = (j-1)*n + i
      EEE(k) = one
      QQQ(k) = phi_top(i)
      SSS(k) = phi_top(i)
    ENDDO

! Solve system using Conjugate Gradient Method

! Preparation Phase

! Compute Residual Vector
    RR(:) = zero ! Stays zero at boundaries for Dirichlet BC
! Interior Nodes
    DO i = 2, n-1
      DO j = 2,m-1
        k = (j-1)*n + i
        RR(k) = QQQ(k) - EEE(k)*SSS(k) - FFF(k)*SSS(k+1) - HHH(k)*SSS(k+n) - &
                DDD(k)*SSS(k-1) - BBB(k)*SSS(k-n)
      ENDDO
    ENDDO
    d(:) = RR(:)

! Compute L2NORM 
    delta = zero
    DO k = 1,neqn
      delta = delta + RR(k)*RR(k)
    ENDDO

! Start of Iteration loop
    
    DO iter = 1,max_iter

! Compute vector q
      q(:) = zero 
      DO i = 2, n-1
        DO j = 2,m-1
          k = (j-1)*n + i
          q(k) =  EEE(k)*d(k) + FFF(k)*d(k+1) + HHH(k)*d(k+n) + &
                  DDD(k)*d(k-1) + BBB(k)*d(k-n)
        ENDDO
      ENDDO

! Compute dT*q
      dtq = zero
      DO k = 1,neqn
        dtq = dtq + d(k)*q(k)
      ENDDO

      alpha = delta/dtq

! Update Solution
      SSS(:) = SSS(:) + alpha*d(:)

! Compute Residual Vector
      RR(:) = zero ! Stays zero at boundaries for Dirichlet BC
! Interior Nodes
      DO i = 2, n-1
        DO j = 2,m-1
          k = (j-1)*n + i
          RR(k) = QQQ(k) - EEE(k)*SSS(k) - FFF(k)*SSS(k+1) - HHH(k)*SSS(k+n) - &
                DDD(k)*SSS(k-1) - BBB(k)*SSS(k-n)
        ENDDO
      ENDDO

! Compute L2NORM
      delta0 = delta
      delta = zero
      DO k = 1,neqn
        delta = delta + RR(k)**2
      ENDDO

      beta = delta/delta0
      d(:) = RR(:) + beta*d(:)

      res = sqrt(MAX(zero,delta))
      print *, iter, res
      WRITE(20,*) iter,res

      IF(res < tol) EXIT

    ENDDO !Iteration loop
  
    time_end = etime(ta)
    
! Print results
    WRITE(10,*) 'VARIABLES = "X", "Y", "PHIAN", "PHI", "Error"'
    WRITE(10,*) 'ZONE I=', n, ', J=', m, ', DATAPACKING=POINT'
    DO j = 1,m
      DO i = 1,n
         k = (j-1)*N+i
         phi(i,j) = SSS(k)
         phia(i,j) = 500.0d0*exp(-50.0*((one-x(i))**2+y(j)**2)) + 100.0d0*x(i)*(one-y(j))
         WRITE(10,10) x(i),y(j),phia(i,j),phi(i,j),phia(i,j)-phi(i,j)
      ENDDO
    ENDDO

    print *, "time taken = ", time_end-time_start

 10 format(5(1x,F13.6))

    CLOSE(unit=10)
    CLOSE(unit=20)

    END Program main

