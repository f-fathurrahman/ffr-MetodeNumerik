! Start of Main Program
! This is a program to solve a second order linear elliptic PDE
! using the finite-difference method.
! The solver is the CG Method with ILU0 pre-conditioning.
! Exercise 4.5

    program main

    USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 161, m = 161, neqn = n*m
    INTEGER(4), PARAMETER :: max_iter = 100000
    INTEGER(4) :: i,j,k,iter
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), PARAMETER :: length = 1.0D0
    REAL(8) :: tol = 1.0d-6
    REAL(8), DIMENSION(n) :: x,phi_bot,phi_top
    REAL(8), DIMENSION(m) :: y,phi_left,phi_right
    REAL(8), DIMENSION(:), ALLOCATABLE :: BBB,DDD,EEE,FFF,HHH,QQQ,SSS,YY,RR,RRH
    REAL(8), DIMENSION(:), ALLOCATABLE :: q,d
    REAL(8), DIMENSION(:), ALLOCATABLE :: g,e,f,b,c
    REAL(8), DIMENSION(n,m) :: S,phi,phia
    REAL(8) :: dx,dx2,dy,dy2,res1
    REAL(8) :: ae,aw,as,an,ao,source,res,fde,delta,dtq,alpha,beta,delta0,deltam,deltam0
    REAL(4) :: ta(2),time_start,time_end
    

! Open output file
    OPEN(unit=10, file="Exercise4_5.out", status = "unknown")
    OPEN(unit=20, file="Exercise4_5.rsl", status = "unknown")

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
    DO j = 2,m
      phi(1,j) = phi_left(j)
    ENDDO
    DO j = 2,m
      phi(n,j)= phi_right(j)
    ENDDO
    DO i = 1,n
      phi(i,1) = phi_bot(i)
    ENDDO
    DO i = 2,n-1
      phi(i,m) = phi_top(i)
    ENDDO


! Link Coefficients
    aw = one/dx2
    ae = one/dx2
    as = one/dy2
    an = one/dy2
    ao = -two/dx2-two/dy2

    time_start = etime(ta)
    
! Define coefficient matrix [A], and Right hand Vector [Q]

    ALLOCATE(BBB(neqn),DDD(neqn),EEE(neqn),FFF(neqn),HHH(neqn),QQQ(neqn),SSS(neqn))
    ALLOCATE(RR(neqn),RRH(neqn),YY(neqn))
    ALLOCATE(q(neqn),d(neqn))
    ALLOCATE(g(neqn),e(neqn),f(neqn),b(neqn),c(neqn))

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

! Compute Modified Residual
! First Compute [L] & [U]
    g(:)=zero;e(:)=zero;c(:)=zero;b(:)=zero;f(:)=zero
    g(1) = EEE(1) 
    c(1) = zero
    b(1) = zero
    e(1) = FFF(1)
    f(1) = HHH(1)
    DO k = 2,neqn
      IF(k > n)THEN
        b(k) = BBB(k)/g(k-n)
      ENDIF
      c(k) = DDD(k)/g(k-1)
      IF(k > n)THEN
        g(k) = EEE(k) - b(k)*f(k-n) - c(k)*e(k-1)
      ELSE
        g(k) = EEE(k) - c(k)*e(k-1)
      ENDIF
      f(k) = HHH(k)
      e(k) = FFF(k)
    ENDDO

    RRH(:) = zero
! Forward substitution: [L]*[Y] = [R]
    YY(1) = RR(1)
    DO k = 2,neqn
      IF(k > n)THEN
        YY(k) = (RR(k) - c(k)*YY(k-1) - b(k)*YY(k-n))
      ELSE
        YY(k) = (RR(k) - c(k)*YY(k-1))
      ENDIF
    ENDDO
! Backward Substitution: [U][Rhat] = [Y]
    RRH(:) = zero
    RRH(neqn) = YY(neqn)/g(neqn)
    DO k = neqn-1,1,-1
      IF(k > neqn-n)THEN
        RRH(k) = (YY(k) - e(k)*RRH(k+1))/g(k)
      ELSE
        RRH(k) = (YY(k) - e(k)*RRH(k+1) - f(k)*RRH(k+n))/g(k)
      ENDIF
    ENDDO

    d(:) = RRH(:)

! Compute L2NORM 
    delta = zero
    deltam = zero
    DO k = 1,neqn
      delta = delta + RR(k)*RR(k)
      deltam = deltam + RR(k)*RRH(k)
    ENDDO
    res1 = sqrt(MAX(zero,delta))

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

      alpha = deltam/dtq

! Update Solution
      SSS(:) = SSS(:) + alpha*d(:)

! Compute New Residual Vector
      RR(:) = zero ! Stays zero at boundaries for Dirichlet BC
! Interior Nodes
      DO i = 2, n-1
        DO j = 2,m-1
          k = (j-1)*n + i
          RR(k) = QQQ(k) - EEE(k)*SSS(k) - FFF(k)*SSS(k+1) - HHH(k)*SSS(k+n) - &
                DDD(k)*SSS(k-1) - BBB(k)*SSS(k-n)
        ENDDO
      ENDDO

!Compute New Modified Residual Vector
    RRH(:) = zero
! Forward substitution: [L]*[Y] = [R]
    YY(1) = RR(1)
    DO k = 2,neqn
      IF(k > n)THEN
        YY(k) = (RR(k) - c(k)*YY(k-1) - b(k)*YY(k-n))
      ELSE
        YY(k) = (RR(k) - c(k)*YY(k-1))
      ENDIF
    ENDDO
! Backward Substitution: [U][Rhat] = [Y]
    RRH(:) = zero
    RRH(neqn) = YY(neqn)/g(neqn)
    DO k = neqn-1,1,-1
      IF(k > neqn-n)THEN
        RRH(k) = (YY(k) - e(k)*RRH(k+1))/g(k)
      ELSE
        RRH(k) = (YY(k) - e(k)*RRH(k+1) - f(k)*RRH(k+n))/g(k)
      ENDIF
    ENDDO

! Compute L2NORM
      deltam0 = deltam
      delta = zero
      deltam = zero
      DO k = 1,neqn
        delta = delta + RR(k)**2
        deltam = deltam + RR(k)*RRH(k)
      ENDDO

      beta = deltam/deltam0
      d(:) = RRH(:) + beta*d(:)

      res = sqrt(MAX(zero,delta))/res1
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

