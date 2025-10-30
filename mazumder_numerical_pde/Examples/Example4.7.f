! Start of Main Program
! PCG Method with point Jacobi pre-conditioning.
! Example 4.7

    program main

    USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 161, m = 161, neqn = n*m
    INTEGER(4), PARAMETER :: max_iter = 1000000
    INTEGER(4) :: i,j,k,iter
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), PARAMETER :: aa = 0.5d0, bb = 0.5d0, cc = 10.0d0, &
                               dd = 1.0d0, ee = 2.0d0
    REAL(8), PARAMETER :: length = 1.0D0
    REAL(8) :: tol = 1.0d-6
    REAL(8), DIMENSION(n) :: x,phi_bot,phi_top
    REAL(8), DIMENSION(m) :: y,phi_left,phi_right
    REAL(8), DIMENSION(:), ALLOCATABLE :: BBB,DDD,EEE,FFF,HHH,QQQ,SSS,YY,DEL,RR
    REAL(8), DIMENSION(:), ALLOCATABLE :: q,d
    REAL(8), DIMENSION(n,m) :: S,phi,phian
    REAL(8) :: dx,dx2,dy,dy2,res1
    REAL(8) :: ae,aw,as,an,ao,source,res,fde,delta,dtq,alpha,beta,delta0
    REAL(4) :: ta(2),time_start,time_end
    

! Open output file
    OPEN(unit=10, file="Example4_7.out", status = "unknown")
    OPEN(unit=20, file="Example4_7.rsl", status = "unknown")

    dx = length/(n-1); dy = length/(m-1)
    dx2 = dx*dx; dy2 = dy*dy
    
    DO i = 1,n
      x(i) = (i-1)*dx
      phi_bot(i) = bb*bb*sinh(-cc*bb)+(x(i)-aa)*(x(i)-aa)*sinh(cc*(x(i)-aa))+dd
      phi_top(i) = (one-bb)*(one-bb)*sinh(cc*(one-bb))+(x(i)-aa)* &
                   (x(i)-aa)*sinh(cc*(x(i)-aa))+dd*exp(ee*x(i))
    ENDDO
    DO j = 1,m
      y(j) = (j-1)*dy
      phi_left(j) = aa*aa*sinh(-cc*aa)+(y(j)-bb)*(y(j)-bb)*sinh(cc*(y(j)-bb))+dd
      phi_right(j) = (one-aa)*(one-aa)*sinh(cc*(one-aa))+(y(j)-bb)* &
                     (y(j)-bb)*sinh(cc*(y(j)-bb))+dd*exp(ee*y(j))
    ENDDO
    DO i = 1,n
      DO j = 1,m
        S(i,j) = 2.0d0*sinh(cc*(x(i)-aa))+4.0d0*(x(i)-aa)*cc*cosh(cc*(x(i)-aa))+ &
                 (x(i)-aa)*(x(i)-aa)*(cc*cc)*sinh(cc*(x(i)-aa)) + &
                 2.0d0*sinh(cc*(y(j)-bb))+4.0d0*(y(j)-bb)*cc*cosh(cc*(y(j)-bb))+ &
                 (y(j)-bb)*(y(j)-bb)*(cc*cc)*sinh(cc*(y(j)-bb)) + &
                 dd*ee*ee*(y(j)*y(j)+x(i)*x(i))*exp(ee*x(i)*y(j))
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

! Point Jacobi Pre-conditioning
    QQQ(:) = QQQ(:)/EEE(:)
    BBB(:) = BBB(:)/EEE(:)
    HHH(:) = HHH(:)/EEE(:)
    DDD(:) = DDD(:)/EEE(:)
    FFF(:) = FFF(:)/EEE(:)
    EEE(:) = one
    

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
         phian(i,j) = ((x(i)-aa)**2)*sinh(cc*(x(i)-aa)) + &
                      ((y(j)-bb)**2)*sinh(cc*(y(j)-bb))+dd*exp(ee*x(i)*y(j))
         WRITE(10,10) x(i),y(j),phian(i,j),phi(i,j),phian(i,j)-phi(i,j)
      ENDDO
    ENDDO

    print *, "time taken = ", time_end-time_start

 10 format(5(1x,F13.6))

    CLOSE(unit=10)
    CLOSE(unit=20)

    END Program main

