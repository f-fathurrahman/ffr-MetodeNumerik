! Start of Main Program
! Newton's method with Stone's strongly implicit solver
! Example 8.10

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 101, m = 101, neqn = n*m
    INTEGER(4), PARAMETER :: max_iter = 1000, max_nt = 100
    INTEGER(4) :: i,j,k,iter,nt
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), PARAMETER :: length = 1.0D0
    REAL(8) :: tol_in = 1.0d-6, tol_out = 1.0D-6, alpha = 0.92d0
    REAL(8), DIMENSION(n) :: x,phi_bot,phi_top
    REAL(8), DIMENSION(m) :: y,phi_left,phi_right
    REAL(8), DIMENSION(:), ALLOCATABLE :: BBB,DDD,EEE,FFF,HHH,QQQ,SSS,YY,DEL,RR
    REAL(8), DIMENSION(:), ALLOCATABLE :: d,e,f,b,c
    REAL(8), DIMENSION(n,m) :: S,phi,fun
    REAL(8) :: dx,dx2,dy,dy2
    REAL(8) :: sumnt,resnt,adv,dif,res
    REAL(8) :: ae,aw,as,an,ao
    
! Open output file
    OPEN(unit=10, file="Example8_10.out", status = "unknown")
    OPEN(unit=20, file="Example8_10_stone.rsl", status = "unknown")
    OPEN(unit=30, file="Example8_10_newton.rsl", status = "unknown")

    dx = length/(n-1); dy = length/(m-1)
    dx2 = dx*dx; dy2 = dy*dy
    
    DO i = 1,n
      x(i) = (i-1)*dx
      phi_bot(i) = zero
      phi_top(i) = one
    ENDDO
    DO j = 1,m
      y(j) = (j-1)*dy
      phi_left(j) = zero
      phi_right(j) = one
    ENDDO
    
! Initial Guess
    phi(:,:) = half

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

! Allocate memory
    ALLOCATE(BBB(neqn),DDD(neqn),EEE(neqn),FFF(neqn),HHH(neqn),QQQ(neqn),SSS(neqn))
    ALLOCATE(YY(neqn),DEL(neqn),RR(neqn))
    ALLOCATE(d(neqn),e(neqn),f(neqn),b(neqn),c(neqn))
    
! Start of Newton's Iteration loop (Outer Loop)

    DO nt = 1, max_nt

! Compute function
    DO j = 2,m-1
      DO i = 2,n-1
         adv = (phi(i+1,j)*phi(i+1,j)-phi(i-1,j)*phi(i-1,j))/(two*dx) + &
               (phi(i,j+1)*phi(i,j+1)-phi(i,j-1)*phi(i,j-1))/(two*dy)
         dif = half*(phi(i+1,j)*phi(i+1,j)-two*phi(i,j)*phi(i,j)+phi(i-1,j)*phi(i-1,j))/dx2 + &
               half*(phi(i,j+1)*phi(i,j+1)-two*phi(i,j)*phi(i,j)+phi(i,j-1)*phi(i,j-1))/dy2 + &
               (phi(i+1,j)-two*phi(i,j)+phi(i-1,j))/dx2 + &
               (phi(i,j+1)-two*phi(i,j)+phi(i,j-1))/dy2
         fun(i,j) = adv-dif
      ENDDO
    ENDDO

! Compute residual of non-linear equation
    sumnt = zero
    DO j = 2,m-1
      DO i = 2,n-1
         sumnt = sumnt + fun(i,j)**2
      ENDDO
    ENDDO

    resnt = SQRT(MAX(zero,sumnt))
    WRITE(*,*) "Newton Iteration", nt, resnt
    WRITE(20,*) "Newton Iteration Number", nt
    WRITE(30,*) nt,resnt

    IF(resnt < tol_out) EXIT

! Start of Stone's Implicit Method loop
    
! Define coefficient matrix [A], and Right hand Vector [Q]

    EEE(:) = zero; BBB(:) = zero; DDD(:) = zero; FFF(:) = zero; HHH(:) = zero
    QQQ(:) = zero; RR(:) = zero; YY(:) = zero; DEL(:) = zero


! Interior Nodes
    DO i = 2, n-1
      DO j = 2,m-1
        k = (j-1)*n + i
! Link Coefficients
        aw = -phi(i-1,j)/dx - (one+phi(i-1,j))/dx2
        ae = +phi(i+1,j)/dx - (one+phi(i+1,j))/dx2
        as = -phi(i,j-1)/dy - (one+phi(i,j-1))/dy2
        an = +phi(i,j+1)/dy - (one+phi(i,j+1))/dy2
        ao = (one+phi(i,j))*(two/dx2+two/dy2)

! Set up matrices for Stone's Method
        EEE(k) = ao
        BBB(k) = as
        HHH(k) = an
        DDD(k) = aw
        FFF(k) = ae
        QQQ(k) = -fun(i,j)  ! Copy to 1D array
        SSS(k) = zero  ! Copy to 1D array
      ENDDO
    ENDDO

! Left Boundary
    i = 1
    DO j = 1,m
      k = (j-1)*n + i
      EEE(k) = one
      QQQ(k) = -(phi(i,j)-phi_left(j))
      SSS(k) = zero
    ENDDO
! Right Boundary
    i = n
    DO j = 1,m
      k = (j-1)*n + i
      EEE(k) = one
      QQQ(k) = -(phi(i,j)-phi_right(j))
      SSS(k) = zero
    ENDDO
! Bottom Boundary
    j = 1
    DO i = 1,n
      k = (j-1)*n + i
      EEE(k) = one
      QQQ(k) = -(phi(i,j)-phi_bot(i))
      SSS(k) = zero
    ENDDO
! Top Boundary
    j = m
    DO i = 1,n
      k = (j-1)*n + i
      EEE(k) = one
      QQQ(k) = -(phi(i,j)-phi_top(i))
      SSS(k) = zero
    ENDDO

! Solve system using Stone's strongly implicit method

! Compute coefficients of L and U
    d(:)=zero;e(:)=zero;c(:)=zero;b(:)=zero;f(:)=zero
    d(1) = EEE(1) 
    c(1) = zero
    b(1) = zero
    e(1) = FFF(1)/d(1)
    f(1) = HHH(1)/d(1)
    DO k = 2,neqn
      IF(k > n)THEN
        b(k) = BBB(k)/(one+alpha*e(k-n))
      ENDIF
      c(k) = DDD(k)/(one+alpha*f(k-1))
      IF(k > n)THEN
        d(k) = EEE(k) + alpha*(b(k)*e(k-n)+c(k)*f(k-1))-b(k)*f(k-n)- &
               c(k)*e(k-1)
      ELSE
        d(k) = EEE(k) + alpha*(c(k)*f(k-1))-c(k)*e(k-1)
      ENDIF
      f(k) = (HHH(k)-alpha*c(k)*f(k-1))/d(k)
      IF(k > n)THEN
        e(k) = (FFF(k) - alpha*b(k)*e(k-n))/d(k)
      ELSE
        e(k) = FFF(k)/d(k)
      ENDIF
    ENDDO

! Start of Sone's Method Iteration loop (Inner loop)
    iter = 0
 90 iter = iter + 1

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

! Forward substitution: L*Y = R
    YY(1) = RR(1)/d(1)
    DO k = 2,neqn
      IF(k > n)THEN
        YY(k) = (RR(k) - c(k)*YY(k-1) - b(k)*YY(k-n))/d(k)
      ELSE
        YY(k) = (RR(k) - c(k)*YY(k-1))/d(k)
      ENDIF
    ENDDO

! Backward Substitution: U*Delta = Y
    DEL(:) = zero
    DEL(neqn) = YY(neqn)
    DO k = neqn-1,1,-1
      IF(k > neqn-n)THEN
        DEL(k) = YY(k) - e(k)*DEL(k+1)
      ELSE
        DEL(k) = YY(k) - e(k)*DEL(k+1) - f(k)*DEL(k+n)
      ENDIF
    ENDDO

! Update Solution
    SSS(:) = SSS(:) + DEL(:)

! Compute L2NORM
    res = zero
    DO k = 1,neqn
      res = res + RR(k)**2
    ENDDO
    res = sqrt(MAX(zero,res))
    print *, iter, res
    WRITE(20,*) iter,res

    IF(res > tol_in .and. iter < max_iter) GO TO 90

! Update Solution
    DO j = 2,m-1
      DO i = 2,n-1
         k = (j-1)*N+i
         phi(i,j) = phi(i,j) + SSS(k)
      ENDDO
    ENDDO

    ENDDO ! Newton Iteration loop
  
    
! Print results
    WRITE(10,*) 'VARIABLES = "X", "Y", "PHI"'
    WRITE(10,*) 'ZONE I=', n, ', J=', m, ', DATAPACKING=POINT'
    DO j = 1,m
      DO i = 1,n
         WRITE(10,10) x(i),y(j),phi(i,j)
      ENDDO
    ENDDO


 10 format(3(1x,F13.6))

    CLOSE(unit=10)
    CLOSE(unit=20)

    END Program main

