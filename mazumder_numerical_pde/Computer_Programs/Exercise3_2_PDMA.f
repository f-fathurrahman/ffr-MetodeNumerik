! Start of Main Program
! ADI Method with Penta-Diagonal Matrix (PDMA) Solver
! Exercise 3.2

    program main

    USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 81, m = 81
    INTEGER(4), PARAMETER :: max_iter = 100000
    INTEGER(4) :: i,j,iter
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), PARAMETER :: aa = 0.5d0, bb = 0.5d0, cc = 10.0d0, &
                               dd = 1.0d0, ee = 2.0d0
    REAL(8), PARAMETER :: length = 1.0D0
    REAL(8) :: tol = 1.0d-6
    REAL(8), DIMENSION(n) :: x,phi_bot,phi_top
    REAL(8), DIMENSION(n) :: d1,a1,c1,b1,sol1,e1,f1
    REAL(8), DIMENSION(m) :: y,phi_left,phi_right
    REAL(8), DIMENSION(m) :: d2,a2,c2,b2,sol2,e2,f2
    REAL(8), DIMENSION(n,m) :: S,phi,phian
    REAL(8) :: dx,dx2,dy,dy2
    REAL(8) :: ae,aw,as,an,ao,source,res,fde,aee,aww,ass,ann
    REAL(4) :: ta(2),time_start,time_end
    

! Open output file
    OPEN(unit=10, file="Exercise3_2_PDMA.out", status = "unknown")
    OPEN(unit=20, file="Exercise3_2_PDMA.rsl", status = "unknown")

    dx = length/(n-1); dy = length/(m-1)
    dx2 = dx*dx; dy2 = dy*dy
    
    DO i = 1,n
      x(i) = (i-1)*dx
      phi_bot(i) = bb*bb*sinh(-cc*bb)+(x(i)-aa)*(x(i)-aa)*sinh(cc*(x(i)-aa))+dd
      phi_top(i) = (one-bb)*(one-bb)*sinh(cc*(one-bb))+(x(i)-aa)*(x(i)-aa)*sinh(cc*(x(i)-aa))+dd*exp(ee*x(i))
    ENDDO
    DO j = 1,m
      y(j) = (j-1)*dy
      phi_left(j) = aa*aa*sinh(-cc*aa)+(y(j)-bb)*(y(j)-bb)*sinh(cc*(y(j)-bb))+dd
      phi_right(j) = (one-aa)*(one-aa)*sinh(cc*(one-aa))+(y(j)-bb)*(y(j)-bb)*sinh(cc*(y(j)-bb))+dd*exp(ee*y(j))
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


! Link Coefficients (for interior rows/columns)
    aw = -4.0d0/(3.0d0*dx2)
    ae = -4.0d0/(3.0d0*dx2)
    as = -4.0d0/(3.0d0*dy2)
    an = -4.0d0/(3.0d0*dy2)
    ao = 2.5d0/dx2+2.5d0/dy2
    aww = one/(12.0d0*dx2)
    aee = one/(12.0d0*dx2)
    ass = one/(12.0d0*dy2)
    ann = one/(12.0d0*dy2)

! Start of iteration loop

    time_start = etime(ta)
    
    DO iter = 1,max_iter

! Row-wise sweep
      j = 2   ! Second Row
      d1(:) = zero; c1(:) = zero; a1(:) = zero; b1(:) = zero; e1(:) = zero; f1(:) = zero
      d1(1) = one
      b1(1) = phi_left(j)
      DO i = 2,n-1
        d1(i) = 2.0d0/dx2+2.0d0/dy2
        c1(i) = -1.0d0/dx2
        a1(i-1) = -1.0d0/dx2
        b1(i) = -S(i,j)+(phi(i,j+1)+phi(i,j-1))/dy2
      ENDDO
      d1(n) = one
      b1(n) = phi_right(j)

      CALL PENTA(n,e1,a1,d1,c1,f1,b1,sol1)

      phi(:,j) = sol1(:)

      DO j = 3,m-2   ! Interior rows
        
        d1(:) = zero; c1(:) = zero; a1(:) = zero; b1(:) = zero; e1(:) = zero; f1(:) = zero
        d1(1) = one
        b1(1) = phi_left(j)
        d1(2) = 2.0d0/dx2+2.0d0/dy2
        c1(2) = -1.0d0/dx2
        a1(1) = -1.0d0/dx2
        b1(2) = -S(2,j)+(phi(2,j+1)+phi(2,j-1))/dy2
        DO i = 3,n-2
          d1(i) = ao
          c1(i) = ae
          a1(i-1) = aw
          e1(i-2) = aww
          f1(i) = aee
          b1(i) = -S(i,j)-an*phi(i,j+1)-as*phi(i,j-1) &
                  -(phi(i,j+2)+phi(i,j-2))/(12.0d0*dy2)
        ENDDO
        d1(n-1) = 2.0d0/dx2+2.0d0/dy2
        c1(n-1) = -1.0d0/dx2
        a1(n-2) = -1.0d0/dx2
        b1(n-1) = -S(n-1,j)+(phi(n-1,j+1)+phi(n-1,j-1))/dy2
        d1(n) = one
        b1(n) = phi_right(j)

        CALL PENTA(n,e1,a1,d1,c1,f1,b1,sol1)

        phi(:,j) = sol1(:)

      ENDDO

      j = m-1   ! Second-last Row
      d1(:) = zero; c1(:) = zero; a1(:) = zero; b1(:) = zero; e1(:) = zero; f1(:) = zero
      d1(1) = one
      b1(1) = phi_left(j)
      DO i = 2,n-1
        d1(i) = 2.0d0/dx2+2.0d0/dy2
        c1(i) = -1.0d0/dx2
        a1(i-1) = -1.0d0/dx2
        b1(i) = -S(i,j)+(phi(i,j+1)+phi(i,j-1))/dy2
      ENDDO
      d1(n) = one
      b1(n) = phi_right(j)

      CALL PENTA(n,e1,a1,d1,c1,f1,b1,sol1)

      phi(:,j) = sol1(:)

! Column-wise sweep
      i = 2   ! Second Column
      d2(:) = zero; c2(:) = zero; a2(:) = zero; b2(:) = zero; e2(:) = zero; f2(:) = zero
      d2(1) = one
      b2(1) = phi_bot(i)
      DO j = 2,m-1
        d2(j) = 2.0d0/dx2+2.0d0/dy2
        c2(j) = -1.0d0/dy2
        a2(j-1) = -1.0d0/dy2
        b2(j) = -S(i,j)+(phi(i+1,j)+phi(i-1,j))/dx2
      ENDDO
      d2(m) = one
      b2(m) = phi_top(i)

      CALL PENTA(m,e2,a2,d2,c2,f2,b2,sol2)

      phi(i,:) = sol2(:)

      DO i = 3,n-2   ! Interior columns
        
        d2(:) = zero; c2(:) = zero; a2(:) = zero; b2(:) = zero; e2(:) = zero; f2(:) = zero
        d2(1) = one
        b2(1) = phi_bot(i)
        d2(2) = 2.0d0/dx2+2.0d0/dy2
        c2(2) = -1.0d0/dy2
        a2(1) = -1.0d0/dy2
        b2(2) = -S(i,2)+(phi(i+1,2)+phi(i-1,2))/dx2
        DO j = 3,m-2
          d2(j) = ao
          c2(j) = as
          a2(j-1) = an
          e2(j-2) = ass
          f2(j) = ann
          b2(j) = -S(i,j)-ae*phi(i+1,j)-aw*phi(i-1,j) &
                  -(phi(i+2,j)+phi(i-2,j))/(12.0d0*dx2)
        ENDDO
        d2(m-1) = 2.0d0/dx2+2.0d0/dy2
        c2(m-1) = -1.0d0/dy2
        a2(m-2) = -1.0d0/dy2
        b2(m-1) = -S(i,m-1)+(phi(i+1,m-1)+phi(i-1,m-1))/dx2
        d2(m) = one
        b2(m) = phi_top(i)

        CALL PENTA(m,e2,a2,d2,c2,f2,b2,sol2)

        phi(i,:) = sol2(:)

      ENDDO

      i = n-1   ! Second-last Column
      d2(:) = zero; c2(:) = zero; a2(:) = zero; b2(:) = zero; e2(:) = zero; f2(:) = zero
      d2(1) = one
      b2(1) = phi_bot(i)
      DO j = 2,m-1
        d2(j) = 2.0d0/dx2+2.0d0/dy2
        c2(j) = -1.0d0/dy2
        a2(j-1) = -1.0d0/dy2
        b2(j) = -S(i,j)+(phi(i+1,j)+phi(i-1,j))/dx2
      ENDDO
      d2(m) = one
      b2(m) = phi_top(i)

      CALL PENTA(m,e2,a2,d2,c2,f2,b2,sol2)

      phi(i,:) = sol2(:)


! Residual calculation
      res = zero
! Interior nodes      
      DO i = 3,n-2
        DO j = 3,m-2
          source = -S(i,j)
          fde = source-aw*phi(i-1,j)-ae*phi(i+1,j)-as*phi(i,j-1)-an*phi(i,j+1)-ao*phi(i,j) &
                      -(phi(i+2,j)+phi(i-2,j))/(12.0d0*dx2)-(phi(i,j+2)+phi(i,j-2))/(12.0d0*dy2)
          res = res + fde*fde
        ENDDO
      ENDDO
! Bottom Boundary
      j = 2
      DO i = 2,n-1
          source = -S(i,j)
          fde = source+(phi(i-1,j)+phi(i+1,j))/dx2+(phi(i,j-1)+phi(i,j+1))/dy2-phi(i,j)*(2.0d0/dx2+2.0d0/dy2)
          res = res + fde*fde
      ENDDO
! Bottom Boundary
      j = m-1
      DO i = 2,n-1
          source = -S(i,j)
          fde = source+(phi(i-1,j)+phi(i+1,j))/dx2+(phi(i,j-1)+phi(i,j+1))/dy2-phi(i,j)*(2.0d0/dx2+2.0d0/dy2)
          res = res + fde*fde
      ENDDO
! Left Boundary
      i = 2
      DO j = 3,m-2
          source = -S(i,j)
          fde = source+(phi(i-1,j)+phi(i+1,j))/dx2+(phi(i,j-1)+phi(i,j+1))/dy2-phi(i,j)*(2.0d0/dx2+2.0d0/dy2)
          res = res + fde*fde
      ENDDO
! Right Boundary
      i = n-1
      DO j = 3,m-2
          source = -S(i,j)
          fde = source+(phi(i-1,j)+phi(i+1,j))/dx2+(phi(i,j-1)+phi(i,j+1))/dy2-phi(i,j)*(2.0d0/dx2+2.0d0/dy2)
          res = res + fde*fde
      ENDDO


      res = SQRT(MAX(zero,res))
      WRITE(*,*) iter*2,res
      WRITE(20,*) iter*2,res

      IF(res < tol) EXIT

    ENDDO  !iterations
  
    time_end = etime(ta)
    
! Print results
    WRITE(10,*) 'VARIABLES = "X", "Y", "PHIAN", "PHI", "Error"'
    WRITE(10,*) 'ZONE I=', n, ', J=', m, ', DATAPACKING=POINT'
    DO j = 1,m
      DO i = 1,n
         phian(i,j) = ((x(i)-aa)**2)*sinh(cc*(x(i)-aa)) + ((y(j)-bb)**2)*sinh(cc*(y(j)-bb))+dd*exp(ee*x(i)*y(j))
         WRITE(10,10) x(i),y(j),phian(i,j),phi(i,j),phian(i,j)-phi(i,j)
      ENDDO
    ENDDO

    print *, "time taken = ", time_end-time_start

 10 format(5(1x,F13.6))

    CLOSE(unit=10)
    CLOSE(unit=20)

    END Program main

!--------------------------------------------------------------------------

! Pentadiagonal solver

    SUBROUTINE PENTA(n,e,a,d,c,f,b,x)

    INTEGER(4), INTENT(in) :: n
    INTEGER(4) :: i
    REAL(8) :: xmult
    REAL(8), DIMENSION(n) :: a,b,c,d,e,f,x

    DO i = 2,n-1
      xmult = a(i-1)/d(i-1)
      d(i) = d(i) - xmult*c(i-1)
      c(i) = c(i) - xmult*f(i-1)
      b(i) = b(i) - xmult*b(i-1)
      xmult = e(i-1)/d(i-1)
      a(i) = a(i) - xmult*c(i-1)
      d(i+1) = d(i+1) - xmult*f(i-1)
      b(i+1) = b(i+1) - xmult*b(i-1)
    ENDDO
    
    xmult = a(n-1)/d(n-1)
    d(n) = d(n) - xmult*c(n-1)
    x(n) = (b(n) - xmult*b(n-1))/d(n)
    x(n-1) = (b(n-1)-c(n-1)*x(n))/d(n-1)
    DO i = n-2,1,-1
      x(i) = (b(i)-c(i)*x(i+1)-f(i)*x(i+2))/d(i)
    ENDDO

    END SUBROUTINE PENTA


