
! Start of Main Program
! This is a code for solving the 2D Poisson equation using the FDM
! on a parallelogram using coordinate transformation.
! Exercise 2.7

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 101, m = 11, neq = n*m
    INTEGER(4) :: i,j,k
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), PARAMETER :: phi_right = 0.0D0, &
                             & phi_left = 0.0d0, &
                               phi_bot = 1.0d0, &
                               phi_top = 0.0d0, &
                               length = 10.0D0, costh = 0.5d0
    REAL(8), DIMENSION(neq) :: b,xi1,xi2,phi1D,x,y
    REAL(8), DIMENSION(neq,neq) :: a
    REAL(8), DIMENSION(n,m) :: phi2D
    REAL(8) :: dxi1,dxi12,dxi2,dxi22,sinth
    
! Open output file
    OPEN(unit=10, file="Exercise2_7.out", status = "unknown")

    sinth = sqrt(one-costh*costh)
    dxi1 = length/(n-1); dxi2 = one/(m-1)
    dxi12 = dxi1*dxi1; dxi22 = dxi2*dxi2
    
    DO i = 1,n
      xi1(i) = (i-1)*dxi1
    ENDDO
    DO j = 1,m
      xi2(j) = (j-1)*dxi2
    ENDDO


    a(:,:) = zero; b(:) = zero

! Set up coefficient matrix and source vector

! Boundary conditions
    i = 1
    DO j = 2,m-1
      k = (j-1)*n+i
      a(k,k) = one
      b(k) = phi_left
    ENDDO
    i = n
    DO j = 2,m-1
      k = (j-1)*n+i
      a(k,k) = one
      b(k) = phi_right
    ENDDO
    j = 1
    DO i = 1,n
      k = (j-1)*n+i
      a(k,k) = one
      b(k) = phi_bot
    ENDDO
    j = m
    DO i = 1,n
      k = (j-1)*n+i
      a(k,k) = one
      b(k) = phi_top
    ENDDO
    DO i = 2,n-1
      DO j = 2,m-1
        k = (j-1)*n+i
        b(k) = zero
        a(k,k-1) = one/dxi12
        a(k,k+1) = one/dxi12
        a(k,k-N) = one/dxi22
        a(k,k+N) = one/dxi22
        a(k,k+n+1) = -costh/(2.0d0*dxi1*dxi2)
        a(k,k-n-1) = -costh/(2.0d0*dxi1*dxi2)
        a(k,k+n-1) = costh/(2.0d0*dxi1*dxi2)
        a(k,k-n+1) = costh/(2.0d0*dxi1*dxi2)
        a(k,k) = -two/dxi12-two/dxi22
      ENDDO
    ENDDO
     
! Solve system of equations using the Gaussian Elimination

    CALL gauss_elim(neq,a,b,phi1D)

! Print results
    WRITE(10,*) 'VARIABLES = "X", "Y", "PHI"'
    WRITE(10,*) 'ZONE I=', n, ', J=', m, ', DATAPACKING=POINT'
    DO j = 1,m
      DO i = 1,n
         k = (j-1)*n+i
         phi2D(i,j) = phi1D(k)
         x(i) = xi1(i) + xi2(j)*costh
         y(j) = xi2(j)*sinth
         WRITE(*,10) x(i),y(j),phi2D(i,j)
         WRITE(10,10) x(i),y(j),phi2D(i,j)
      ENDDO
    ENDDO

 10 format(3(1x,F13.6))

    CLOSE(unit=10)

    END Program main

!-----------------------------------------------------------------------
  SUBROUTINE gauss_elim(n,a,b,x)
!-----------------------------------------------------------------------
!
! Purpose: Solves a set of linear algebraic equations.
!          Using Gaussian Elimination
!
!          n = number of equations (input)
!          a = nxn coefficient matrix (input)
!          b = right hand vector (input)
!          x = solution vector (output
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   INTEGER(4) :: i,j,k
   INTEGER(4) , INTENT(IN) :: n
   REAL(8) :: sum1,xmult
   REAL(8) :: b(n),a(n,n)
   REAL(8), INTENT(OUT) :: x(n)
!-----------------------------------------------------------------------

   DO k = 1,n-1
     DO i = k+1,n
       xmult = a(i,k)/a(k,k)
       DO j=k+1,n
         a(i,j) = a(i,j)-xmult*a(k,j)
       ENDDO
       a(i,k) = xmult
       b(i) = b(i) - xmult*b(k)
     ENDDO
   ENDDO
   x(n) = b(n)/a(n,n)
   DO i = n-1,1,-1
      sum1 = b(i)
      DO j = i+1,n
         sum1 = sum1 - a(i,j)*x(j)
      ENDDO
      x(i) = sum1/a(i,i)
   ENDDO

   END SUBROUTINE gauss_elim

