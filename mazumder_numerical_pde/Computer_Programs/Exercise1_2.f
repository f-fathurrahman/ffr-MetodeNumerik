! Start of Main Program
! Exercise 1.2

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_nodes = 9, nplot = 1001
    INTEGER(4) :: j,i,il
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0, &
                             & zero= 0.0d0,phi_right=1.0D0, &
                             & length=1.0D0,phi_left = 0.0d0
    REAL(8), DIMENSION(n_nodes) :: x,phi,s,f,acoef
    REAL(8), DIMENSION(n_nodes,n_nodes) :: sm
    REAL(8) :: dx,err,sums,psi,xx
    REAL(8), DIMENSION(nplot) :: phiplot,phi_an
    REAL(8) :: flux_l_ana,flux_r_ana,sum_source_ana,imb_ana
    REAL(8) :: flux_l_num,flux_r_num,sum_source_num,imb_num
    
! Open output file
    OPEN(unit=10, file="Example1_2.out", status = "unknown")
    
    dx = length/(n_nodes-1)
    DO i = 1,n_nodes
      x(i) = (i-1)*dx
      s(i) = cos(x(i))
    ENDDO

    sm(:,:) = zero; f(:) = zero

! Set up stiffness matrix (sm) and load vector (f)

    f(1) = zero   ! Recover a(1) = 0
    sm(1,1) = one
    DO i = 2,n_nodes-1
      f(i) = -(two*s(i)-s(i-1)-s(i+1))/dx
      sm(i,i-1) = -one/dx
      sm(i,i+1) = -one/dx
      sm(i,i) = two/dx
    ENDDO
    f(n_nodes) = one  ! Recover a(n) = 1
    sm(n_nodes,n_nodes) = one

! Solve system of equations using the Gaussian Elimination to determine coefficients

    CALL gauss_elim(n_nodes,sm,f,acoef)

! Compute Solution (phiplot)
    DO i = 2,nplot-1
      xx = (i-1)*one/(nplot-1)
      sums = zero
      DO j = 2, n_nodes-1 
        IF((xx .LE. x(j)) .AND. (xx .GE. x(j-1)))THEN
          psi = (xx - x(j-1))/dx
        ELSEIF((xx .LE. x(j+1)) .AND. (xx .GE. x(j)))THEN
          psi = (x(j+1) - xx)/dx
        ELSE
          psi = zero
        ENDIF
        sums = sums + acoef(j)*psi
      ENDDO
      IF((xx .LE. x(n_nodes)) .AND. (xx .GE. x(n_nodes-1)))THEN
        psi = (xx - x(n_nodes-1))/dx
      ELSE
        psi = zero
      ENDIF
      sums = sums + acoef(n_nodes)*psi
      phiplot(i) = sums
    ENDDO

! Print results for phi
    WRITE(*,10) x(1),phi_left,phi_left,zero
    WRITE(10,10) x(1),phi_left,phi_left,zero
    DO i = 2,nplot-1
      xx = (i-1)*one/(nplot-1)
      phi_an(i) = one-cos(xx)+cos(one)*xx      
      err = (phi_an(i)-phiplot(i))
      WRITE(*,10) xx,phi_an(i),phiplot(i),err
      WRITE(10,10) xx,phi_an(i),phiplot(i),err
    ENDDO
    WRITE(*,10) x(n_nodes),phi_right,phi_right,zero
    WRITE(10,10) x(n_nodes),phi_right,phi_right,zero

! Print results for fluxes
    flux_l_ana = cos(one); flux_r_ana = sin(one)+cos(one); sum_source_ana = sin(one)
    imb_ana = flux_r_ana-flux_l_ana-sum_source_ana
    flux_l_num = acoef(2)/dx; flux_r_num = (one-acoef(n_nodes-1))/dx; sum_source_num = sin(one)
    imb_num = flux_r_num-flux_l_num-sum_source_num
    WRITE(*,11) flux_l_ana,flux_r_ana,sum_source_ana,imb_ana
    WRITE(*,11) flux_l_num,flux_r_num,sum_source_num,imb_num

 10 format(3(F14.7),1x,E14.6)
 11 format(4(F14.6))

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

