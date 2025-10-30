! Start of Main Program
! This is a program to solve a second order linear elliptic PDE
! using the finite-difference method and the geometric multi-grid (GMG)
! method. Only 2 grid levels are used. The smoother is Gauss-Seidel.
! Exercise 4.4

! method = 1 => Gauss-Seidel Smoothing on coarse mesh
! method = 2 => Gauusian Elimination on coarse mesh

    program main

    USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: nff=81, mff=81
    INTEGER(4), PARAMETER :: ncc=41, mcc=41
    INTEGER(4), PARAMETER :: neq = ncc*mcc
    INTEGER(4), PARAMETER :: max_iter=100000, max_iter_c = 1, method = 2
    INTEGER(4) :: iter,i,j,iff,jff,icc,jcc,k
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,pi = 3.1415927d0, &
                             & zero= 0.0d0,length=1.0D0,tol=1.0D-6, &
                             & four=4.0d0, width=1.0d0,tiny=1.0d-40
    REAL(8), DIMENSION(nff) :: xff
    REAL(8), DIMENSION(mff) :: yff
    REAL(8), DIMENSION(nff,mff) :: phiff,phia,rff,errff,sij,phie
    REAL(8), DIMENSION(ncc) :: xcc
    REAL(8), DIMENSION(mcc) :: ycc
    REAL(8), DIMENSION(ncc,mcc) :: phicc,rcc,errcc
    REAL(8) :: acoef(neq,neq),sol(neq),rhs(neq)
    REAL(8) :: dxff,dx2ff,dyff,dy2ff,err_i,sumres,res,fde,diag
    REAL(8) :: dxcc,dx2cc,dycc,dy2cc
    REAL(8) :: ae,aw,an,as
    REAL(4) :: ta(2), time_start, time_end
    

! Open output file
    OPEN(unit=10, file="Exercise4_4.out", status = "unknown")
    OPEN(unit=20, file="Exercise4_4.rsl", status = "unknown")

! Set up Coarse Mesh
    dxcc = length/(ncc-1)
    dycc = width/(mcc-1)
    dx2cc = dxcc*dxcc
    dy2cc = dycc*dycc
    DO i = 1,ncc
      xcc(i) = (i-1)*dxcc
    ENDDO
    DO j = 1,mcc
      ycc(j) = (j-1)*dycc
    ENDDO

! Set up Fine Mesh
    dxff = length/(nff-1)
    dyff = width/(mff-1)
    dx2ff = dxff*dxff
    dy2ff = dyff*dyff
    DO i = 1,nff
      xff(i) = (i-1)*dxff
    ENDDO
    DO j = 1,mff
      yff(j) = (j-1)*dyff
    ENDDO

    phicc(:,:) = zero   ! Initial Guess for phi in coarse mesh
    phiff(:,:) = zero   ! Initial Guess for phi in fine mesh

! Boundary Conditions for coarse mesh
    DO i = 1,ncc
      phicc(i,1) = 100.0d0*xcc(i)+500.0d0*exp(-50.0*(one-xcc(i))**2) !bottom
      phicc(i,mcc) = 500.0d0*exp(-50.0*((one-xcc(i))**2+one)) !top
    ENDDO
    DO j = 1,mcc
      phicc(1,j) = 500.0d0*exp(-50.0*(one+ycc(j)**2)) !left
      phicc(ncc,j) = 100*(one-ycc(j))+500.0d0*exp(-50.0*ycc(j)**2) !right
    ENDDO

! Boundary Conditions for fine mesh
    DO i = 1,nff
      phiff(i,1) = 100.0d0*xff(i)+500.0d0*exp(-50.0*(one-xff(i))**2) !bottom
      phiff(i,mff) = 500.0d0*exp(-50.0*((one-xff(i))**2+one)) !top
    ENDDO
    DO j = 1,mff
      phiff(1,j) = 500.0d0*exp(-50.0*(one+yff(j)**2)) !left
      phiff(nff,j) = 100*(one-yff(j))+500.0d0*exp(-50.0*yff(j)**2) !right
    ENDDO

! Source term calculation: Only fine mesh needed
    DO i = 1,nff
      DO j = 1,mff
        sij(i,j) = 5.0D4*(100.0d0*((one-xff(i))**2+yff(j)**2)-two)* &
                  exp(-50.0d0*((one-xff(i))**2+yff(j)**2))
      ENDDO
    ENDDO

    iter = 0

    time_start = etime(ta)

 11 iter = iter + 1   ! Start iteration

! Gauss-Seidel Smoothing for Fine Mesh
    diag = two/dx2ff + two/dy2ff
    DO j = 2,mff-1
      DO i = 2,nff-1
        phiff(i,j) = ((phiff(i-1,j)+phiff(i+1,j))/dx2ff + &
                   (phiff(i,j-1)+phiff(i,j+1))/dy2ff - &
                   sij(i,j))/diag
      ENDDO
    ENDDO
    
! Calculate Residual for fine mesh. Default of 1 sweep
    rff(:,:) = zero ! Sets residuals to zero at boundaries
    DO i = 2,nff-1
      DO j = 2,mff-1
        rff(i,j) = (phiff(i+1,j)-two*phiff(i,j)+phiff(i-1,j))/dx2ff + &
                   (phiff(i,j+1)-two*phiff(i,j)+phiff(i,j-1))/dy2ff - &
                   sij(i,j)
      ENDDO
    ENDDO

! Obtain coarse mesh residuals from fine mesh (Restriction)
    rcc(:,:) = zero !Sets residuals to zero at boundaries
    DO icc = 2,ncc-1
      iff = 2*icc-1
      DO jcc = 2,mcc-1
        jff = 2*jcc-1
        rcc(icc,jcc) = rff(iff,jff)
      ENDDO
    ENDDO

    IF (method == 1) THEN

! Gauss-Seidel Smoothing for Coarse Mesh (Equation in Correction form)
      diag = two/dx2cc + two/dy2cc
      errcc(:,:) = zero
      DO k = 1,max_iter_c  ! Coarse mesh sweeps
        DO j = 2,mcc-1
          DO i = 2,ncc-1
            errcc(i,j) = (rcc(i,j)+(errcc(i+1,j)+errcc(i-1,j))/dx2cc+(errcc(i,j+1)+errcc(i,j-1))/dy2cc)/diag
          ENDDO
        ENDDO
      ENDDO
    
    ELSE  

! Gaussian Elimination
      acoef(:,:) = zero; sol(:) = zero; rhs(:) = zero
      diag = two/dx2cc + two/dy2cc
      ae = -one/dx2cc
      aw = -one/dx2cc
      an = -one/dy2cc
      as = -one/dy2cc
! Interior Nodes
      DO i = 2,ncc-1
        DO j = 2,mcc-1
          k = (j-1)*ncc+i
          acoef(k,k) = diag
          acoef(k,k+1) = ae
          acoef(k,k-1) = aw
          acoef(k,k+ncc) = an
          acoef(k,k-ncc) = as
          rhs(k) = rcc(i,j)
        ENDDO
      ENDDO
! Left Boundary
      i = 1
      DO j = 1,mcc
        k = (j-1)*ncc+i
        acoef(k,k) = one
        rhs(k) = zero
      ENDDO
! Right Boundary
      i = ncc
      DO j = 1,mcc
        k = (j-1)*ncc+i
        acoef(k,k) = one
        rhs(k) = zero
      ENDDO
! Bottom Boundary
      j = 1
      DO i = 1,ncc
        k = (j-1)*ncc+i
        acoef(k,k) = one
        rhs(k) = zero
      ENDDO
! Top Boundary
      j = mcc
      DO i = 1,ncc
        k = (j-1)*ncc+i
        acoef(k,k) = one
        rhs(k) = zero
      ENDDO


      CALL gauss_elim(neq,acoef,rhs,sol)
   
! Transfer correction to 2D array
      DO i = 1,ncc
        DO j = 1,mcc
          k = (j-1)*ncc+i
          errcc(i,j) = sol(k)
        ENDDO
      ENDDO

    ENDIF

! Obtain fine mesh correction from coarse mesh (Prolongation)

    errff(:,:) = zero
    DO icc = 2,ncc-1
      iff = 2*icc-1
      DO jcc = 2,mcc-1
        jff = 2*jcc-1
        errff(iff,jff) = errcc(icc,jcc)
        errff(iff-1,jff) = (errcc(icc,jcc)+errcc(icc-1,jcc))/two
        errff(iff,jff-1) = (errcc(icc,jcc)+errcc(icc,jcc-1))/two
        errff(iff-1,jff-1) = (errcc(icc,jcc)+errcc(icc-1,jcc) + &
                              errcc(icc,jcc-1)+errcc(icc-1,jcc-1))/four
      ENDDO
    ENDDO

! Update fine mesh solution
    phiff(:,:) = phiff(:,:) + errff(:,:)

!   Calculate L2NORM on fine mesh
    sumres = zero
    DO j = 2,mff-1
      DO i = 2,nff-1
        fde = rff(i,j)
        sumres = sumres + fde*fde
      ENDDO
    ENDDO
    res = sqrt(MAX(tiny,sumres))

    print *, iter, res
    WRITE(20,*) iter, res

    IF(res > tol .AND. iter < max_iter) GO TO 11    ! Convergence check
    
    time_end = etime(ta)

    print *, "CPU time taken =", time_end-time_start

! Print results for fine mesh
    WRITE(10,*) 'VARIABLES = "X", "Y", "PHIAN", "PHI", "Error"'
    WRITE(10,*) 'ZONE I=', nff, ', J=', mff, ', DATAPACKING=POINT'
    DO j = 1,mff
      DO i = 1,nff
        phia(i,j) = 500.0d0*exp(-50.0*((one-xff(i))**2+yff(j)**2)) + 100.0d0*xff(i)*(one-yff(j))
        WRITE(10,10) xff(i),yff(j),phia(i,j),phiff(i,j),phia(i,j)-phiff(i,j)
      ENDDO
    ENDDO

 10 format(5(1x,F13.6))

    CLOSE(unit=10)
    CLOSE(unit=20)

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


