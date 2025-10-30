! Start of Main Program
! This is a program to solve a second order linear elliptic PDE
! using the finite-difference method and the geometric multi-grid (GMG)
! method. Only two levels are used. The smoother is Gauss-Seidel.
! Example 4.9

    program main

    USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: nff=161, mff=161, max_iter=100000
    INTEGER(4), PARAMETER :: ncc=81, mcc=81
    INTEGER(4) :: iter,i,j,iff,jff,icc,jcc,k
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,pi = 3.1415927d0, &
                             & zero= 0.0d0,length=1.0D0,tol=2.0D-8, &
                             & four=4.0d0, width=1.0d0,tiny=1.0d-40
    REAL(8), PARAMETER :: aa = 0.5d0, bb = 0.5d0, cc = 10.0d0, &
                               dd = 1.0d0, ee = 2.0d0
    REAL(8), DIMENSION(nff) :: xff
    REAL(8), DIMENSION(mff) :: yff
    REAL(8), DIMENSION(nff,mff) :: phiff,phia,rff,errff,sij,phie
    REAL(8), DIMENSION(ncc) :: xcc
    REAL(8), DIMENSION(mcc) :: ycc
    REAL(8), DIMENSION(ncc,mcc) :: phicc,rcc,errcc
    REAL(8) :: dxff,dx2ff,dyff,dy2ff,err_i,sumres,res,fde,diag
    REAL(8) :: dxcc,dx2cc,dycc,dy2cc
    REAL(8) :: errmin,errmax
    REAL(4) :: ta(2), time_start, time_end
    LOGICAL :: second
    

! Open output file
    OPEN(unit=10, file="Two_grid.out", status = "unknown")
    OPEN(unit=20, file="Two_grid.rsl", status = "unknown")
    OPEN(unit=30, file="Two_grid_e1.out", status = "unknown")
    OPEN(unit=40, file="Two_grid_e2.out", status = "unknown")
    WRITE(30,*) 'VARIABLES = "X", "Y", "Error"'
    WRITE(30,*) 'ZONE I=', nff, ', J=', mff, ', DATAPACKING=POINT'
    WRITE(40,*) 'VARIABLES = "X", "Y", "Error"'
    WRITE(40,*) 'ZONE I=', nff, ', J=', mff, ', DATAPACKING=POINT'

    second = .false.

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

 12 CONTINUE

    phicc(:,:) = zero   ! Initial Guess for phi in coarse mesh
    phiff(:,:) = zero   ! Initial Guess for phi in fine mesh

! Boundary Conditions for coarse mesh
    DO i = 1,ncc
      phicc(i,1) = bb*bb*sinh(-cc*bb)+(xcc(i)-aa)*(xcc(i)-aa)*sinh(cc*(xcc(i)-aa))+dd !bottom
      phicc(i,mcc) = (one-bb)*(one-bb)*sinh(cc*(one-bb))+ &
                     (xcc(i)-aa)*(xcc(i)-aa)*sinh(cc*(xcc(i)-aa))+dd*exp(ee*xcc(i))  !top
    ENDDO
    DO j = 1,mcc
      phicc(1,j) = aa*aa*sinh(-cc*aa)+(ycc(j)-bb)*(ycc(j)-bb)*sinh(cc*(ycc(j)-bb))+dd !left
      phicc(ncc,j) = (one-aa)*(one-aa)*sinh(cc*(one-aa))+ &
                     (ycc(j)-bb)*(ycc(j)-bb)*sinh(cc*(ycc(j)-bb))+dd*exp(ee*ycc(j)) !right
    ENDDO

! Boundary Conditions for fine mesh
    DO i = 1,nff
      phiff(i,1) = bb*bb*sinh(-cc*bb)+(xff(i)-aa)*(xff(i)-aa)*sinh(cc*(xff(i)-aa))+dd !bottom
      phiff(i,mff) = (one-bb)*(one-bb)*sinh(cc*(one-bb))+ &
                     (xff(i)-aa)*(xff(i)-aa)*sinh(cc*(xff(i)-aa))+dd*exp(ee*xff(i))  !top
    ENDDO
    DO j = 1,mff
      phiff(1,j) = aa*aa*sinh(-cc*aa)+(yff(j)-bb)*(yff(j)-bb)*sinh(cc*(yff(j)-bb))+dd !left
      phiff(nff,j) = (one-aa)*(one-aa)*sinh(cc*(one-aa))+ &
                     (yff(j)-bb)*(yff(j)-bb)*sinh(cc*(yff(j)-bb))+dd*exp(ee*yff(j)) !right
    ENDDO

! Source term calculation: Only fine mesh needed
    DO i = 1,nff
      DO j = 1,mff
        sij(i,j) = 2.0d0*sinh(cc*(xff(i)-aa))+4.0d0*(xff(i)-aa)*cc*cosh(cc*(xff(i)-aa))+ &
                 (xff(i)-aa)*(xff(i)-aa)*(cc*cc)*sinh(cc*(xff(i)-aa)) + &
                 2.0d0*sinh(cc*(yff(j)-bb))+4.0d0*(yff(j)-bb)*cc*cosh(cc*(yff(j)-bb))+ &
                 (yff(j)-bb)*(yff(j)-bb)*(cc*cc)*sinh(cc*(yff(j)-bb)) + &
                 dd*ee*ee*(yff(j)*yff(j)+xff(i)*xff(i))*exp(ee*xff(i)*yff(j))
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
    
    IF(second .and. iter == 1)THEN
      errmin = 1.0D10
      errmax = -1.0D10
      DO j = 1,mff
        DO i = 1,nff
          errmin = MIN(errmin,phie(i,j)-phiff(i,j))
          errmax = MAX(errmax,phie(i,j)-phiff(i,j))
          WRITE(30,*) xff(i),yff(j),phie(i,j)-phiff(i,j)
        ENDDO
      ENDDO
      print *, errmin, errmax
    ENDIF

! Calculate Residual for fine mesh
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

! Gauss-Seidel Smoothing for Coarse Mesh (Equation in Correction form)
    diag = two/dx2cc + two/dy2cc
    errcc(:,:) = zero
    DO k = 1,3
    DO j = 2,mcc-1
      DO i = 2,ncc-1
        !errcc(i,j) = rcc(i,j)/diag
        errcc(i,j) = (rcc(i,j)+(errcc(i+1,j)+errcc(i-1,j))/dx2cc+(errcc(i,j+1)+errcc(i,j-1))/dy2cc)/diag
      ENDDO
    ENDDO
    ENDDO

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

    IF(second .and. iter == 1)THEN
      errmin = 1.0D10
      errmax = -1.0D10
      DO j = 1,mff
        DO i = 1,nff
          errmin = MIN(errmin,phie(i,j)-phiff(i,j))
          errmax = MAX(errmax,phie(i,j)-phiff(i,j))
          WRITE(40,*) xff(i),yff(j),phie(i,j)-phiff(i,j)
        ENDDO
      ENDDO
      print *, errmin, errmax
      STOP
    ENDIF

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

    phie(:,:) = phiff(:,:)
    second = .true.
    GO TO 12

    time_end = etime(ta)

    print *, "CPU time taken =", time_end-time_start

! Print results for fine mesh
    WRITE(10,*) 'VARIABLES = "X", "Y", "PHIAN", "PHI", "Error"'
    WRITE(10,*) 'ZONE I=', nff, ', J=', mff, ', DATAPACKING=POINT'
    DO j = 1,mff
      DO i = 1,nff
        phia(i,j) = ((xff(i)-aa)**2)*sinh(cc*(xff(i)-aa)) + &
                    ((yff(j)-bb)**2)*sinh(cc*(yff(j)-bb)) + dd*exp(ee*xff(i)*yff(j))
        WRITE(10,10) xff(i),yff(j),phia(i,j),phiff(i,j),phia(i,j)-phiff(i,j)
      ENDDO
    ENDDO

 10 format(5(1x,F13.6))

    CLOSE(unit=10)
    CLOSE(unit=20)

    END Program main

