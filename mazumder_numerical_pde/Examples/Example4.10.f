! Start of Main Program
! This is a program to solve a second order linear elliptic PDE
! using the finite-difference method and the V-cycle geometric multi-grid (GMG)
! method. Number of levels = "lev". The smoother is Gauss-Seidel.

    program main

    USE DFPORT

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: lev = 5
    INTEGER(4), PARAMETER :: max_iter=100000
    INTEGER(4), PARAMETER :: nff=161, mff=161
    INTEGER(4) :: n(lev), m(lev)
    INTEGER(4) :: iter,i,j,iff,jff,l,k,kmax
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,pi = 3.1415927d0, &
                             & zero= 0.0d0,length=1.0D0,tol=1.0D-6, &
                             & four=4.0d0, width=1.0d0,tiny=1.0d-40
    REAL(8), PARAMETER :: aa = 0.5d0, bb = 0.5d0, cc = 10.0d0, &
                               dd = 1.0d0, ee = 2.0d0
    REAL(8), DIMENSION(nff,lev) :: x
    REAL(8), DIMENSION(mff,lev) :: y
    REAL(8), DIMENSION(lev) :: dx,dy,dx2,dy2
    REAL(8), DIMENSION(nff,mff,lev) :: res,err
    REAL(8), DIMENSION(nff,mff) :: phi,phia,sij
    REAL(8) :: sumres,res2,fde,diag
    REAL(4) :: ta(2), time_start, time_end
    LOGICAL :: second
    

! Open output file
    OPEN(unit=10, file="multi_grid.out", status = "unknown")
    OPEN(unit=20, file="multi_grid.rsl", status = "unknown")

! Number of grid points at each grid level
    n(1) = nff; m(1) = mff
    DO l = 2,lev
      n(l) = (n(l-1)-1)/2 + 1
      m(l) = (m(l-1)-1)/2 + 1
    ENDDO
   
! Set up Grids
    DO l = 1,lev
      dx(l) = length/(n(l)-1)
      dy(l) = width/(m(l)-1)
      dx2(l) = dx(l)*dx(l)
      dy2(l) = dy(l)*dy(l)
      DO i = 1,n(l)
        x(i,l) = (i-1)*dx(l)
      ENDDO
      DO j = 1,m(l)
        y(j,l) = (j-1)*dy(l)
      ENDDO
    ENDDO
   
    phi(:,:) = zero   ! Initial Guess for phi in fine mesh

! Boundary Conditions for fine mesh
    DO i = 1,nff
      phi(i,1) = bb*bb*sinh(-cc*bb)+(x(i,1)-aa)*(x(i,1)-aa)*sinh(cc*(x(i,1)-aa))+dd !bottom
      phi(i,mff) = (one-bb)*(one-bb)*sinh(cc*(one-bb))+ &
                     (x(i,1)-aa)*(x(i,1)-aa)*sinh(cc*(x(i,1)-aa))+dd*exp(ee*x(i,1))  !top
    ENDDO
    DO j = 1,mff
      phi(1,j) = aa*aa*sinh(-cc*aa)+(y(j,1)-bb)*(y(j,1)-bb)*sinh(cc*(y(j,1)-bb))+dd !left
      phi(nff,j) = (one-aa)*(one-aa)*sinh(cc*(one-aa))+ &
                     (y(j,1)-bb)*(y(j,1)-bb)*sinh(cc*(y(j,1)-bb))+dd*exp(ee*y(j,1)) !right
    ENDDO

! Source term calculation: Only fine mesh needed
    DO i = 1,nff
      DO j = 1,mff
        sij(i,j) = 2.0d0*sinh(cc*(x(i,1)-aa))+4.0d0*(x(i,1)-aa)*cc*cosh(cc*(x(i,1)-aa))+ &
                   (x(i,1)-aa)*(x(i,1)-aa)*(cc*cc)*sinh(cc*(x(i,1)-aa)) + &
                   2.0d0*sinh(cc*(y(j,1)-bb))+4.0d0*(y(j,1)-bb)*cc*cosh(cc*(y(j,1)-bb))+ &
                   (y(j,1)-bb)*(y(j,1)-bb)*(cc*cc)*sinh(cc*(y(j,1)-bb)) + &
                   dd*ee*ee*(y(j,1)*y(j,1)+x(i,1)*x(i,1))*exp(ee*x(i,1)*y(j,1))
      ENDDO
    ENDDO

    iter = 0

    time_start = etime(ta)

 11 iter = iter + 1   ! Start iteration

! Gauss-Seidel Smoothing for Finest Mesh
    diag = two/dx2(1) + two/dy2(1)
    DO j = 2,mff-1
      DO i = 2,nff-1
        phi(i,j) = ((phi(i-1,j)+phi(i+1,j))/dx2(1) + &
                   (phi(i,j-1)+phi(i,j+1))/dy2(1) - &
                   sij(i,j))/diag
      ENDDO
    ENDDO
    
 
! Calculate Residual for fine mesh
    res(:,:,1) = zero ! Sets residuals to zero at boundaries
    DO i = 2,nff-1
      DO j = 2,mff-1
        res(i,j,1) = (phi(i+1,j)-two*phi(i,j)+phi(i-1,j))/dx2(1) + &
                     (phi(i,j+1)-two*phi(i,j)+phi(i,j-1))/dy2(1) - &
                     sij(i,j)
      ENDDO
    ENDDO

! Obtain coarse mesh residuals from fine mesh (Restriction)
    DO l = 2,lev
    
    res(:,:,l) = zero !Sets residuals to zero at boundaries
    DO i = 2,n(l)-1
      iff = 2*i-1
      DO j = 2,m(l)-1
        jff = 2*j-1
        res(i,j,l) = res(iff,jff,l-1)
      ENDDO
    ENDDO

! Gauss-Seidel Smoothing for Coarse Mesh (Equation in Correction form)

    diag = two/dx2(l) + two/dy2(l)
    err(:,:,l) = zero
    DO k = 1,3   ! Inner sweeps for course mesh correction
      DO i = 2,n(l)-1
        DO j = 2,m(l)-1
          err(i,j,l) = (res(i,j,l)+(err(i+1,j,l)+err(i-1,j,l))/dx2(l)+ &
                       (err(i,j+1,l)+err(i,j-1,l))/dy2(l))/diag
        ENDDO
      ENDDO
    ENDDO

!   Residual Calculation

    res(:,:,l) = zero
    DO i = 2,n(l)-1
      iff = 2*i-1
      DO j = 2,m(l)-1
        jff = 2*j-1
        res(i,j,l) = res(iff,jff,l-1)+(err(i+1,j,l)-two*err(i,j,l)+err(i-1,j,l))/dx2(l) + &
                     (err(i,j+1,l)-two*err(i,j,l)+err(i,j-1,l))/dy2(l)
      ENDDO
    ENDDO

    ENDDO ! Recursive Restriction

! Obtain fine mesh correction from coarse mesh (Prolongation)

    err(:,:,1) = zero
    DO l = lev,2,-1

    DO i = 2,n(l)-1
      iff = 2*i-1
      DO j = 2,m(l)-1
        jff = 2*j-1

        err(iff,jff,l-1) = err(iff,jff,l-1) + err(i,j,l)
        err(iff-1,jff,l-1) = err(iff-1,jff,l-1) + (err(i,j,l)+err(i-1,j,l))/two
        err(iff,jff-1,l-1) = err(iff,jff-1,l-1) + (err(i,j,l)+err(i,j-1,l))/two
        err(iff-1,jff-1,l-1) = err(iff-1,jff-1,l-1)+ (err(i,j,l)+err(i-1,j,l) + &
                                err(i,j-1,l)+err(i-1,j-1,l))/four
      ENDDO
    ENDDO
     
    ENDDO  ! Recursive Prolongation

! Update fine mesh solution
    phi(:,:) = phi(:,:) + err(:,:,1)

!   Calculate L2NORM on fine mesh
    sumres = zero
    DO j = 2,mff-1
      DO i = 2,nff-1
        fde = res(i,j,1)
        sumres = sumres + fde*fde
      ENDDO
    ENDDO
    res2 = sqrt(MAX(tiny,sumres))

    print *, iter, res2
    WRITE(20,*) iter, res2
    IF(res2 > tol .AND. iter < max_iter) GO TO 11    ! Convergence check

    time_end = etime(ta)

    print *, "CPU time taken =", time_end-time_start

! Print results for fine mesh
    WRITE(10,*) 'VARIABLES = "X", "Y", "PHIAN", "PHI", "Error"'
    WRITE(10,*) 'ZONE I=', nff, ', J=', mff, ', DATAPACKING=POINT'
    DO j = 1,mff
      DO i = 1,nff
        phia(i,j) = ((x(i,1)-aa)**2)*sinh(cc*(x(i,1)-aa)) + &
                    ((y(j,1)-bb)**2)*sinh(cc*(y(j,1)-bb)) + dd*exp(ee*x(i,1)*y(j,1))
        WRITE(10,10) x(i,1),y(j,1),phia(i,j),phi(i,j),phia(i,j)-phi(i,j)
      ENDDO
    ENDDO

 10 format(5(1x,F13.6))

    CLOSE(unit=10)
    CLOSE(unit=20)

    END Program main

