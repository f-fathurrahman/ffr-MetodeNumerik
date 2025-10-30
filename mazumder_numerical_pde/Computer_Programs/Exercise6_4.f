! Start of Main Program
! This is a program to solve the 2D advection-diffusion equation
! in cylidrical coordinates using Stone's method.
! Exercise 6.4

    program main


! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: nc = 200, mc = 40
    INTEGER(4), PARAMETER :: max_iter = 1000
    INTEGER(4) :: i,j,iter,neqn,k
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0, alpha = 0.92d0
    REAL(8), PARAMETER :: length = 0.5D0, height = 0.1d0, &
                           phi_in = 0.0d0, phi_top = 1.0d0
    REAL(8), PARAMETER :: tol = 1.0d-6
    REAL(8), DIMENSION(nc) :: z
    REAL(8), DIMENSION(nc) :: d1,a1,c1,b1,sol1
    REAL(8), DIMENSION(mc) :: r,ff
    REAL(8), DIMENSION(mc) :: d2,a2,c2,b2,sol2
    REAL(8), DIMENSION(nc,mc) :: phi, phi_o
    REAL(8), DIMENSION(:), ALLOCATABLE :: BBB,DDD,EEE,FFF,HHH,QQQ,SSS,YY,DEL,RR
    REAL(8), DIMENSION(:), ALLOCATABLE :: d,e,f,b,c
    REAL(8) :: dz,dr,dzdr,drdz,dr2
    REAL(8) :: res,Pe,res1
    REAL(8) :: flux,sum_flux,flux_l,flux_r,flux_b,flux_t,pi,imb
    

! Open output file
    OPEN(unit=10, file="Exercise6_4.out", status = "unknown")
    OPEN(unit=20, file="Exercise6_4.rsl", status = "unknown")
    OPEN(unit=30, file="Exercise6_4_flux.out", status = "unknown")

    pi = 4.0d0*datan(one)
    Pe = 10.0d0
    neqn = nc*mc

! Mesh
    dz = length/nc; dr = height/mc; dr2 = half*dr
    drdz = dr/dz; dzdr = dz/dr

    z(1) = half*dz
    DO i = 2,nc
      z(i) = z(i-1)+dz
    ENDDO
    r(1) = half*dr
    ff(1) = 2.0d0*(one-(r(1)/height)**2)/height
    DO j = 2,mc
      r(j) = r(j-1)+dr
      ff(j) = 2.0d0*(one-(r(j)/height)**2)/height
    ENDDO


    phi(:,:) = zero  ! Set initial guess 

! Define coefficient matrix [A], and Right hand Vector [Q]

    ALLOCATE(BBB(neqn),DDD(neqn),EEE(neqn),FFF(neqn),HHH(neqn),QQQ(neqn),SSS(neqn))
    ALLOCATE(YY(neqn),DEL(neqn),RR(neqn))
    ALLOCATE(d(neqn),e(neqn),f(neqn),b(neqn),c(neqn))

    EEE(:) = zero; BBB(:) = zero; DDD(:) = zero; FFF(:) = zero; HHH(:) = zero
    QQQ(:) = zero

! Interior Nodes
    DO i = 2, nc-1
      DO j = 2,mc-1
        k = (j-1)*nc + i
        EEE(k) = Pe*ff(j)*r(j)*dr + 2.0d0*drdz*r(j) + dzdr*(r(j)+dr2) + dzdr*(r(j)-dr2)
        BBB(k) = -(r(j)-dr2)*dzdr !South
        HHH(k) = -(r(j)+dr2)*dzdr !North
        DDD(k) = -(Pe*ff(j)*r(j)*dr+r(j)*drdz) !West
        FFF(k) = -r(j)*drdz !East
        QQQ(k) = zero  
        SSS(k) = phi(i,j)  ! Copy to 1D array
      ENDDO
    ENDDO

! Left Boundary (inlet)
    i = 1
    DO j = 2,mc-1
      k = (j-1)*nc + i
      EEE(k) = Pe*ff(j)*r(j)*dr + 4.0d0*drdz*r(j) + dzdr*(r(j)+dr2) + dzdr*(r(j)-dr2)
      BBB(k) = -(r(j)-dr2)*dzdr !South
      HHH(k) = -(r(j)+dr2)*dzdr !North
      DDD(k) = zero !West
      FFF(k) = -4.0d0*r(j)*drdz/3.0d0 !East
      QQQ(k) = (Pe*ff(j)*r(j)*dr+8.0d0*r(j)*drdz/3.0d0)*phi_in  
      SSS(k) = phi(i,j)  ! Copy to 1D array
    ENDDO
! Right Boundary
    i = nc
    DO j = 2,mc-1
      k = (j-1)*nc + i
      EEE(k) = Pe*ff(j)*r(j)*dr + drdz*r(j) + dzdr*(r(j)+dr2) + dzdr*(r(j)-dr2)
      BBB(k) = -(r(j)-dr2)*dzdr !South
      HHH(k) = -(r(j)+dr2)*dzdr !North
      DDD(k) = -(Pe*ff(j)*r(j)*dr+r(j)*drdz) !West
      FFF(k) = zero !East
      QQQ(k) = zero  
      SSS(k) = phi(i,j)  ! Copy to 1D array
    ENDDO
! Bottom Boundary
    j = 1
    DO i = 2,nc-1
      k = (j-1)*nc + i
      EEE(k) = Pe*ff(j)*r(j)*dr + 2.0d0*drdz*r(j) + dzdr*(r(j)+dr2)
      BBB(k) = zero !South
      HHH(k) = -(r(j)+dr2)*dzdr !North
      DDD(k) = -(Pe*ff(j)*r(j)*dr+r(j)*drdz) !West
      FFF(k) = -r(j)*drdz !East
      QQQ(k) = zero  
      SSS(k) = phi(i,j)  ! Copy to 1D array
    ENDDO
! Top Boundary
    j = mc
    DO i = 2,nc-1
      k = (j-1)*nc + i
      EEE(k) = Pe*ff(j)*r(j)*dr + 2.0d0*drdz*r(j) + 3.0d0*dzdr*(r(j)+dr2) + dzdr*(r(j)-dr2)
      BBB(k) = -(r(j)-dr2)*dzdr-dzdr*(r(j)+dr2)/3.0d0 !South
      HHH(k) = zero !North
      DDD(k) = -(Pe*ff(j)*r(j)*dr+r(j)*drdz) !West
      FFF(k) = -r(j)*drdz !East
      QQQ(k) = 8.0d0*(r(j)+dr2)*dzdr*phi_top/3.0d0
      SSS(k) = phi(i,j)  ! Copy to 1D array
    ENDDO
! Inlet bottom corner
    i = 1; j = 1
    k = (j-1)*nc + i
    EEE(k) = Pe*ff(j)*r(j)*dr + 4.0d0*drdz*r(j) + dzdr*(r(j)+dr2)
    BBB(k) = zero !South
    HHH(k) = -(r(j)+dr2)*dzdr !North
    DDD(k) = zero !West
    FFF(k) = -4.0d0*r(j)*drdz/3.0d0 !East
    QQQ(k) = (Pe*ff(j)*r(j)*dr + 8.0d0*drdz*r(j)/3.0d0)*phi_in  
    SSS(k) = phi(i,j)  ! Copy to 1D array
! Inlet top corner
    i = 1; j = mc
    k = (j-1)*nc + i
    EEE(k) = Pe*ff(j)*r(j)*dr + 4.0d0*drdz*r(j) + 3.0d0*dzdr*(r(j)+dr2) + dzdr*(r(j)-dr2)
    BBB(k) = -(r(j)-dr2)*dzdr-dzdr*(r(j)+dr2)/3.0d0 !South
    HHH(k) = zero !North
    DDD(k) = zero !West
    FFF(k) = -4.0d0*r(j)*drdz/3.0d0 !East
    QQQ(k) = (Pe*ff(j)*r(j)*dr + 8.0d0*drdz*r(j)/3.0d0)*phi_in + 8.0d0*(r(j)+dr2)*dzdr*phi_top/3.0d0  
    SSS(k) = phi(i,j)  ! Copy to 1D array
! Outlet bottom corner
    i = nc; j = 1
    k = (j-1)*nc + i
    EEE(k) = Pe*ff(j)*r(j)*dr + drdz*r(j) + dzdr*(r(j)+dr2)
    BBB(k) = zero !South
    HHH(k) = -(r(j)+dr2)*dzdr !North
    DDD(k) = -(Pe*ff(j)*r(j)*dr+r(j)*drdz) !West
    FFF(k) = zero !East
    QQQ(k) = zero  
    SSS(k) = phi(i,j)  ! Copy to 1D array
! Outlet top corner
    i = nc; j = mc
    k = (j-1)*nc + i
    EEE(k) = Pe*ff(j)*r(j)*dr + drdz*r(j) + 3.0d0*dzdr*(r(j)+dr2) + dzdr*(r(j)-dr2)
    BBB(k) = -(r(j)-dr2)*dzdr-dzdr*(r(j)+dr2)/3.0d0 !South
    HHH(k) = zero !North
    DDD(k) = -(Pe*ff(j)*r(j)*dr + drdz*r(j)) !West
    FFF(k) = zero !East
    QQQ(k) = 8.0d0*(r(j)+dr2)*dzdr*phi_top/3.0d0  
    SSS(k) = phi(i,j)  ! Copy to 1D array

! Solve system using Stone's strongly implicit method

! Compute coefficients of L and U
    d(:)=zero;e(:)=zero;c(:)=zero;b(:)=zero;f(:)=zero
    d(1) = EEE(1) 
    c(1) = zero
    b(1) = zero
    e(1) = FFF(1)/d(1)
    f(1) = HHH(1)/d(1)
    DO k = 2,neqn
      IF(k > nc)THEN
        b(k) = BBB(k)/(one+alpha*e(k-nc))
      ENDIF
      c(k) = DDD(k)/(one+alpha*f(k-1))
      IF(k > nc)THEN
        d(k) = EEE(k) + alpha*(b(k)*e(k-nc)+c(k)*f(k-1))-b(k)*f(k-nc)- &
               c(k)*e(k-1)
      ELSE
        d(k) = EEE(k) + alpha*(c(k)*f(k-1))-c(k)*e(k-1)
      ENDIF
      f(k) = (HHH(k)-alpha*c(k)*f(k-1))/d(k)
      IF(k > nc)THEN
        e(k) = (FFF(k) - alpha*b(k)*e(k-nc))/d(k)
      ELSE
        e(k) = FFF(k)/d(k)
      ENDIF
    ENDDO

! Start of Iteration loop
    iter = 0
 90 iter = iter + 1

! Compute Residual Vector
    RR(:) = zero ! Stays zero at boundaries for Dirichlet BC
! Interior Nodes
    DO i = 2, nc-1
      DO j = 2,mc-1
        k = (j-1)*nc + i
        RR(k) = QQQ(k) - EEE(k)*SSS(k) - FFF(k)*SSS(k+1) - HHH(k)*SSS(k+nc) - &
                DDD(k)*SSS(k-1) - BBB(k)*SSS(k-nc)
      ENDDO
    ENDDO
! Left (inlet) Boundary Nodes
    i = 1
    DO j = 2,mc-1
      k = (j-1)*nc + i
      RR(k) = QQQ(k) - EEE(k)*SSS(k) - HHH(k)*SSS(k+nc) - &
                FFF(k)*SSS(k+1) - BBB(k)*SSS(k-nc)
    ENDDO
! Right (outlet) Boundary Nodes
    i = nc
    DO j = 2,mc-1
      k = (j-1)*nc + i
      RR(k) = QQQ(k) - EEE(k)*SSS(k) - HHH(k)*SSS(k+nc) - &
                DDD(k)*SSS(k-1) - BBB(k)*SSS(k-nc)
    ENDDO
! Bottom Boundary Nodes
    j = 1
    DO i = 2,nc-1
      k = (j-1)*nc + i
      RR(k) = QQQ(k) - EEE(k)*SSS(k) - FFF(k)*SSS(k+1) - HHH(k)*SSS(k+nc) - &
                DDD(k)*SSS(k-1)
    ENDDO
! Top Boundary Nodes
    j = mc
    DO i = 2,nc-1
      k = (j-1)*nc + i
      RR(k) = QQQ(k) - EEE(k)*SSS(k) - FFF(k)*SSS(k+1) - &
                DDD(k)*SSS(k-1) - BBB(k)*SSS(k-nc)
    ENDDO
! Inlet bottom corner
   i = 1; j = 1
   k = (j-1)*nc + i
   RR(k) = QQQ(k) - EEE(k)*SSS(k) - FFF(k)*SSS(k+1) - HHH(k)*SSS(k+nc)
! Inlet top corner
   i = 1; j = mc
   k = (j-1)*nc + i
   RR(k) = QQQ(k) - EEE(k)*SSS(k) - FFF(k)*SSS(k+1) - BBB(k)*SSS(k-nc)
! Outlet bottom corner
   i = nc; j = 1
   k = (j-1)*nc + i
   RR(k) = QQQ(k) - EEE(k)*SSS(k) - HHH(k)*SSS(k+nc) - DDD(k)*SSS(k-1) 
! Outlet top corner
   i = nc; j = mc
   k = (j-1)*nc + i
   RR(k) = QQQ(k) - EEE(k)*SSS(k) - DDD(k)*SSS(k-1) - BBB(k)*SSS(k-nc)


! Forward substitution: L Y = R
    YY(1) = RR(1)/d(1)
    DO k = 2,neqn
      IF(k > nc)THEN
        YY(k) = (RR(k) - c(k)*YY(k-1) - b(k)*YY(k-nc))/d(k)
      ELSE
        YY(k) = (RR(k) - c(k)*YY(k-1))/d(k)
      ENDIF
    ENDDO

! Backward Substitution: U Delta = Y
    DEL(:) = zero
    DEL(neqn) = YY(neqn)
    DO k = neqn-1,1,-1
      IF(k > neqn-nc)THEN
        DEL(k) = YY(k) - e(k)*DEL(k+1)
      ELSE
        DEL(k) = YY(k) - e(k)*DEL(k+1) - f(k)*DEL(k+nc)
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
    IF(res > tol .and. iter < max_iter) GO TO 90

! Transform to 2D array
    DO i = 1,nc
      DO j = 1,mc
        k = (j-1)*nc + i
        phi(i,j) = SSS(k)
      ENDDO
    ENDDO

    
! Print results

! Contour Data
    WRITE(10,*) 'VARIABLES = "Z", "R", "PHI"'
    WRITE(10,*) 'ZONE I=', nc+2, ', J=', mc+2, ', DATAPACKING=POINT'

!Bottom Row
    WRITE(10,10) zero,zero,phi_in
    DO i = 1,nc
      WRITE(10,10) z(i),zero,(9.0d0*phi(i,1)-phi(i,2))/8.0d0
    ENDDO
    WRITE(10,10) length,zero,(9.0d0*phi(nc,1)-phi(nc,2))/8.0d0

!Interior Rows
    DO j = 1,mc
      WRITE(10,10) zero, r(j), phi_in
      DO i = 1,nc
         WRITE(10,10) z(i),r(j),phi(i,j)
      ENDDO
      WRITE(10,10) length, r(j), phi(nc,j)
    ENDDO

!Top Row
    WRITE(10,10) zero,height,phi_top
    DO i = 1,nc
      WRITE(10,10) z(i),height,phi_top
    ENDDO
    WRITE(10,10) length,height,phi_top

 10 format(3(1x,F13.6))

! Fluxes
! Inlet
  sum_flux = zero
  WRITE(30,*) "Inlet"
  DO j = 1,mc
    flux = Pe*ff(j)*phi_in - (9.0d0*phi(1,j)-8.0d0*phi_in-phi(2,j))/(3.0d0*dz)
    sum_flux = sum_flux + flux*dr*r(j)
    WRITE(30,*) r(j), ABS(flux)
  ENDDO
  flux_l = sum_flux*two*pi
! Outlet
  sum_flux = zero
  WRITE(30,*) "Outlet"
  DO j = 1,mc
    flux = Pe*ff(j)*phi(nc,j)
    sum_flux = sum_flux + flux*dr*r(j)
    WRITE(30,*) r(j), ABS(flux)
  ENDDO
  flux_r = sum_flux*two*pi
! Top
  sum_flux = zero
  WRITE(30,*) "Top"
  DO i = 1,nc
    flux = +(9.0d0*phi(i,mc)-8.0d0*phi_top-phi(i,mc-1))/(3.0d0*dr)
    sum_flux = sum_flux + flux*(height)*dz
    WRITE(30,*) z(i), ABS(flux)
  ENDDO
  flux_t = sum_flux*two*pi
  flux_b = zero
  imb = flux_r-flux_l+flux_t-flux_b
  print *, flux_r,flux_l,flux_t,flux_b,imb
  WRITE(30,*) "Imbalance"
  WRITE(30,*) flux_r,flux_l,flux_t,flux_b,imb
  

    CLOSE(unit=10)
    CLOSE(unit=20)

    END Program main

